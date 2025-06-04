use crate::field_config::ConfigRef;
use core::ops::IndexMut;

use ark_ff::Zero;

use ark_std::{
    borrow::ToOwned,
    cfg_iter, cfg_iter_mut, log2,
    ops::{Add, AddAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign},
    rand, vec,
    vec::Vec,
};
use crypto_bigint::{Int, Random};
#[cfg(feature = "parallel")]
use rayon::iter::*;

use super::{swap_bits, MultilinearExtension};
use crate::poly::ArithErrors;
use crate::{
    field::conversion::FieldMap,
    poly_f::mle::DenseMultilinearExtension as DenseMultilinearExtensionF,
    sparse_matrix::SparseMatrix,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseMultilinearExtension<const I: usize> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: Vec<Int<I>>,
    /// Number of variables
    pub num_vars: usize,
}

impl<const I: usize> DenseMultilinearExtension<I> {
    pub fn from_evaluations_slice(num_vars: usize, evaluations: &[Int<I>]) -> Self {
        Self::from_evaluations_vec(num_vars, evaluations.to_vec())
    }

    pub fn evaluate(&self, point: &[Int<I>]) -> Option<Int<I>> {
        if point.len() == self.num_vars {
            Some(self.fixed_variables(point)[0])
        } else {
            None
        }
    }

    pub fn from_evaluations_vec(num_vars: usize, evaluations: Vec<Int<I>>) -> Self {
        // assert that the number of variables matches the size of evaluations
        assert!(
            evaluations.len() <= 1 << num_vars,
            "The size of evaluations should not exceed 2^num_vars. \n eval len: {:?}. num vars: {num_vars}", evaluations.len()
        );

        if evaluations.len() != 1 << num_vars {
            let mut evaluations = evaluations;
            evaluations.resize(1 << num_vars, Int::<I>::zero());
            return Self {
                num_vars,
                evaluations,
            };
        }

        Self {
            num_vars,
            evaluations,
        }
    }

    /// Returns the dense MLE from the given matrix, without modifying the original matrix.
    pub fn from_matrix(matrix: &SparseMatrix<Int<I>>) -> Self {
        let n_vars: usize = (log2(matrix.nrows()) + log2(matrix.ncols())) as usize; // n_vars = s + s'

        // Matrices might need to get padded before turned into an MLE
        let padded_rows = matrix.n_rows.next_power_of_two();
        let padded_cols = matrix.n_cols.next_power_of_two();

        // build dense vector representing the sparse padded matrix
        let mut v = vec![Int::<I>::ZERO; padded_rows * padded_cols];

        for (row_i, row) in matrix.coeffs.iter().enumerate() {
            for (val, col_i) in row {
                v[(padded_cols * row_i) + *col_i] = *val;
            }
        }

        // convert the dense vector into a mle
        Self::from_slice(n_vars, &v)
    }

    /// Takes n_vars and a dense slice and returns its dense MLE.
    pub fn from_slice(n_vars: usize, v: &[Int<I>]) -> Self {
        let v_padded: Vec<Int<I>> = if v.len() != (1 << n_vars) {
            // pad to 2^n_vars
            [
                v.to_owned(),
                ark_std::iter::repeat(Int::<I>::zero())
                    .take((1 << n_vars) - v.len())
                    .collect(),
            ]
            .concat()
        } else {
            v.to_owned()
        };
        DenseMultilinearExtension::from_evaluations_vec(n_vars, v_padded)
    }

    pub fn relabel_in_place(&mut self, mut a: usize, mut b: usize, k: usize) {
        // enforce order of a and b
        if a > b {
            ark_std::mem::swap(&mut a, &mut b);
        }
        if a == b || k == 0 {
            return;
        }
        assert!(b + k <= self.num_vars, "invalid relabel argument");
        assert!(a + k <= b, "overlapped swap window is not allowed");
        for i in 0..self.evaluations.len() {
            let j = swap_bits(i, a, b, k);
            if i < j {
                self.evaluations.swap(i, j);
            }
        }
    }
}

impl<'cfg, const I: usize, const N: usize> FieldMap<'cfg, N> for DenseMultilinearExtension<I> {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = DenseMultilinearExtensionF<'cfg, N>;
    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        DenseMultilinearExtensionF::from_evaluations_vec(
            self.num_vars,
            self.evaluations.map_to_field(config),
            config,
        )
    }
}

impl<const I: usize> MultilinearExtension<I> for DenseMultilinearExtension<I> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn rand<Rn: rand::Rng>(num_vars: usize, rng: &mut Rn) -> Self {
        Self::from_evaluations_vec(
            num_vars,
            (0..1 << num_vars).map(|_| Int::<I>::random(rng)).collect(),
        )
    }

    fn relabel(&self, a: usize, b: usize, k: usize) -> Self {
        let mut copy = self.clone();
        copy.relabel_in_place(a, b, k);
        copy
    }

    fn fix_variables(&mut self, partial_point: &[Int<I>]) {
        assert!(
            partial_point.len() <= self.num_vars,
            "too many partial points"
        );

        let poly = &mut self.evaluations;
        let nv = self.num_vars;
        let dim = partial_point.len();

        for i in 1..dim + 1 {
            let r = partial_point[i - 1];
            for b in 0..1 << (nv - i) {
                let left = poly[b << 1];
                let right = poly[(b << 1) + 1];
                let a = right - left;
                if !a.is_zero() {
                    poly[b] = left + r * a;
                } else {
                    poly[b] = left;
                };
            }
        }

        self.evaluations.truncate(1 << (nv - dim));
        self.num_vars = nv - dim;
    }

    fn fixed_variables(&self, partial_point: &[Int<I>]) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point);
        res
    }

    fn to_evaluations(&self) -> Vec<Int<I>> {
        self.evaluations.to_vec()
    }
}

impl<const I: usize> Zero for DenseMultilinearExtension<I> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: vec![Int::<I>::zero()],
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations[0].is_zero()
    }
}

impl<const I: usize> Add for DenseMultilinearExtension<I> {
    type Output = DenseMultilinearExtension<I>;

    fn add(self, other: DenseMultilinearExtension<I>) -> Self {
        &self + &other
    }
}

impl<'a, const I: usize> Add<&'a DenseMultilinearExtension<I>> for &DenseMultilinearExtension<I> {
    type Output = DenseMultilinearExtension<I>;

    fn add(self, rhs: &'a DenseMultilinearExtension<I>) -> Self::Output {
        if rhs.is_zero() {
            return self.clone();
        }

        if self.is_zero() {
            return rhs.clone();
        }

        assert_eq!(
            self.num_vars, rhs.num_vars,
            "trying to add two dense MLEs with different numbers of variables"
        );

        let result = cfg_iter!(self.evaluations)
            .zip(cfg_iter!(rhs.evaluations))
            .map(|(a, b)| *a + *b)
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result)
    }
}

impl<const I: usize> AddAssign for DenseMultilinearExtension<I> {
    fn add_assign(&mut self, rhs: DenseMultilinearExtension<I>) {
        self.add_assign(&rhs);
    }
}

impl<'a, const I: usize> AddAssign<&'a DenseMultilinearExtension<I>>
    for DenseMultilinearExtension<I>
{
    fn add_assign(&mut self, other: &'a DenseMultilinearExtension<I>) {
        if self.is_zero() {
            *self = other.clone();
            return;
        }

        if other.is_zero() {
            return;
        }

        assert_eq!(
            self.num_vars, other.num_vars,
            "trying to add two dense MLEs with different numbers of variables"
        );

        cfg_iter_mut!(self.evaluations)
            .zip(cfg_iter!(other.evaluations))
            .for_each(|(a, b)| a.add_assign(b));
    }
}

impl<const I: usize> AddAssign<(Int<I>, &DenseMultilinearExtension<I>)>
    for DenseMultilinearExtension<I>
{
    fn add_assign(&mut self, (r, other): (Int<I>, &DenseMultilinearExtension<I>)) {
        if self.is_zero() {
            *self = other.clone();

            cfg_iter_mut!(self.evaluations).for_each(|a| *a = *a * (r));

            return;
        }

        if other.is_zero() {
            return;
        }

        assert_eq!(
            self.num_vars, other.num_vars,
            "trying to add two dense MLEs with different numbers of variables"
        );

        cfg_iter_mut!(self.evaluations)
            .zip(cfg_iter!(other.evaluations))
            .for_each(|(a, b)| a.add_assign(&(r * b)));
    }
}

impl<const I: usize> Neg for DenseMultilinearExtension<I> {
    type Output = DenseMultilinearExtension<I>;

    fn neg(mut self) -> Self {
        cfg_iter_mut!(self.evaluations).for_each(|a| *a = Int::<I>::zero() - *a);

        self
    }
}

impl<const I: usize> Sub for DenseMultilinearExtension<I> {
    type Output = Self;

    fn sub(self, other: DenseMultilinearExtension<I>) -> Self {
        &self - &other
    }
}

impl<'a, const I: usize> Sub<&'a DenseMultilinearExtension<I>> for &DenseMultilinearExtension<I> {
    type Output = DenseMultilinearExtension<I>;

    fn sub(self, rhs: &'a DenseMultilinearExtension<I>) -> Self::Output {
        if rhs.is_zero() {
            return self.clone();
        }

        if self.is_zero() {
            return rhs.clone().neg();
        }

        assert_eq!(
            self.num_vars, rhs.num_vars,
            "trying to subtract two dense MLEs with different numbers of variables"
        );

        let result = cfg_iter!(self.evaluations)
            .zip(cfg_iter!(rhs.evaluations))
            .map(|(a, b)| *a - *b)
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result)
    }
}

impl<const I: usize> SubAssign for DenseMultilinearExtension<I> {
    fn sub_assign(&mut self, other: DenseMultilinearExtension<I>) {
        self.sub_assign(&other);
    }
}

impl<'a, const I: usize> SubAssign<&'a DenseMultilinearExtension<I>>
    for DenseMultilinearExtension<I>
{
    fn sub_assign(&mut self, rhs: &'a DenseMultilinearExtension<I>) {
        if self.is_zero() {
            *self = rhs.clone().neg();
            return;
        }

        if rhs.is_zero() {
            return;
        }

        assert_eq!(
            self.num_vars, rhs.num_vars,
            "trying to subtract two dense MLEs with different numbers of variables"
        );

        cfg_iter_mut!(self.evaluations)
            .zip(cfg_iter!(rhs.evaluations))
            .for_each(|(a, b)| *a = *a - *b);
    }
}

impl<const I: usize> Mul<Int<I>> for DenseMultilinearExtension<I> {
    type Output = DenseMultilinearExtension<I>;

    fn mul(mut self, rhs: Int<I>) -> DenseMultilinearExtension<I> {
        self.evaluations.iter_mut().for_each(|x| *x = *x * rhs);

        self
    }
}

impl<const I: usize> MulAssign<Int<I>> for DenseMultilinearExtension<I> {
    fn mul_assign(&mut self, rhs: Int<I>) {
        self.evaluations.iter_mut().for_each(|x| *x = *x * rhs);
    }
}

impl<const I: usize> Sub<Int<I>> for DenseMultilinearExtension<I> {
    type Output = DenseMultilinearExtension<I>;

    fn sub(mut self, rhs: Int<I>) -> DenseMultilinearExtension<I> {
        self.evaluations.iter_mut().for_each(|x| *x = *x - rhs);

        self
    }
}

impl<const I: usize> Add<Int<I>> for DenseMultilinearExtension<I> {
    type Output = DenseMultilinearExtension<I>;

    fn add(mut self, rhs: Int<I>) -> DenseMultilinearExtension<I> {
        self.evaluations.iter_mut().for_each(|x| *x += &rhs);

        self
    }
}

impl<const I: usize> Index<usize> for DenseMultilinearExtension<I> {
    type Output = Int<I>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.evaluations[index]
    }
}

impl<const I: usize> IndexMut<usize> for DenseMultilinearExtension<I> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.evaluations[index]
    }
}

unsafe impl<const I: usize> Send for DenseMultilinearExtension<I> {}

unsafe impl<const I: usize> Sync for DenseMultilinearExtension<I> {}

/// This function build the eq(x, r) polynomial for any given r.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r<const I: usize>(
    r: &[Int<I>],
) -> Result<DenseMultilinearExtension<I>, ArithErrors> {
    let evals = build_eq_x_r_vec(r)?;
    let mle = DenseMultilinearExtension::from_evaluations_vec(r.len(), evals);

    Ok(mle)
}

/// This function build the eq(x, r) polynomial for any given r, and output the
/// evaluation of eq(x, r) in its vector form.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r_vec<const I: usize>(r: &[Int<I>]) -> Result<Vec<Int<I>>, ArithErrors> {
    // we build eq(x,r) from its evaluations
    // we want to evaluate eq(x,r) over x \in {0, 1}^num_vars
    // for example, with num_vars = 4, x is a binary vector of 4, then
    //  0 0 0 0 -> (1-r0)   * (1-r1)    * (1-r2)    * (1-r3)
    //  1 0 0 0 -> r0       * (1-r1)    * (1-r2)    * (1-r3)
    //  0 1 0 0 -> (1-r0)   * r1        * (1-r2)    * (1-r3)
    //  1 1 0 0 -> r0       * r1        * (1-r2)    * (1-r3)
    //  ....
    //  1 1 1 1 -> r0       * r1        * r2        * r3
    // we will need 2^num_var evaluations

    let mut eval = Vec::new();
    build_eq_x_r_helper(r, &mut eval)?;

    Ok(eval)
}

/// A helper function to build eq(x, r) recursively.
/// This function takes `r.len()` steps, and for each step it requires a maximum
/// `r.len()-1` multiplications.
fn build_eq_x_r_helper<const I: usize>(
    r: &[Int<I>],
    buf: &mut Vec<Int<I>>,
) -> Result<(), ArithErrors> {
    if r.is_empty() {
        return Err(ArithErrors::InvalidParameters("r length is 0".into()));
    } else if r.len() == 1 {
        // initializing the buffer with [1-r_0, r_0]
        buf.push(Int::<I>::ONE - r[0]);
        buf.push(r[0]);
    } else {
        build_eq_x_r_helper(&r[1..], buf)?;

        // suppose at the previous step we received [b_1, ..., b_k]
        // for the current step we will need
        // if x_0 = 0:   (1-r0) * [b_1, ..., b_k]
        // if x_0 = 1:   r0 * [b_1, ..., b_k]
        // let mut res = vec![];
        // for &b_i in buf.iter() {
        //     let tmp = r[0] * b_i;
        //     res.push(b_i - tmp);
        //     res.push(tmp);
        // }
        // *buf = res;

        let mut res = vec![Int::<I>::zero(); buf.len() << 1];
        cfg_iter_mut!(res).enumerate().for_each(|(i, val)| {
            let bi = buf[i >> 1];
            let tmp = r[0] * bi;
            if (i & 1) == 0 {
                *val = bi - tmp;
            } else {
                *val = tmp;
            }
        });
        *buf = res;
    }

    Ok(())
}
