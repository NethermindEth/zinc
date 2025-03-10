use core::ops::IndexMut;

use ark_ff::{UniformRand, Zero};

use ark_std::{
    borrow::ToOwned,
    cfg_iter, cfg_iter_mut, log2,
    ops::{Add, AddAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign},
    rand,
    vec::*,
};
#[cfg(feature = "parallel")]
use rayon::iter::*;

use super::{swap_bits, MultilinearExtension};
use crate::{
    field::conversion::FieldMap, field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension as DenseMultilinearExtensionF,
    poly_z::polynomials::ArithErrors, sparse_matrix::SparseMatrix,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseMultilinearExtension {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: Vec<i64>,
    /// Number of variables
    pub num_vars: usize,
}

impl DenseMultilinearExtension {
    pub fn from_evaluations_slice(num_vars: usize, evaluations: &[i64]) -> Self {
        Self::from_evaluations_vec(num_vars, evaluations.to_vec())
    }

    pub fn evaluate(&self, point: &[i64]) -> Option<i64> {
        if point.len() == self.num_vars {
            Some(self.fixed_variables(point)[0])
        } else {
            None
        }
    }

    pub fn to_random_field<const N: usize>(
        &self,
        config: *const FieldConfig<N>,
    ) -> DenseMultilinearExtensionF<N> {
        let evaluations = self
            .evaluations
            .iter()
            .map(|x| x.map_to_field(config))
            .collect();
        DenseMultilinearExtensionF::from_evaluations_vec(self.num_vars, evaluations, config)
    }

    pub fn from_evaluations_vec(num_vars: usize, evaluations: Vec<i64>) -> Self {
        // assert that the number of variables matches the size of evaluations
        assert!(
            evaluations.len() <= 1 << num_vars,
            "The size of evaluations should not exceed 2^num_vars. \n eval len: {:?}. num vars: {num_vars}", evaluations.len()
        );

        if evaluations.len() != 1 << num_vars {
            let mut evaluations = evaluations;
            evaluations.resize(1 << num_vars, 0i64);
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
    pub fn from_matrix(matrix: &SparseMatrix<i64>) -> Self {
        let n_vars: usize = (log2(matrix.nrows()) + log2(matrix.ncols())) as usize; // n_vars = s + s'

        // Matrices might need to get padded before turned into an MLE
        let padded_rows = matrix.n_rows.next_power_of_two();
        let padded_cols = matrix.n_cols.next_power_of_two();

        // build dense vector representing the sparse padded matrix
        let mut v = vec![0i64; padded_rows * padded_cols];

        for (row_i, row) in matrix.coeffs.iter().enumerate() {
            for (val, col_i) in row {
                v[(padded_cols * row_i) + *col_i] = *val;
            }
        }

        // convert the dense vector into a mle
        Self::from_slice(n_vars, &v)
    }

    /// Takes n_vars and a dense slice and returns its dense MLE.
    pub fn from_slice(n_vars: usize, v: &[i64]) -> Self {
        let v_padded: Vec<i64> = if v.len() != (1 << n_vars) {
            // pad to 2^n_vars
            [
                v.to_owned(),
                ark_std::iter::repeat(0i64)
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

impl MultilinearExtension for DenseMultilinearExtension {
    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn rand<Rn: rand::Rng>(num_vars: usize, rng: &mut Rn) -> Self {
        Self::from_evaluations_vec(
            num_vars,
            (0..1 << num_vars).map(|_| i64::rand(rng)).collect(),
        )
    }

    fn relabel(&self, a: usize, b: usize, k: usize) -> Self {
        let mut copy = self.clone();
        copy.relabel_in_place(a, b, k);
        copy
    }

    fn fix_variables(&mut self, partial_point: &[i64]) {
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

    fn fixed_variables(&self, partial_point: &[i64]) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point);
        res
    }

    fn to_evaluations(&self) -> Vec<i64> {
        self.evaluations.to_vec()
    }
}

impl Zero for DenseMultilinearExtension {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: vec![0i64],
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations[0].is_zero()
    }
}

impl Add for DenseMultilinearExtension {
    type Output = DenseMultilinearExtension;

    fn add(self, other: DenseMultilinearExtension) -> Self {
        &self + &other
    }
}

impl<'a> Add<&'a DenseMultilinearExtension> for &DenseMultilinearExtension {
    type Output = DenseMultilinearExtension;

    fn add(self, rhs: &'a DenseMultilinearExtension) -> Self::Output {
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

impl AddAssign for DenseMultilinearExtension {
    fn add_assign(&mut self, rhs: DenseMultilinearExtension) {
        self.add_assign(&rhs);
    }
}

impl<'a> AddAssign<&'a DenseMultilinearExtension> for DenseMultilinearExtension {
    fn add_assign(&mut self, other: &'a DenseMultilinearExtension) {
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

impl AddAssign<(i64, &DenseMultilinearExtension)> for DenseMultilinearExtension {
    fn add_assign(&mut self, (r, other): (i64, &DenseMultilinearExtension)) {
        if self.is_zero() {
            *self = other.clone();

            cfg_iter_mut!(self.evaluations).for_each(|a| a.mul_assign(&r));

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

impl Neg for DenseMultilinearExtension {
    type Output = DenseMultilinearExtension;

    fn neg(mut self) -> Self {
        cfg_iter_mut!(self.evaluations).for_each(|a| *a = a.neg());

        self
    }
}

impl Sub for DenseMultilinearExtension {
    type Output = Self;

    fn sub(self, other: DenseMultilinearExtension) -> Self {
        &self - &other
    }
}

impl<'a> Sub<&'a DenseMultilinearExtension> for &DenseMultilinearExtension {
    type Output = DenseMultilinearExtension;

    fn sub(self, rhs: &'a DenseMultilinearExtension) -> Self::Output {
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

impl SubAssign for DenseMultilinearExtension {
    fn sub_assign(&mut self, other: DenseMultilinearExtension) {
        self.sub_assign(&other);
    }
}

impl<'a> SubAssign<&'a DenseMultilinearExtension> for DenseMultilinearExtension {
    fn sub_assign(&mut self, rhs: &'a DenseMultilinearExtension) {
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
            .for_each(|(a, b)| a.sub_assign(b));
    }
}

impl Mul<i64> for DenseMultilinearExtension {
    type Output = DenseMultilinearExtension;

    fn mul(mut self, rhs: i64) -> DenseMultilinearExtension {
        self.evaluations.iter_mut().for_each(|x| *x *= rhs);

        self
    }
}

impl MulAssign<i64> for DenseMultilinearExtension {
    fn mul_assign(&mut self, rhs: i64) {
        self.evaluations.iter_mut().for_each(|x| *x *= rhs);
    }
}

impl Sub<i64> for DenseMultilinearExtension {
    type Output = DenseMultilinearExtension;

    fn sub(mut self, rhs: i64) -> DenseMultilinearExtension {
        self.evaluations.iter_mut().for_each(|x| *x -= rhs);

        self
    }
}

impl Add<i64> for DenseMultilinearExtension {
    type Output = DenseMultilinearExtension;

    fn add(mut self, rhs: i64) -> DenseMultilinearExtension {
        self.evaluations.iter_mut().for_each(|x| *x += &rhs);

        self
    }
}

impl Index<usize> for DenseMultilinearExtension {
    type Output = i64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.evaluations[index]
    }
}

impl IndexMut<usize> for DenseMultilinearExtension {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.evaluations[index]
    }
}

unsafe impl Send for DenseMultilinearExtension {}

unsafe impl Sync for DenseMultilinearExtension {}

/// This function build the eq(x, r) polynomial for any given r.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r(r: &[i64]) -> Result<DenseMultilinearExtension, ArithErrors> {
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
pub fn build_eq_x_r_vec(r: &[i64]) -> Result<Vec<i64>, ArithErrors> {
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
fn build_eq_x_r_helper(r: &[i64], buf: &mut Vec<i64>) -> Result<(), ArithErrors> {
    if r.is_empty() {
        return Err(ArithErrors::InvalidParameters("r length is 0".to_string()));
    } else if r.len() == 1 {
        // initializing the buffer with [1-r_0, r_0]
        buf.push(1i64 - r[0]);
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

        let mut res = vec![0i64; buf.len() << 1];
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
