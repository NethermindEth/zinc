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
use crate::sparse_matrix::SparseMatrix;
use crate::{field::RandomField, field_config::FieldConfig};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseMultilinearExtension<const N: usize> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: Vec<RandomField<N>>,
    /// Number of variables
    pub num_vars: usize,
    /// Field in which the MLE is operating
    pub config: *const FieldConfig<N>,
}

impl<const N: usize> DenseMultilinearExtension<N> {
    pub fn from_evaluations_slice(
        num_vars: usize,
        evaluations: &[RandomField<N>],
        config: *const FieldConfig<N>,
    ) -> Self {
        Self::from_evaluations_vec(num_vars, evaluations.to_vec(), config)
    }

    pub fn evaluate(
        &self,
        point: &[RandomField<N>],
        config: *const FieldConfig<N>,
    ) -> Option<RandomField<N>> {
        if point.len() == self.num_vars {
            Some(self.fixed_variables(point, config)[0])
        } else {
            None
        }
    }

    pub fn from_evaluations_vec(
        num_vars: usize,
        evaluations: Vec<RandomField<N>>,
        config: *const FieldConfig<N>,
    ) -> Self {
        // assert that the number of variables matches the size of evaluations
        assert!(
            evaluations.len() <= 1 << num_vars,
            "The size of evaluations should not exceed 2^num_vars. \n eval len: {:?}. num vars: {num_vars}", evaluations.len()
        );

        if evaluations.len() != 1 << num_vars {
            let mut evaluations = evaluations;
            evaluations.resize(
                1 << num_vars,
                RandomField::new_unchecked(config, 0u32.into()),
            );
            return Self {
                num_vars,
                evaluations,
                config,
            };
        }

        Self {
            num_vars,
            evaluations,
            config,
        }
    }

    /// Returns the dense MLE from the given matrix, without modifying the original matrix.
    pub fn from_matrix(
        matrix: &SparseMatrix<RandomField<N>>,
        config: *const FieldConfig<N>,
    ) -> Self {
        let n_vars: usize = (log2(matrix.nrows()) + log2(matrix.ncols())) as usize; // n_vars = s + s'

        // Matrices might need to get padded before turned into an MLE
        let padded_rows = matrix.n_rows.next_power_of_two();
        let padded_cols = matrix.n_cols.next_power_of_two();

        // build dense vector representing the sparse padded matrix
        let mut v = vec![RandomField::<N>::zero(); padded_rows * padded_cols];

        for (row_i, row) in matrix.coeffs.iter().enumerate() {
            for (val, col_i) in row {
                v[(padded_rows * *col_i) + row_i] = *val;
            }
        }

        // convert the dense vector into a mle
        Self::from_slice(n_vars, &v, config)
    }

    /// Takes n_vars and a dense slice and returns its dense MLE.
    pub fn from_slice(n_vars: usize, v: &[RandomField<N>], config: *const FieldConfig<N>) -> Self {
        let v_padded: Vec<RandomField<N>> = if v.len() != (1 << n_vars) {
            // pad to 2^n_vars
            [
                v.to_owned(),
                ark_std::iter::repeat(RandomField::<N>::zero())
                    .take((1 << n_vars) - v.len())
                    .collect(),
            ]
            .concat()
        } else {
            v.to_owned()
        };
        DenseMultilinearExtension::<N>::from_evaluations_vec(n_vars, v_padded, config)
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

impl<const N: usize> MultilinearExtension<N> for DenseMultilinearExtension<N> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn rand<Rn: rand::Rng>(num_vars: usize, config: *const FieldConfig<N>, rng: &mut Rn) -> Self {
        Self::from_evaluations_vec(
            num_vars,
            (0..1 << num_vars)
                .map(|_| RandomField::<N>::rand(rng))
                .collect(),
            config,
        )
    }

    fn relabel(&self, a: usize, b: usize, k: usize) -> Self {
        let mut copy = self.clone();
        copy.relabel_in_place(a, b, k);
        copy
    }

    fn fix_variables(&mut self, partial_point: &[RandomField<N>], _config: *const FieldConfig<N>) {
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

    fn fixed_variables(
        &self,
        partial_point: &[RandomField<N>],
        config: *const FieldConfig<N>,
    ) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point, config);
        res
    }

    fn to_evaluations(&self) -> Vec<RandomField<N>> {
        self.evaluations.to_vec()
    }
}

impl<const N: usize> Zero for DenseMultilinearExtension<N> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: vec![RandomField::<N>::zero()],
            config: std::ptr::null(),
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations[0].is_zero()
    }
}

impl<const N: usize> Add for DenseMultilinearExtension<N> {
    type Output = DenseMultilinearExtension<N>;

    fn add(self, other: DenseMultilinearExtension<N>) -> Self {
        &self + &other
    }
}

impl<'a, const N: usize> Add<&'a DenseMultilinearExtension<N>> for &DenseMultilinearExtension<N> {
    type Output = DenseMultilinearExtension<N>;

    fn add(self, rhs: &'a DenseMultilinearExtension<N>) -> Self::Output {
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
        assert_eq!(
            self.config, rhs.config,
            "trying to add two dense MLEs in different fields"
        );

        let result = cfg_iter!(self.evaluations)
            .zip(cfg_iter!(rhs.evaluations))
            .map(|(a, b)| *a + *b)
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result, self.config)
    }
}

impl<const N: usize> AddAssign for DenseMultilinearExtension<N> {
    fn add_assign(&mut self, rhs: DenseMultilinearExtension<N>) {
        self.add_assign(&rhs);
    }
}

impl<'a, const N: usize> AddAssign<&'a DenseMultilinearExtension<N>>
    for DenseMultilinearExtension<N>
{
    fn add_assign(&mut self, other: &'a DenseMultilinearExtension<N>) {
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
        assert_eq!(
            self.config, other.config,
            "trying to add two dense MLEs in different fields"
        );

        cfg_iter_mut!(self.evaluations)
            .zip(cfg_iter!(other.evaluations))
            .for_each(|(a, b)| a.add_assign(b));
    }
}

impl<const N: usize> AddAssign<(RandomField<N>, &DenseMultilinearExtension<N>)>
    for DenseMultilinearExtension<N>
{
    fn add_assign(&mut self, (r, other): (RandomField<N>, &DenseMultilinearExtension<N>)) {
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

        assert_eq!(
            self.config, other.config,
            "trying to add two dense MLEs in different fields"
        );

        cfg_iter_mut!(self.evaluations)
            .zip(cfg_iter!(other.evaluations))
            .for_each(|(a, b)| a.add_assign(&(r * b)));
    }
}

impl<const N: usize> Neg for DenseMultilinearExtension<N> {
    type Output = DenseMultilinearExtension<N>;

    fn neg(mut self) -> Self {
        cfg_iter_mut!(self.evaluations).for_each(|a| *a = a.neg());

        self
    }
}

impl<const N: usize> Sub for DenseMultilinearExtension<N> {
    type Output = DenseMultilinearExtension<N>;

    fn sub(self, other: DenseMultilinearExtension<N>) -> Self {
        &self - &other
    }
}

impl<'a, const N: usize> Sub<&'a DenseMultilinearExtension<N>> for &DenseMultilinearExtension<N> {
    type Output = DenseMultilinearExtension<N>;

    fn sub(self, rhs: &'a DenseMultilinearExtension<N>) -> Self::Output {
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
        assert_eq!(
            self.config, rhs.config,
            "trying to add two dense MLEs in different fields"
        );
        let result = cfg_iter!(self.evaluations)
            .zip(cfg_iter!(rhs.evaluations))
            .map(|(a, b)| *a - *b)
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result, self.config)
    }
}

impl<const N: usize> SubAssign for DenseMultilinearExtension<N> {
    fn sub_assign(&mut self, other: DenseMultilinearExtension<N>) {
        self.sub_assign(&other);
    }
}

impl<'a, const N: usize> SubAssign<&'a DenseMultilinearExtension<N>>
    for DenseMultilinearExtension<N>
{
    fn sub_assign(&mut self, rhs: &'a DenseMultilinearExtension<N>) {
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

impl<const N: usize> Mul<RandomField<N>> for DenseMultilinearExtension<N> {
    type Output = DenseMultilinearExtension<N>;

    fn mul(mut self, rhs: RandomField<N>) -> DenseMultilinearExtension<N> {
        self.evaluations.iter_mut().for_each(|x| *x *= rhs);

        self
    }
}

impl<const N: usize> MulAssign<RandomField<N>> for DenseMultilinearExtension<N> {
    fn mul_assign(&mut self, rhs: RandomField<N>) {
        self.evaluations.iter_mut().for_each(|x| *x *= rhs);
    }
}

impl<const N: usize> Sub<RandomField<N>> for DenseMultilinearExtension<N> {
    type Output = DenseMultilinearExtension<N>;

    fn sub(mut self, rhs: RandomField<N>) -> DenseMultilinearExtension<N> {
        self.evaluations.iter_mut().for_each(|x| *x -= rhs);

        self
    }
}

impl<const N: usize> Add<RandomField<N>> for DenseMultilinearExtension<N> {
    type Output = DenseMultilinearExtension<N>;

    fn add(mut self, rhs: RandomField<N>) -> DenseMultilinearExtension<N> {
        self.evaluations.iter_mut().for_each(|x| *x += &rhs);

        self
    }
}

impl<const N: usize> Index<usize> for DenseMultilinearExtension<N> {
    type Output = RandomField<N>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.evaluations[index]
    }
}

impl<const N: usize> IndexMut<usize> for DenseMultilinearExtension<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.evaluations[index]
    }
}

unsafe impl<const N: usize> Send for DenseMultilinearExtension<N> {}

unsafe impl<const N: usize> Sync for DenseMultilinearExtension<N> {}
