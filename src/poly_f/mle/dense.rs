use core::ops::IndexMut;

use ark_ff::Zero;
use ark_std::{
    borrow::ToOwned,
    cfg_iter, cfg_iter_mut, log2,
    ops::{Add, AddAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign},
    rand, vec,
    vec::Vec,
};
use crypto_bigint::Random;
#[cfg(feature = "parallel")]
use rayon::iter::*;

use super::{MultilinearExtension, swap_bits};
use crate::{field::RandomField, sparse_matrix::SparseMatrix, traits::ConfigReference};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseMultilinearExtension<C: ConfigReference> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: Vec<RandomField<C>>,
    /// Number of variables
    pub num_vars: usize,
    /// Field in which the MLE is operating
    pub config: C,
}

impl<C: ConfigReference> DenseMultilinearExtension<C> {
    pub fn from_evaluations_slice(
        num_vars: usize,
        evaluations: &[RandomField<C>],
        config: C,
    ) -> Self {
        Self::from_evaluations_vec(num_vars, evaluations.to_vec(), config)
    }

    pub fn evaluate(&self, point: &[RandomField<C>], config: C) -> Option<RandomField<C>> {
        if point.len() == self.num_vars {
            Some(self.fixed_variables(point, config)[0].clone())
        } else {
            None
        }
    }

    pub fn from_evaluations_vec(
        num_vars: usize,
        evaluations: Vec<RandomField<C>>,
        config: C,
    ) -> Self {
        // assert that the number of variables matches the size of evaluations
        assert!(
            evaluations.len() <= 1 << num_vars,
            "The size of evaluations should not exceed 2^num_vars. \n eval len: {:?}. num vars: {num_vars}",
            evaluations.len()
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
    pub fn from_matrix(matrix: &SparseMatrix<RandomField<C>>, config: C) -> Self {
        let n_vars: usize = (log2(matrix.nrows()) + log2(matrix.ncols())) as usize; // n_vars = s + s'

        // Matrices might need to get padded before turned into an MLE
        let padded_rows = matrix.n_rows.next_power_of_two();
        let padded_cols = matrix.n_cols.next_power_of_two();

        // build dense vector representing the sparse padded matrix
        let mut v = vec![RandomField::zero(); padded_rows * padded_cols];

        for (row_i, row) in matrix.coeffs.iter().enumerate() {
            for (val, col_i) in row {
                v[(padded_rows * *col_i) + row_i] = val.clone();
            }
        }

        // convert the dense vector into a mle
        Self::from_slice(n_vars, &v, config)
    }

    /// Takes n_vars and a dense slice and returns its dense MLE.
    pub fn from_slice(n_vars: usize, v: &[RandomField<C>], config: C) -> Self {
        let v_padded: Vec<RandomField<C>> = if v.len() != (1 << n_vars) {
            // pad to 2^n_vars
            [
                v.to_owned(),
                ark_std::iter::repeat_n(RandomField::zero(), (1 << n_vars) - v.len()).collect(),
            ]
            .concat()
        } else {
            v.to_owned()
        };
        Self::from_evaluations_vec(n_vars, v_padded, config)
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

impl<C: ConfigReference> MultilinearExtension<C> for DenseMultilinearExtension<C> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn rand<Rn: rand::Rng>(num_vars: usize, config: C, rng: &mut Rn) -> Self {
        Self::from_evaluations_vec(
            num_vars,
            (0..1 << num_vars)
                .map(|_| RandomField::random(rng))
                .collect(),
            config,
        )
    }

    fn relabel(&self, a: usize, b: usize, k: usize) -> Self {
        let mut copy = self.clone();
        copy.relabel_in_place(a, b, k);
        copy
    }

    fn fix_variables(&mut self, partial_point: &[RandomField<C>], _config: C) {
        assert!(
            partial_point.len() <= self.num_vars,
            "too many partial points"
        );

        let poly = &mut self.evaluations;
        let nv = self.num_vars;
        let dim = partial_point.len();

        for i in 1..dim + 1 {
            let r = partial_point[i - 1].clone();
            for b in 0..1 << (nv - i) {
                let left = poly[b << 1].clone();
                let right = poly[(b << 1) + 1].clone();
                let a = right - left.clone();
                if !a.is_zero() {
                    poly[b] = left + r.clone() * a;
                } else {
                    poly[b] = left;
                };
            }
        }

        self.evaluations.truncate(1 << (nv - dim));
        self.num_vars = nv - dim;
    }

    fn fixed_variables(&self, partial_point: &[RandomField<C>], config: C) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point, config);
        res
    }

    fn to_evaluations(&self) -> Vec<RandomField<C>> {
        self.evaluations.to_vec()
    }
}

impl<C: ConfigReference> Zero for DenseMultilinearExtension<C> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: vec![RandomField::zero()],
            config: C::NONE,
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations[0].is_zero()
    }
}

impl<C: ConfigReference> Add for DenseMultilinearExtension<C> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<C: ConfigReference> Add for &DenseMultilinearExtension<C> {
    type Output = DenseMultilinearExtension<C>;

    fn add(self, rhs: Self) -> Self::Output {
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
            .map(|(a, b)| a.clone() + b.clone())
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result, self.config)
    }
}

impl<C: ConfigReference> AddAssign for DenseMultilinearExtension<C> {
    fn add_assign(&mut self, rhs: Self) {
        self.add_assign(&rhs);
    }
}

impl<C: ConfigReference> AddAssign<&Self> for DenseMultilinearExtension<C> {
    fn add_assign(&mut self, other: &Self) {
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

impl<C: ConfigReference> AddAssign<(RandomField<C>, &Self)> for DenseMultilinearExtension<C> {
    fn add_assign(&mut self, (r, other): (RandomField<C>, &Self)) {
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
            .for_each(|(a, b)| a.add_assign(&(r.clone() * b)));
    }
}

impl<C: ConfigReference> Neg for DenseMultilinearExtension<C> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        cfg_iter_mut!(self.evaluations).for_each(|a| *a = a.clone().neg());

        self
    }
}

impl<C: ConfigReference> Sub for DenseMultilinearExtension<C> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        &self - &other
    }
}

impl<C: ConfigReference> Sub for &DenseMultilinearExtension<C> {
    type Output = DenseMultilinearExtension<C>;

    fn sub(self, rhs: Self) -> Self::Output {
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
            .map(|(a, b)| a.clone() - b.clone())
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result, self.config)
    }
}

impl<C: ConfigReference> SubAssign for DenseMultilinearExtension<C> {
    fn sub_assign(&mut self, other: Self) {
        self.sub_assign(&other);
    }
}

impl<C: ConfigReference> SubAssign<&Self> for DenseMultilinearExtension<C> {
    fn sub_assign(&mut self, rhs: &Self) {
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
            .for_each(|(a, b)| a.sub_assign(b.clone()));
    }
}

impl<C: ConfigReference> Mul<RandomField<C>> for DenseMultilinearExtension<C> {
    type Output = Self;

    fn mul(mut self, rhs: RandomField<C>) -> Self::Output {
        self.evaluations.iter_mut().for_each(|x| *x *= &rhs);

        self
    }
}

impl<C: ConfigReference> MulAssign<RandomField<C>> for DenseMultilinearExtension<C> {
    fn mul_assign(&mut self, rhs: RandomField<C>) {
        self.evaluations.iter_mut().for_each(|x| *x *= &rhs);
    }
}

impl<C: ConfigReference> Sub<RandomField<C>> for DenseMultilinearExtension<C> {
    type Output = Self;

    fn sub(mut self, rhs: RandomField<C>) -> Self::Output {
        self.evaluations.iter_mut().for_each(|x| *x -= rhs.clone());

        self
    }
}

impl<C: ConfigReference> Add<RandomField<C>> for DenseMultilinearExtension<C> {
    type Output = Self;

    fn add(mut self, rhs: RandomField<C>) -> Self::Output {
        self.evaluations.iter_mut().for_each(|x| *x += &rhs);

        self
    }
}

impl<C: ConfigReference> Index<usize> for DenseMultilinearExtension<C> {
    type Output = RandomField<C>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.evaluations[index]
    }
}

impl<C: ConfigReference> IndexMut<usize> for DenseMultilinearExtension<C> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.evaluations[index]
    }
}

unsafe impl<C: ConfigReference> Send for DenseMultilinearExtension<C> {}

unsafe impl<C: ConfigReference> Sync for DenseMultilinearExtension<C> {}
