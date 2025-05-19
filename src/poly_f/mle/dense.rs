use crate::field_config::ConfigPtr;
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
use crate::field::RandomField;
use crate::sparse_matrix::SparseMatrix;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseMultilinearExtension<'cfg, const N: usize> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: Vec<RandomField<'cfg, N>>,
    /// Number of variables
    pub num_vars: usize,
    /// Field in which the MLE is operating
    pub config: ConfigPtr<'cfg, N>,
}

impl<'cfg, const N: usize> DenseMultilinearExtension<'cfg, N> {
    pub type Field = RandomField<'cfg, N>;
    pub type Cfg = ConfigPtr<'cfg, N>;

    pub fn from_evaluations_slice(
        num_vars: usize,
        evaluations: &[Self::Field],
        config: Self::Cfg,
    ) -> Self {
        Self::from_evaluations_vec(num_vars, evaluations.to_vec(), config)
    }

    pub fn evaluate(&self, point: &[Self::Field], config: Self::Cfg) -> Option<Self::Field> {
        if point.len() == self.num_vars {
            Some(self.fixed_variables(point, config)[0])
        } else {
            None
        }
    }

    pub fn from_evaluations_vec(
        num_vars: usize,
        evaluations: Vec<Self::Field>,
        config: Self::Cfg,
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
                Self::Field::new_unchecked(config, 0u32.into()),
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
    pub fn from_matrix(matrix: &SparseMatrix<Self::Field>, config: Self::Cfg) -> Self {
        let n_vars: usize = (log2(matrix.nrows()) + log2(matrix.ncols())) as usize; // n_vars = s + s'

        // Matrices might need to get padded before turned into an MLE
        let padded_rows = matrix.n_rows.next_power_of_two();
        let padded_cols = matrix.n_cols.next_power_of_two();

        // build dense vector representing the sparse padded matrix
        let mut v = vec![Self::Field::zero(); padded_rows * padded_cols];

        for (row_i, row) in matrix.coeffs.iter().enumerate() {
            for (val, col_i) in row {
                v[(padded_rows * *col_i) + row_i] = *val;
            }
        }

        // convert the dense vector into a mle
        Self::from_slice(n_vars, &v, config)
    }

    /// Takes n_vars and a dense slice and returns its dense MLE.
    pub fn from_slice(n_vars: usize, v: &[Self::Field], config: Self::Cfg) -> Self {
        let v_padded: Vec<Self::Field> = if v.len() != (1 << n_vars) {
            // pad to 2^n_vars
            [
                v.to_owned(),
                ark_std::iter::repeat(Self::Field::zero())
                    .take((1 << n_vars) - v.len())
                    .collect(),
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

impl<'cfg, const N: usize> MultilinearExtension<'cfg, N> for DenseMultilinearExtension<'cfg, N> {
    type Field = RandomField<'cfg, N>;
    type Cfg = ConfigPtr<'cfg, N>;
    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn rand<Rn: rand::Rng>(num_vars: usize, config: Self::Cfg, rng: &mut Rn) -> Self {
        Self::from_evaluations_vec(
            num_vars,
            (0..1 << num_vars).map(|_| Self::Field::rand(rng)).collect(),
            config,
        )
    }

    fn relabel(&self, a: usize, b: usize, k: usize) -> Self {
        let mut copy = self.clone();
        copy.relabel_in_place(a, b, k);
        copy
    }

    fn fix_variables(&mut self, partial_point: &[Self::Field], _config: Self::Cfg) {
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

    fn fixed_variables(&self, partial_point: &[Self::Field], config: Self::Cfg) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point, config);
        res
    }

    fn to_evaluations(&self) -> Vec<Self::Field> {
        self.evaluations.to_vec()
    }
}

impl<const N: usize> Zero for DenseMultilinearExtension<'_, N> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: vec![RandomField::<N>::zero()],
            config: ConfigPtr::NONE,
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations[0].is_zero()
    }
}

impl<const N: usize> Add for DenseMultilinearExtension<'_, N> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<'cfg, const N: usize> Add for &DenseMultilinearExtension<'cfg, N> {
    type Output = DenseMultilinearExtension<'cfg, N>;

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
            .map(|(a, b)| *a + *b)
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result, self.config)
    }
}

impl<const N: usize> AddAssign for DenseMultilinearExtension<'_, N> {
    fn add_assign(&mut self, rhs: Self) {
        self.add_assign(&rhs);
    }
}

impl<const N: usize> AddAssign<&Self> for DenseMultilinearExtension<'_, N> {
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

impl<'cfg, const N: usize> AddAssign<(RandomField<'cfg, N>, &Self)>
    for DenseMultilinearExtension<'cfg, N>
{
    fn add_assign(&mut self, (r, other): (RandomField<'cfg, N>, &Self)) {
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

impl<const N: usize> Neg for DenseMultilinearExtension<'_, N> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        cfg_iter_mut!(self.evaluations).for_each(|a| *a = a.neg());

        self
    }
}

impl<const N: usize> Sub for DenseMultilinearExtension<'_, N> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        &self - &other
    }
}

impl<'cfg, const N: usize> Sub for &DenseMultilinearExtension<'cfg, N> {
    type Output = DenseMultilinearExtension<'cfg, N>;

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
            .map(|(a, b)| *a - *b)
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result, self.config)
    }
}

impl<const N: usize> SubAssign for DenseMultilinearExtension<'_, N> {
    fn sub_assign(&mut self, other: Self) {
        self.sub_assign(&other);
    }
}

impl<const N: usize> SubAssign<&Self> for DenseMultilinearExtension<'_, N> {
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
            .for_each(|(a, b)| a.sub_assign(b));
    }
}

impl<'cfg, const N: usize> Mul<RandomField<'cfg, N>> for DenseMultilinearExtension<'cfg, N> {
    type Output = Self;

    fn mul(mut self, rhs: RandomField<'cfg, N>) -> Self::Output {
        self.evaluations.iter_mut().for_each(|x| *x *= rhs);

        self
    }
}

impl<'cfg, const N: usize> MulAssign<RandomField<'cfg, N>> for DenseMultilinearExtension<'cfg, N> {
    fn mul_assign(&mut self, rhs: RandomField<'cfg, N>) {
        self.evaluations.iter_mut().for_each(|x| *x *= rhs);
    }
}

impl<'cfg, const N: usize> Sub<RandomField<'cfg, N>> for DenseMultilinearExtension<'cfg, N> {
    type Output = Self;

    fn sub(mut self, rhs: RandomField<'cfg, N>) -> Self::Output {
        self.evaluations.iter_mut().for_each(|x| *x -= rhs);

        self
    }
}

impl<'cfg, const N: usize> Add<RandomField<'cfg, N>> for DenseMultilinearExtension<'cfg, N> {
    type Output = Self;

    fn add(mut self, rhs: RandomField<'cfg, N>) -> Self::Output {
        self.evaluations.iter_mut().for_each(|x| *x += &rhs);

        self
    }
}

impl<'cfg, const N: usize> Index<usize> for DenseMultilinearExtension<'cfg, N> {
    type Output = RandomField<'cfg, N>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.evaluations[index]
    }
}

impl<const N: usize> IndexMut<usize> for DenseMultilinearExtension<'_, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.evaluations[index]
    }
}

unsafe impl<const N: usize> Send for DenseMultilinearExtension<'_, N> {}

unsafe impl<const N: usize> Sync for DenseMultilinearExtension<'_, N> {}
