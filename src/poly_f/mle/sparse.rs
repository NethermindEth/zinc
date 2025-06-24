use ark_ff::{UniformRand, Zero};
use ark_std::{
    cfg_iter,
    collections::BTreeMap,
    log2,
    ops::{Add, AddAssign, Index, Neg, Sub, SubAssign},
    rand::Rng,
    vec,
    vec::Vec,
};
use hashbrown::HashMap;
#[cfg(feature = "parallel")]
use rayon::iter::*;

use super::{swap_bits, MultilinearExtension};
use crate::{
    sparse_matrix::SparseMatrix,
    traits::{ConfigReference, Field, FieldMap, Integer},
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseMultilinearExtension<F: Field> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: BTreeMap<usize, F>,
    /// Number of variables
    pub num_vars: usize,
    zero: F,
    /// Field in which the MLE is operating
    pub config: F::Cr,
}
impl<F: Field> SparseMultilinearExtension<F> {
    pub fn from_evaluations<'a>(
        num_vars: usize,
        evaluations: impl IntoIterator<Item = &'a (usize, F)>,
        config: F::Cr,
    ) -> Self
    where
        F: 'a,
    {
        let bit_mask = 1 << num_vars;

        let evaluations: Vec<_> = evaluations
            .into_iter()
            .map(|(i, v): &(usize, F)| {
                assert!(*i < bit_mask, "index out of range");
                (*i, *v)
            })
            .collect();
        Self {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars,
            zero: F::zero(),
            config,
        }
    }
    pub fn evaluate(&self, point: &[F], config: F::Cr) -> F {
        assert!(point.len() == self.num_vars);
        self.fixed_variables(point, config)[0]
    }
    /// Outputs an `l`-variate multilinear extension where value of evaluations
    /// are sampled uniformly at random. The number of nonzero entries is
    /// `num_nonzero_entries` and indices of those nonzero entries are
    /// distributed uniformly at random.
    ///
    /// Note that this function uses rejection sampling. As number of nonzero
    /// entries approach `2 ^ num_vars`, sampling will be very slow due to
    /// large number of collisions.
    pub fn rand_with_config<Rn: Rng>(
        num_vars: usize,
        num_nonzero_entries: usize,
        config: F::Cr,
        rng: &mut Rn,
    ) -> Self {
        assert!(num_nonzero_entries <= 1 << num_vars);

        let mut map = HashMap::new();
        for _ in 0..num_nonzero_entries {
            let mut index = usize::rand(rng) & ((1 << num_vars) - 1);
            while map.contains_key(&index) {
                index = usize::rand(rng) & ((1 << num_vars) - 1);
            }
            map.entry(index).or_insert(F::rand(rng));
        }
        let mut buf = Vec::new();
        for (arg, v) in map.iter() {
            if *v != F::zero() {
                buf.push((*arg, *v));
            }
        }
        let evaluations = hashmap_to_treemap(&map);
        Self {
            num_vars,
            evaluations,
            zero: F::zero(),
            config,
        }
    }

    /// Returns the sparse MLE from the given matrix, without modifying the original matrix.
    pub fn from_matrix(m: &SparseMatrix<F>, config: F::Cr) -> Self {
        let n_rows = m.n_rows.next_power_of_two();
        let n_cols = m.n_cols.next_power_of_two();
        let n_vars: usize = (log2(n_rows * n_cols)) as usize; // n_vars = s + s'

        // build the sparse vec representing the sparse matrix
        let total_elements: usize = m.coeffs.iter().map(|row| row.len()).sum();
        let mut v: Vec<(usize, F)> = Vec::with_capacity(total_elements);

        for (row_i, row) in m.coeffs.iter().enumerate() {
            for (val, col_i) in row {
                let index = row_i * n_cols + col_i;
                v.push((index, *val));
            }
        }

        // convert the sparse vector into a mle
        Self::from_sparse_slice(n_vars, &v, config)
    }

    /// Takes n_vars and a sparse slice and returns its sparse MLE.
    pub fn from_sparse_slice(n_vars: usize, v: &[(usize, F)], config: F::Cr) -> Self {
        SparseMultilinearExtension::<F>::from_evaluations(n_vars, v, config)
    }

    /// Takes n_vars and a dense slice and returns its sparse MLE.
    pub fn from_slice(n_vars: usize, v: &[F], config: F::Cr) -> Self {
        let v_sparse = v
            .iter()
            .enumerate()
            .map(|(i, v_i)| (i, *v_i))
            .collect::<Vec<(usize, F)>>();
        SparseMultilinearExtension::<F>::from_evaluations(n_vars, &v_sparse, config)
    }
}

impl<F: Field> MultilinearExtension<F> for SparseMultilinearExtension<F> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }
    /// Outputs an `l`-variate multilinear extension where value of evaluations
    /// are sampled uniformly at random. The number of nonzero entries is
    /// `sqrt(2^num_vars)` and indices of those nonzero entries are distributed
    /// uniformly at random.
    fn rand<Rn: ark_std::rand::Rng>(num_vars: usize, config: F::Cr, rng: &mut Rn) -> Self {
        Self::rand_with_config(num_vars, 1 << (num_vars / 2), config, rng)
    }

    fn relabel(&self, mut a: usize, mut b: usize, k: usize) -> Self {
        if a > b {
            // swap
            core::mem::swap(&mut a, &mut b);
        }
        // sanity check
        assert!(
            a + k < self.num_vars && b + k < self.num_vars,
            "invalid relabel argument"
        );
        if a == b || k == 0 {
            return self.clone();
        }
        assert!(a + k <= b, "overlapped swap window is not allowed");
        let ev: Vec<_> = cfg_iter!(self.evaluations)
            .map(|(&i, &v)| (swap_bits(i, a, b, k), v))
            .collect();
        Self {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: F::zero(),
            config: self.config,
        }
    }

    fn fix_variables(&mut self, partial_point: &[F], config: F::Cr) {
        let dim = partial_point.len();
        assert!(dim <= self.num_vars, "invalid partial point dimension");

        let mut window = ark_std::log2(self.evaluations.len()) as usize;
        if window == 0 {
            window = 1;
        }
        let mut point = partial_point;
        let mut last = treemap_to_hashmap(&self.evaluations);

        // batch evaluation
        while !point.is_empty() {
            let focus_length = if point.len() > window {
                window
            } else {
                point.len()
            };
            let focus = &point[..focus_length];
            point = &point[focus_length..];
            let pre = precompute_eq(focus, config);
            let dim = focus.len();
            let mut result = HashMap::new();
            for src_entry in last.iter() {
                let old_idx = *src_entry.0;
                let gz = pre[old_idx & ((1 << dim) - 1)];
                let new_idx = old_idx >> dim;
                let dst_entry = result.entry(new_idx).or_insert(F::zero());
                *dst_entry += gz * src_entry.1;
            }
            last = result;
        }
        let evaluations = hashmap_to_treemap(&last);

        self.evaluations = evaluations;
        self.num_vars -= dim;
        self.zero = F::zero();
    }

    fn fixed_variables(&self, partial_point: &[F], config: F::Cr) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point, config);
        res
    }

    fn to_evaluations(&self) -> Vec<F> {
        let mut evaluations: Vec<_> = (0..1 << self.num_vars).map(|_| F::zero()).collect();
        self.evaluations
            .iter()
            .map(|(&i, &v)| {
                evaluations[i] = v;
            })
            .last();
        evaluations
    }
}
impl<F: Field> Zero for SparseMultilinearExtension<F> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: tuples_to_treemap(&Vec::new()),
            zero: F::zero(),
            config: F::Cr::NONE,
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations.is_empty()
    }
}

impl<F: Field> Add for SparseMultilinearExtension<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        &self + &other
    }
}

impl<F: Field> Add for &SparseMultilinearExtension<F> {
    type Output = SparseMultilinearExtension<F>;

    fn add(self, rhs: Self) -> Self::Output {
        // handle zero case
        if self.is_zero() {
            return rhs.clone();
        }
        if rhs.is_zero() {
            return self.clone();
        }

        assert_eq!(
            rhs.num_vars, self.num_vars,
            "trying to add non-zero polynomial with different number of variables"
        );

        assert_eq!(
            rhs.config, self.config,
            "trying to add two sparse MLEs in different fields"
        );
        // simply merge the evaluations
        let mut evaluations = HashMap::new();
        for (&i, &v) in self.evaluations.iter().chain(rhs.evaluations.iter()) {
            *evaluations.entry(i).or_insert(F::zero()) += &v;
        }
        let evaluations: Vec<_> = evaluations
            .into_iter()
            .filter(|(_, v)| !v.is_zero())
            .collect();

        Self::Output {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars: self.num_vars,
            zero: F::zero(),
            config: self.config,
        }
    }
}

impl<F: Field> AddAssign for SparseMultilinearExtension<F> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + &other;
    }
}
impl<F: Field> AddAssign<&Self> for SparseMultilinearExtension<F> {
    fn add_assign(&mut self, rhs: &Self) {
        *self = &*self + rhs;
    }
}
impl<F: Field> AddAssign<(F, &Self)> for SparseMultilinearExtension<F> {
    fn add_assign(&mut self, (r, other): (F, &Self)) {
        if !self.is_zero() && !other.is_zero() {
            assert_eq!(
                other.num_vars, self.num_vars,
                "trying to add non-zero polynomial with different number of variables"
            );
            assert_eq!(
                other.config, self.config,
                "trying to add two MLEs in different fields"
            );
        }
        let ev: Vec<_> = cfg_iter!(other.evaluations)
            .map(|(i, v)| (*i, r * v))
            .collect();
        let other = Self {
            num_vars: other.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: F::zero(),
            config: self.config,
        };
        *self += &other;
    }
}
impl<F: Field> Neg for SparseMultilinearExtension<F> {
    type Output = Self;

    fn neg(self) -> Self {
        let ev: Vec<_> = cfg_iter!(self.evaluations)
            .map(|(i, v)| (*i, -*v))
            .collect();
        Self::Output {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: F::zero(),
            config: self.config,
        }
    }
}
impl<F: Field> Sub for SparseMultilinearExtension<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}
impl<F: Field> Sub for &SparseMultilinearExtension<F> {
    type Output = SparseMultilinearExtension<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: Self) -> Self::Output {
        self + &rhs.clone().neg()
    }
}
impl<F: Field> SubAssign for SparseMultilinearExtension<F> {
    fn sub_assign(&mut self, other: Self) {
        *self = &*self - &other;
    }
}
impl<F: Field> SubAssign<&Self> for SparseMultilinearExtension<F> {
    fn sub_assign(&mut self, rhs: &Self) {
        *self = &*self - rhs;
    }
}

impl<F: Field> Index<usize> for SparseMultilinearExtension<F> {
    type Output = F;

    /// Returns the evaluation of the polynomial at a point represented by
    /// index.
    ///
    /// Index represents a vector in {0,1}^`num_vars` in little endian form. For
    /// example, `0b1011` represents `P(1,1,0,1)`
    ///
    /// For Sparse multilinear polynomial, Lookup_evaluation takes log time to
    /// the size of polynomial.
    fn index(&self, index: usize) -> &Self::Output {
        if let Some(v) = self.evaluations.get(&index) {
            v
        } else {
            &self.zero
        }
    }
}

/// Utilities
fn tuples_to_treemap<F: Field>(tuples: &[(usize, F)]) -> BTreeMap<usize, F> {
    BTreeMap::from_iter(tuples.iter().map(|(i, v)| (*i, *v)))
}

fn treemap_to_hashmap<F: Field>(map: &BTreeMap<usize, F>) -> HashMap<usize, F> {
    HashMap::from_iter(map.iter().map(|(i, v)| (*i, *v)))
}
fn hashmap_to_treemap<F: Field>(map: &HashMap<usize, F>) -> BTreeMap<usize, F> {
    BTreeMap::from_iter(map.iter().map(|(i, v)| (*i, *v)))
}

// precompute  f(x) = eq(g,x)
fn precompute_eq<F: Field>(g: &[F], config: F::Cr) -> Vec<F> {
    let dim = g.len();
    let mut dp = vec![F::zero(); 1 << dim];
    dp[0] = <F::I as FieldMap<F>>::map_to_field(&F::I::one(), config) - g[0];
    dp[1] = g[0];
    for i in 1..dim {
        for b in 0..1 << i {
            let prev = dp[b];
            dp[b + (1 << i)] = prev * g[i];
            dp[b] = prev - dp[b + (1 << i)];
        }
    }
    dp
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;
    use crate::{
        field::RandomField,
        field_config::{ConfigRef, FieldConfig},
    };

    // Function to convert usize to a binary vector of Ring elements.
    fn usize_to_binary_vector<F: Field>(n: usize, dimensions: usize, config: F::Cr) -> Vec<F> {
        let mut bits = Vec::with_capacity(dimensions);
        let mut current = n;

        for _ in 0..dimensions {
            if (current & 1) == 1 {
                bits.push(1u64.map_to_field(config));
            } else {
                bits.push(0u64.map_to_field(config));
            }
            current >>= 1;
        }
        bits
    }

    // Wrapper function to generate a boolean hypercube.
    fn boolean_hypercube<F: Field>(dimensions: usize, config: F::Cr) -> Vec<Vec<F>> {
        let max_val = 1 << dimensions; // 2^dimensions
        (0..max_val)
            .map(|i| usize_to_binary_vector(i, dimensions, config))
            .collect()
    }

    fn vec_cast<F: Field>(v: &[usize], config: F::Cr) -> Vec<F> {
        v.iter().map(|c| (*c as u64).map_to_field(config)).collect()
    }

    fn matrix_cast<F: Field>(m: &[Vec<usize>], config: F::Cr) -> SparseMatrix<F> {
        let n_rows = m.len();
        let n_cols = m[0].len();
        let mut coeffs = Vec::with_capacity(n_rows);
        for row in m.iter() {
            let mut row_coeffs = Vec::with_capacity(n_cols);
            for (col_i, &val) in row.iter().enumerate() {
                if val != 0 {
                    row_coeffs.push(((val as u64).map_to_field(config), col_i));
                }
            }
            coeffs.push(row_coeffs);
        }
        SparseMatrix {
            n_rows,
            n_cols,
            coeffs,
        }
    }

    fn get_test_z<F: Field>(input: usize, config: F::Cr) -> Vec<F> {
        vec_cast(
            &[
                input, // io
                1,
                input * input * input + input + 5, // x^3 + x + 5
                input * input,                     // x^2
                input * input * input,             // x^2 * x
                input * input * input + input,     // x^3 + x
            ],
            config,
        )
    }

    #[test]
    fn test_matrix_to_mle() {
        const N: usize = 1;
        let config = FieldConfig::new(293u32.into());
        let config_ptr: ConfigRef<1> = ConfigRef::from(&config);
        let A = matrix_cast::<RandomField<N>>(
            &[
                vec![2, 3, 4, 4],
                vec![4, 11, 14, 14],
                vec![2, 8, 17, 17],
                vec![420, 4, 2, 0],
            ],
            config_ptr,
        );

        let A_mle = SparseMultilinearExtension::from_matrix(&A, config_ptr);
        assert_eq!(A_mle.evaluations.len(), 15); // 15 non-zero elements
        assert_eq!(A_mle.num_vars, 4); // 4x4 matrix, thus 2bit x 2bit, thus 2^4=16 evals

        let A = matrix_cast::<RandomField<N>>(
            &[
                vec![2, 3, 4, 4, 1],
                vec![4, 11, 14, 14, 2],
                vec![2, 8, 17, 17, 3],
                vec![420, 4, 2, 0, 4],
                vec![420, 4, 2, 0, 5],
            ],
            config_ptr,
        );
        let A_mle = SparseMultilinearExtension::from_matrix(&A, config_ptr);
        assert_eq!(A_mle.evaluations.len(), 23); // 23 non-zero elements
        assert_eq!(A_mle.num_vars, 6); // 5x5 matrix, thus 3bit x 3bit, thus 2^6=64 evals
    }

    #[test]
    fn test_vec_to_mle() {
        const N: usize = 1;
        let config = FieldConfig::new(293u32.into());
        let config = ConfigRef::from(&config);
        let z: Vec<RandomField<N>> = get_test_z(3, config);

        let n_vars = 3;
        let z_mle = SparseMultilinearExtension::from_slice(n_vars, &z, config);

        // check that the z_mle evaluated over the boolean hypercube equals the vec z_i values
        let bhc = boolean_hypercube(z_mle.num_vars, config);
        for (i, z_i) in z.iter().enumerate() {
            let s_i = &bhc[i];
            assert_eq!(z_mle.evaluate(s_i, config), z_i.clone());
        }
        // for the rest of elements of the boolean hypercube, expect it to evaluate to zero
        for s_i in bhc.iter().take(1 << z_mle.num_vars).skip(z.len()) {
            assert_eq!(z_mle.fixed_variables(s_i, config)[0], RandomField::zero());
        }
    }
}
