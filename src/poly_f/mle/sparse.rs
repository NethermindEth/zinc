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

use super::{MultilinearExtension, swap_bits};
use crate::{
    field::RandomField,
    sparse_matrix::SparseMatrix,
    traits::{BigInteger, ConfigReference, FieldMap, InSameField},
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseMultilinearExtension<C: ConfigReference> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: BTreeMap<usize, RandomField<C>>,
    /// Number of variables
    pub num_vars: usize,
    zero: RandomField<C>,
}
impl<C: ConfigReference> SparseMultilinearExtension<C> {
    pub fn from_evaluations<'a>(
        num_vars: usize,
        evaluations: impl IntoIterator<Item = &'a (usize, RandomField<C>)>,
    ) -> Self
    where
        RandomField<C>: 'a,
    {
        let bit_mask = 1 << num_vars;

        let evaluations: Vec<_> = evaluations
            .into_iter()
            .map(|(i, v): &(usize, RandomField<C>)| {
                assert!(*i < bit_mask, "index out of range");
                (*i, v.clone())
            })
            .collect();
        let config = match evaluations.first() {
            None => None,
            Some(v) => v.1.config(),
        };

        Self {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars,
            zero: match config {
                None => RandomField::zero(),
                Some(config) => 0u64.map_to_field(config),
            },
        }
    }
    pub fn evaluate(&self, point: &[RandomField<C>], config: C) -> RandomField<C> {
        assert!(point.len() == self.num_vars);
        self.fixed_variables(point, config)[0].clone()
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
        config: C,
        rng: &mut Rn,
    ) -> Self {
        assert!(num_nonzero_entries <= 1 << num_vars);

        let mut map = HashMap::new();
        for _ in 0..num_nonzero_entries {
            let mut index = usize::rand(rng) & ((1 << num_vars) - 1);
            while map.contains_key(&index) {
                index = usize::rand(rng) & ((1 << num_vars) - 1);
            }
            map.entry(index)
                .or_insert(RandomField::rand_with_config(rng, config));
        }

        let zero = 0u64.map_to_field(config);

        let mut buf = Vec::new();
        for (arg, v) in map.iter() {
            if v != &zero {
                buf.push((*arg, v.clone()));
            }
        }
        let evaluations = hashmap_to_treemap(&map);
        Self {
            num_vars,
            evaluations,
            zero: RandomField::zero(),
        }
    }

    /// Returns the sparse MLE from the given matrix, without modifying the original matrix.
    pub fn from_matrix(m: &SparseMatrix<RandomField<C>>) -> Self {
        let n_rows = m.n_rows.next_power_of_two();
        let n_cols = m.n_cols.next_power_of_two();
        let n_vars: usize = (log2(n_rows * n_cols)) as usize; // n_vars = s + s'

        // build the sparse vec representing the sparse matrix
        let total_elements: usize = m.coeffs.iter().map(|row| row.len()).sum();
        let mut v: Vec<(usize, RandomField<C>)> = Vec::with_capacity(total_elements);

        for (row_i, row) in m.coeffs.iter().enumerate() {
            for (val, col_i) in row {
                let index = row_i * n_cols + col_i;
                v.push((index, val.clone()));
            }
        }

        // convert the sparse vector into a mle
        Self::from_sparse_slice(n_vars, &v)
    }

    /// Takes n_vars and a sparse slice and returns its sparse MLE.
    pub fn from_sparse_slice(n_vars: usize, v: &[(usize, RandomField<C>)]) -> Self {
        SparseMultilinearExtension::from_evaluations(n_vars, v)
    }

    /// Takes n_vars and a dense slice and returns its sparse MLE.
    pub fn from_slice(n_vars: usize, v: &[RandomField<C>]) -> Self {
        let v_sparse = v
            .iter()
            .enumerate()
            .map(|(i, v_i)| (i, v_i.clone()))
            .collect::<Vec<(usize, RandomField<C>)>>();
        SparseMultilinearExtension::from_evaluations(n_vars, &v_sparse)
    }
}

impl<C: ConfigReference> MultilinearExtension<C> for SparseMultilinearExtension<C> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }
    /// Outputs an `l`-variate multilinear extension where value of evaluations
    /// are sampled uniformly at random. The number of nonzero entries is
    /// `sqrt(2^num_vars)` and indices of those nonzero entries are distributed
    /// uniformly at random.
    fn rand<Rn: ark_std::rand::Rng>(num_vars: usize, config: C, rng: &mut Rn) -> Self {
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
            .map(|(&i, v)| (swap_bits(i, a, b, k), v.clone()))
            .collect();
        Self {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: RandomField::zero(),
        }
    }

    fn fix_variables(&mut self, partial_point: &[RandomField<C>], config: C) {
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
                let gz = pre[old_idx & ((1 << dim) - 1)].clone();
                let new_idx = old_idx >> dim;
                let dst_entry = result.entry(new_idx).or_insert(RandomField::zero());
                *dst_entry += gz * src_entry.1;
            }
            last = result;
        }
        let evaluations = hashmap_to_treemap(&last);

        self.evaluations = evaluations;
        self.num_vars -= dim;
        self.zero = RandomField::zero();
    }

    fn fixed_variables(&self, partial_point: &[RandomField<C>], config: C) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point, config);
        res
    }

    fn to_evaluations(&self) -> Vec<RandomField<C>> {
        let mut evaluations: Vec<_> = (0..1 << self.num_vars)
            .map(|_| RandomField::zero())
            .collect();
        self.evaluations
            .iter()
            .map(|(&i, v)| {
                evaluations[i] = v.clone();
            })
            .last();
        evaluations
    }
}

impl<C: ConfigReference> InSameField for SparseMultilinearExtension<C> {
    fn is_in_same_field(&self, other: &Self) -> bool {
        if self.evaluations.is_empty() || other.evaluations.is_empty() {
            return true;
        }

        unsafe {
            self.evaluations
                .first_key_value()
                .unwrap_unchecked()
                .1
                .is_in_same_field(other.evaluations.first_key_value().unwrap_unchecked().1)
        }
    }
}

impl<C: ConfigReference> Zero for SparseMultilinearExtension<C> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: tuples_to_treemap(&Vec::new()),
            zero: RandomField::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations.is_empty()
    }
}

impl<C: ConfigReference> Add for SparseMultilinearExtension<C> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        &self + &other
    }
}

impl<C: ConfigReference> Add for &SparseMultilinearExtension<C> {
    type Output = SparseMultilinearExtension<C>;

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

        assert!(
            self.is_in_same_field(rhs),
            "trying to add two sparse MLEs in different fields"
        );

        // simply merge the evaluations
        let mut evaluations = HashMap::new();
        for (&i, v) in self.evaluations.iter().chain(rhs.evaluations.iter()) {
            *evaluations.entry(i).or_insert(self.zero.clone()) += v;
        }
        let evaluations: Vec<_> = evaluations
            .into_iter()
            .filter(|(_, v)| !v.is_zero())
            .collect();

        Self::Output {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars: self.num_vars,
            zero: self.zero.clone(),
        }
    }
}

impl<C: ConfigReference> AddAssign for SparseMultilinearExtension<C> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + &other;
    }
}
impl<C: ConfigReference> AddAssign<&Self> for SparseMultilinearExtension<C> {
    fn add_assign(&mut self, rhs: &Self) {
        *self = &*self + rhs;
    }
}
impl<C: ConfigReference> AddAssign<(RandomField<C>, &Self)> for SparseMultilinearExtension<C> {
    fn add_assign(&mut self, (r, other): (RandomField<C>, &Self)) {
        if !self.is_zero() && !other.is_zero() {
            assert_eq!(
                other.num_vars, self.num_vars,
                "trying to add non-zero polynomial with different number of variables"
            );
            assert!(
                self.is_in_same_field(other),
                "trying to add two MLEs in different fields"
            );
        }
        let ev: Vec<_> = cfg_iter!(other.evaluations)
            .map(|(i, v)| (*i, r.clone() * v))
            .collect();
        let other = Self {
            num_vars: other.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: RandomField::zero(),
        };
        *self += &other;
    }
}
impl<C: ConfigReference> Neg for SparseMultilinearExtension<C> {
    type Output = Self;

    fn neg(self) -> Self {
        let ev: Vec<_> = cfg_iter!(self.evaluations)
            .map(|(i, v)| (*i, -v.clone()))
            .collect();
        Self::Output {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: self.zero,
        }
    }
}
impl<C: ConfigReference> Sub for SparseMultilinearExtension<C> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}
impl<C: ConfigReference> Sub for &SparseMultilinearExtension<C> {
    type Output = SparseMultilinearExtension<C>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: Self) -> Self::Output {
        self + &rhs.clone().neg()
    }
}
impl<C: ConfigReference> SubAssign for SparseMultilinearExtension<C> {
    fn sub_assign(&mut self, other: Self) {
        *self = &*self - &other;
    }
}
impl<C: ConfigReference> SubAssign<&Self> for SparseMultilinearExtension<C> {
    fn sub_assign(&mut self, rhs: &Self) {
        *self = &*self - rhs;
    }
}

impl<C: ConfigReference> Index<usize> for SparseMultilinearExtension<C> {
    type Output = RandomField<C>;

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
fn tuples_to_treemap<C: ConfigReference>(
    tuples: &[(usize, RandomField<C>)],
) -> BTreeMap<usize, RandomField<C>> {
    BTreeMap::from_iter(tuples.iter().map(|(i, v)| (*i, v.clone())))
}

fn treemap_to_hashmap<C: ConfigReference>(
    map: &BTreeMap<usize, RandomField<C>>,
) -> HashMap<usize, RandomField<C>> {
    HashMap::from_iter(map.iter().map(|(i, v)| (*i, v.clone())))
}
fn hashmap_to_treemap<C: ConfigReference>(
    map: &HashMap<usize, RandomField<C>>,
) -> BTreeMap<usize, RandomField<C>> {
    BTreeMap::from_iter(map.iter().map(|(i, v)| (*i, v.clone())))
}

// precompute  f(x) = eq(g,x)
fn precompute_eq<C: ConfigReference>(g: &[RandomField<C>], config: C) -> Vec<RandomField<C>> {
    let dim = g.len();
    let mut dp = vec![RandomField::zero(); 1 << dim];
    dp[0] = C::B::one().map_to_field(config) - g[0].clone();
    dp[1] = g[0].clone();
    for i in 1..dim {
        for b in 0..1 << i {
            let prev = dp[b].clone();
            dp[b + (1 << i)] = prev.clone() * &g[i];
            dp[b] = prev - dp[b + (1 << i)].clone();
        }
    }
    dp
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;
    use crate::{
        field::{ConfigRef, FieldConfig, RandomField},
        traits::Config,
    };

    // Function to convert usize to a binary vector of Ring elements.
    fn usize_to_binary_vector<C: ConfigReference>(
        n: usize,
        dimensions: usize,
        config: C,
    ) -> Vec<RandomField<C>> {
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
    fn boolean_hypercube<C: ConfigReference>(
        dimensions: usize,
        config: C,
    ) -> Vec<Vec<RandomField<C>>> {
        let max_val = 1 << dimensions; // 2^dimensions
        (0..max_val)
            .map(|i| usize_to_binary_vector(i, dimensions, config))
            .collect()
    }

    fn vec_cast<C: ConfigReference>(v: &[usize], config: C) -> Vec<RandomField<C>> {
        v.iter().map(|c| (*c as u64).map_to_field(config)).collect()
    }

    fn matrix_cast<C: ConfigReference>(
        m: &[Vec<usize>],
        config: C,
    ) -> SparseMatrix<RandomField<C>> {
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

    fn get_test_z<C: ConfigReference>(input: usize, config: C) -> Vec<RandomField<C>> {
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
        let config_ptr: ConfigRef<N> = ConfigRef::from(&config);
        let A = matrix_cast(
            &[
                vec![2, 3, 4, 4],
                vec![4, 11, 14, 14],
                vec![2, 8, 17, 17],
                vec![420, 4, 2, 0],
            ],
            config_ptr,
        );

        let A_mle = SparseMultilinearExtension::from_matrix(&A);
        assert_eq!(A_mle.evaluations.len(), 15); // 15 non-zero elements
        assert_eq!(A_mle.num_vars, 4); // 4x4 matrix, thus 2bit x 2bit, thus 2^4=16 evals

        let A = matrix_cast(
            &[
                vec![2, 3, 4, 4, 1],
                vec![4, 11, 14, 14, 2],
                vec![2, 8, 17, 17, 3],
                vec![420, 4, 2, 0, 4],
                vec![420, 4, 2, 0, 5],
            ],
            config_ptr,
        );
        let A_mle = SparseMultilinearExtension::from_matrix(&A);
        assert_eq!(A_mle.evaluations.len(), 23); // 23 non-zero elements
        assert_eq!(A_mle.num_vars, 6); // 5x5 matrix, thus 3bit x 3bit, thus 2^6=64 evals
    }

    #[test]
    fn test_vec_to_mle() {
        const N: usize = 1;
        let config = FieldConfig::new(293u32.into());
        let config = ConfigRef::<N>::from(&config);
        let z = get_test_z(3, config);

        let n_vars = 3;
        let z_mle = SparseMultilinearExtension::from_slice(n_vars, &z);

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
