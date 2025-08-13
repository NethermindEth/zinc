use ark_ff::UniformRand;
use ark_std::{
    collections::BTreeMap,
    log2,
    ops::{Add, AddAssign, Index, Neg, Sub, SubAssign},
    rand::Rng,
    vec,
    vec::Vec,
};
use hashbrown::HashMap;
use num_traits::Zero;

use super::{MultilinearExtension, swap_bits};
use crate::{sparse_matrix::SparseMatrix, traits::Integer};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseMultilinearExtension<I> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: BTreeMap<usize, I>,
    /// Number of variables
    pub num_vars: usize,
    zero: I,
}

impl<I: Integer> SparseMultilinearExtension<I> {
    pub fn from_evaluations<'a>(
        num_vars: usize,
        evaluations: impl IntoIterator<Item = &'a (usize, I)>,
    ) -> Self
    where
        I: 'a,
    {
        let bit_mask = 1 << num_vars;

        let evaluations: Vec<_> = evaluations
            .into_iter()
            .map(|(i, v): &(usize, I)| {
                assert!(*i < bit_mask, "index out of range");
                (*i, v.clone())
            })
            .collect();
        Self {
            evaluations: tuples_to_treemap::<I>(&evaluations),
            num_vars,
            zero: I::ZERO,
        }
    }
    pub fn evaluate(&self, point: &[I]) -> I {
        assert_eq!(point.len(), self.num_vars);
        self.fixed_variables(point)[0].clone()
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
        rng: &mut Rn,
    ) -> Self {
        assert!(num_nonzero_entries <= 1 << num_vars);

        let mut map = HashMap::new();
        for _ in 0..num_nonzero_entries {
            let mut index = usize::rand(rng) & ((1 << num_vars) - 1);
            while map.contains_key(&index) {
                index = usize::rand(rng) & ((1 << num_vars) - 1);
            }
            map.entry(index).or_insert(I::random(rng));
        }
        let mut buf = Vec::new();
        for (arg, v) in map.iter() {
            if *v != I::ZERO {
                buf.push((*arg, v.clone()));
            }
        }
        let evaluations = hashmap_to_treemap(&map);
        Self {
            num_vars,
            evaluations,
            zero: I::ZERO,
        }
    }

    /// Returns the sparse MLE from the given matrix, without modifying the original matrix.
    pub fn from_matrix(m: &SparseMatrix<I>) -> Self {
        let n_rows = m.n_rows.next_power_of_two();
        let n_cols = m.n_cols.next_power_of_two();
        let n_vars: usize = (log2(n_rows * n_cols)) as usize; // n_vars = s + s'

        // build the sparse vec representing the sparse matrix
        let total_elements: usize = m.coeffs.iter().map(|row| row.len()).sum();
        let mut v: Vec<(usize, I)> = Vec::with_capacity(total_elements);

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
    pub fn from_sparse_slice(n_vars: usize, v: &[(usize, I)]) -> Self {
        SparseMultilinearExtension::from_evaluations(n_vars, v)
    }

    /// Takes n_vars and a dense slice and returns its sparse MLE.
    pub fn from_slice(n_vars: usize, v: &[I]) -> Self {
        let v_sparse = v
            .iter()
            .enumerate()
            .map(|(i, v_i)| (i, v_i.clone()))
            .collect::<Vec<(usize, I)>>();
        SparseMultilinearExtension::from_evaluations(n_vars, &v_sparse)
    }
}

impl<I: Integer> MultilinearExtension<I> for SparseMultilinearExtension<I> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }
    /// Outputs an `l`-variate multilinear extension where value of evaluations
    /// are sampled uniformly at random. The number of nonzero entries is
    /// `sqrt(2^num_vars)` and indices of those nonzero entries are distributed
    /// uniformly at random.
    fn rand<Rn: ark_std::rand::Rng>(num_vars: usize, rng: &mut Rn) -> Self {
        Self::rand_with_config(num_vars, 1 << (num_vars / 2), rng)
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
        let ev: Vec<_> = self
            .evaluations
            .iter()
            .map(|(&i, v)| (swap_bits(i, a, b, k), v.clone()))
            .collect();
        Self {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: I::ZERO,
        }
    }

    fn fix_variables(&mut self, partial_point: &[I]) {
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
            let pre = precompute_eq(focus);
            let dim = focus.len();
            let mut result = HashMap::new();
            for src_entry in last.iter() {
                let old_idx = *src_entry.0;
                let gz = pre[old_idx & ((1 << dim) - 1)].clone();
                let new_idx = old_idx >> dim;
                let dst_entry = result.entry(new_idx).or_insert(I::ZERO);
                *dst_entry += &(gz * src_entry.1);
            }
            last = result;
        }
        let evaluations = hashmap_to_treemap(&last);

        self.evaluations = evaluations;
        self.num_vars -= dim;
    }

    fn fixed_variables(&self, partial_point: &[I]) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point);
        res
    }

    fn to_evaluations(&self) -> Vec<I> {
        let mut evaluations: Vec<_> = (0..1 << self.num_vars).map(|_| I::ZERO).collect();
        self.evaluations
            .iter()
            .map(|(&i, v)| {
                evaluations[i] = v.clone();
            })
            .last();
        evaluations
    }
}
impl<I: Integer> Zero for SparseMultilinearExtension<I> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: tuples_to_treemap(&Vec::new()),
            zero: I::ZERO,
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations.is_empty()
    }
}
impl<I: Integer> Add for SparseMultilinearExtension<I> {
    type Output = SparseMultilinearExtension<I>;

    fn add(self, other: SparseMultilinearExtension<I>) -> Self {
        &self + &other
    }
}
impl<'a, I: Integer> Add<&'a SparseMultilinearExtension<I>> for &SparseMultilinearExtension<I> {
    type Output = SparseMultilinearExtension<I>;

    fn add(self, rhs: &'a SparseMultilinearExtension<I>) -> Self::Output {
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

        // simply merge the evaluations
        let mut evaluations = HashMap::<_, I>::new();
        for (&i, v) in self.evaluations.iter().chain(rhs.evaluations.iter()) {
            *evaluations.entry(i).or_insert(I::ZERO) += v;
        }
        let evaluations: Vec<_> = evaluations
            .into_iter()
            .filter(|(_, v)| !<I as Zero>::is_zero(v))
            .collect();

        Self::Output {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars: self.num_vars,
            zero: I::ZERO,
        }
    }
}

impl<I: Integer> AddAssign for SparseMultilinearExtension<I> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + &other;
    }
}
impl<'a, I: Integer> AddAssign<&'a SparseMultilinearExtension<I>>
    for SparseMultilinearExtension<I>
{
    fn add_assign(&mut self, rhs: &'a SparseMultilinearExtension<I>) {
        *self = &*self + rhs;
    }
}
impl<I: Integer> AddAssign<(I, &SparseMultilinearExtension<I>)> for SparseMultilinearExtension<I> {
    fn add_assign(&mut self, (r, other): (I, &SparseMultilinearExtension<I>)) {
        if !self.is_zero() && !other.is_zero() {
            assert_eq!(
                other.num_vars, self.num_vars,
                "trying to add non-zero polynomial with different number of variables"
            );
        }
        let ev: Vec<_> = other
            .evaluations
            .iter()
            .map(|(i, v)| (*i, r.clone() * v))
            .collect();
        let other = Self {
            num_vars: other.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: I::ZERO,
        };
        *self += &other;
    }
}
impl<I: Integer> Neg for SparseMultilinearExtension<I> {
    type Output = SparseMultilinearExtension<I>;

    fn neg(self) -> Self {
        let ev: Vec<_> = self
            .evaluations
            .iter()
            .map(|(i, v)| (*i, I::ZERO - v))
            .collect();
        Self::Output {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: I::ZERO,
        }
    }
}
impl<I: Integer> Sub for SparseMultilinearExtension<I> {
    type Output = SparseMultilinearExtension<I>;

    fn sub(self, other: SparseMultilinearExtension<I>) -> Self {
        &self - &other
    }
}
impl<'a, I: Integer> Sub<&'a SparseMultilinearExtension<I>> for &SparseMultilinearExtension<I> {
    type Output = SparseMultilinearExtension<I>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: &'a SparseMultilinearExtension<I>) -> Self::Output {
        self + &rhs.clone().neg()
    }
}
impl<I: Integer> SubAssign for SparseMultilinearExtension<I> {
    fn sub_assign(&mut self, other: SparseMultilinearExtension<I>) {
        *self = &*self - &other;
    }
}
impl<'a, I: Integer> SubAssign<&'a SparseMultilinearExtension<I>>
    for SparseMultilinearExtension<I>
{
    fn sub_assign(&mut self, rhs: &'a SparseMultilinearExtension<I>) {
        *self = &*self - rhs;
    }
}
impl<I: Integer> Index<usize> for SparseMultilinearExtension<I> {
    type Output = I;

    /// Returns the evaluation of the polynomial at a point represented by
    /// index.
    ///
    /// Index represents a vector in {0,1}^`num_vars` in little endian form. For
    /// example, `0b1011` represents `P(1,1,0,1)`
    ///
    /// For Sparse multilinear polynomial, Lookup_evaluation takes log time to
    /// the size of polynomial.
    fn index(&self, index: usize) -> &Self::Output {
        self.evaluations.get(&index).unwrap_or(&self.zero)
    }
}

/// Utilities
fn tuples_to_treemap<I: Integer>(tuples: &[(usize, I)]) -> BTreeMap<usize, I> {
    BTreeMap::from_iter(tuples.iter().map(|(i, v)| (*i, v.clone())))
}

fn treemap_to_hashmap<I: Integer>(map: &BTreeMap<usize, I>) -> HashMap<usize, I> {
    HashMap::from_iter(map.iter().map(|(i, v)| (*i, v.clone())))
}
fn hashmap_to_treemap<I: Integer>(map: &HashMap<usize, I>) -> BTreeMap<usize, I> {
    BTreeMap::from_iter(map.iter().map(|(i, v)| (*i, v.clone())))
}

// precompute  f(x) = eq(g,x)
fn precompute_eq<I: Integer>(g: &[I]) -> Vec<I> {
    let dim = g.len();
    let mut dp = vec![I::ZERO; 1 << dim];
    dp[0] = I::one() - &g[0];
    dp[1] = g[0].clone();
    for i in 1..dim {
        for b in 0..1 << i {
            let prev = dp[b].clone();
            dp[b + (1 << i)] = prev.clone() * g[i].clone();
            dp[b] = prev - &dp[b + (1 << i)];
        }
    }
    dp
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use num_traits::ConstZero;

    use super::*;
    use crate::field::Int;
    // Function to convert usize to a binary vector of Ring elements.
    fn usize_to_binary_vector<I: Integer>(n: usize, dimensions: usize) -> Vec<I> {
        let mut bits = Vec::with_capacity(dimensions);
        let mut current = n;

        for _ in 0..dimensions {
            if (current & 1) == 1 {
                bits.push(I::one());
            } else {
                bits.push(I::ZERO);
            }
            current >>= 1;
        }
        bits
    }

    // Wrapper function to generate a boolean hypercube.
    fn boolean_hypercube<const N: usize>(dimensions: usize) -> Vec<Vec<Int<N>>> {
        let max_val = 1 << dimensions; // 2^dimensions
        (0..max_val)
            .map(|i| usize_to_binary_vector(i, dimensions))
            .collect()
    }

    fn get_test_z<const N: usize>(input: i64) -> Vec<Int<N>> {
        [
            input, // io
            1,
            input.pow(3) + input + 5, // x^3 + x + 5
            input.pow(2),             // x^2
            input.pow(3),             // x^2 * x
            input.pow(3) + input,     // x^3 + x
        ]
        .iter()
        .cloned()
        .map(Int::from)
        .collect()
    }
    fn matrix_cast<const N: usize>(m: &[Vec<usize>]) -> SparseMatrix<Int<N>> {
        let n_rows = m.len();
        let n_cols = m[0].len();
        let mut coeffs = Vec::with_capacity(n_rows);
        for row in m.iter() {
            let mut row_coeffs = Vec::with_capacity(n_cols);
            for (col_i, &val) in row.iter().enumerate() {
                if val != 0 {
                    row_coeffs.push((Int::from(val as i64), col_i));
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
    #[test]
    fn test_matrix_to_mle() {
        const N: usize = 1;

        let A = matrix_cast::<N>(&[
            vec![2, 3, 4, 4],
            vec![4, 11, 14, 14],
            vec![2, 8, 17, 17],
            vec![420, 4, 2, 0],
        ]);

        let A_mle = SparseMultilinearExtension::from_matrix(&A);
        assert_eq!(A_mle.evaluations.len(), 15); // 15 non-zero elements
        assert_eq!(A_mle.num_vars, 4); // 4x4 matrix, thus 2bit x 2bit, thus 2^4=16 evals

        let A = matrix_cast::<N>(&[
            vec![2, 3, 4, 4, 1],
            vec![4, 11, 14, 14, 2],
            vec![2, 8, 17, 17, 3],
            vec![420, 4, 2, 0, 4],
            vec![420, 4, 2, 0, 5],
        ]);
        let A_mle = SparseMultilinearExtension::from_matrix(&A);
        assert_eq!(A_mle.evaluations.len(), 23); // 23 non-zero elements
        assert_eq!(A_mle.num_vars, 6); // 5x5 matrix, thus 3bit x 3bit, thus 2^6=64 evals
    }

    #[test]
    fn test_vec_to_mle() {
        let z = get_test_z(3);
        const N: usize = 2;
        let n_vars = 3;
        let z_mle = SparseMultilinearExtension::from_slice(n_vars, &z);

        // check that the z_mle evaluated over the boolean hypercube equals the vec z_i values
        let bhc = boolean_hypercube(z_mle.num_vars);
        for (i, z_i) in z.iter().enumerate() {
            let s_i = &bhc[i];
            assert_eq!(z_mle.evaluate(s_i), z_i.clone());
        }
        // for the rest of elements of the boolean hypercube, expect it to evaluate to zero
        for s_i in bhc.iter().take(1 << z_mle.num_vars).skip(z.len()) {
            assert_eq!(z_mle.fixed_variables(s_i,)[0], Int::<N>::ZERO);
        }
    }
}
