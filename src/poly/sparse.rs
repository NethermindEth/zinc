use std::{
    collections::{BTreeMap, HashMap},
    ops::{Add, AddAssign, Index, Neg, SubAssign},
};

use ark_std::log2;
use num_traits::Zero;
use rand_core::RngCore;

use crate::{
    poly::mle::{MultilinearExtension, swap_bits},
    sparse_matrix::SparseMatrix,
    traits::Ring,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseMultilinearExtension<F> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: BTreeMap<usize, F>,
    /// Number of variables
    pub num_vars: usize,
    zero: F,
}

impl<F: Ring> SparseMultilinearExtension<F> {
    pub fn from_evaluations<'a>(
        num_vars: usize,
        evaluations: impl IntoIterator<Item = &'a (usize, F)>,
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
            zero: F::ZERO,
        }
    }

    pub fn evaluate(&self, point: &[F]) -> F {
        assert_eq!(point.len(), self.num_vars);
        self.fixed_variables(point)[0]
    }

    /// Returns the sparse MLE from the given matrix, without modifying the original matrix.
    pub fn from_matrix(m: &SparseMatrix<F>) -> Self {
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
        Self::from_sparse_slice(n_vars, &v)
    }

    /// Takes n_vars and a sparse slice and returns its sparse MLE.
    pub fn from_sparse_slice(n_vars: usize, v: &[(usize, F)]) -> Self {
        SparseMultilinearExtension::<F>::from_evaluations(n_vars, v)
    }

    /// Takes n_vars and a dense slice and returns its sparse MLE.
    pub fn from_slice(n_vars: usize, v: &[F]) -> Self {
        let v_sparse = v
            .iter()
            .enumerate()
            .map(|(i, v_i)| (i, *v_i))
            .collect::<Vec<(usize, F)>>();
        SparseMultilinearExtension::<F>::from_evaluations(n_vars, &v_sparse)
    }
}

impl<F: Ring> MultilinearExtension<F> for SparseMultilinearExtension<F> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn rand<R: RngCore + ?Sized>(num_vars: usize, rng: &mut R) -> Self {
        // choose about sqrt(2^n) non-zero entries uniformly at random
        let nnz = 1 << (num_vars / 2);
        assert!(nnz <= 1 << num_vars);
        let mut map = HashMap::new();
        let mask = (1usize << num_vars) - 1;
        for _ in 0..nnz {
            let mut index = (rng.next_u64() as usize) & mask;
            while map.contains_key(&index) {
                index = (rng.next_u64() as usize) & mask;
            }
            map.entry(index).or_insert(F::random(rng));
        }
        let evaluations = hashmap_to_treemap(&map);
        Self {
            num_vars,
            evaluations,
            zero: F::ZERO,
        }
    }

    fn relabel(&self, mut a: usize, mut b: usize, k: usize) -> Self {
        if a > b {
            core::mem::swap(&mut a, &mut b);
        }
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
            .map(|(&i, v)| (swap_bits(i, a, b, k), *v))
            .collect();
        Self {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: F::ZERO,
        }
    }

    fn fix_variables(&mut self, partial_point: &[F]) {
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
            for (old_idx, val) in last.iter() {
                let gz = pre[*old_idx & ((1 << dim) - 1)];
                let new_idx = *old_idx >> dim;
                let dst_entry = result.entry(new_idx).or_insert(F::ZERO);
                *dst_entry += gz * val;
            }
            last = result;
        }
        let evaluations = hashmap_to_treemap(&last);
        self.evaluations = evaluations;
        self.num_vars -= dim;
        self.zero = F::ZERO;
    }

    fn fixed_variables(&self, partial_point: &[F]) -> Self {
        let mut res = self.clone();
        res.fix_variables(partial_point);
        res
    }

    fn to_evaluations(&self) -> Vec<F> {
        let mut evaluations: Vec<_> = (0..1 << self.num_vars).map(|_| F::ZERO).collect();
        for (&i, v) in self.evaluations.iter() {
            evaluations[i] = *v;
        }
        evaluations
    }
}

impl<F: Ring> Add for SparseMultilinearExtension<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<F: Ring> Add<&Self> for SparseMultilinearExtension<F> {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        // handle zero case by simple merge
        let mut merged: HashMap<usize, F> = treemap_to_hashmap(&self.evaluations);
        for (&i, v) in rhs.evaluations.iter() {
            let e = merged.entry(i).or_insert(F::ZERO);
            *e += v;
        }
        let evaluations: Vec<_> = merged.into_iter().filter(|(_, v)| !v.is_zero()).collect();
        Self {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars: self.num_vars,
            zero: F::ZERO,
        }
    }
}

impl<F: Ring> AddAssign<&Self> for SparseMultilinearExtension<F> {
    fn add_assign(&mut self, rhs: &Self) {
        let mut merged: HashMap<usize, F> = treemap_to_hashmap(&self.evaluations);
        for (&i, v) in rhs.evaluations.iter() {
            let e = merged.entry(i).or_insert(F::ZERO);
            *e += v;
        }
        self.evaluations = hashmap_to_treemap(&merged);
    }
}

impl<F: Ring> SubAssign<&Self> for SparseMultilinearExtension<F> {
    fn sub_assign(&mut self, rhs: &Self) {
        *self += (F::ZERO - F::one(), rhs)
    }
}

impl<F: Ring> AddAssign<(F, &Self)> for SparseMultilinearExtension<F> {
    fn add_assign(&mut self, (r, other): (F, &Self)) {
        let mut scaled = other.clone();
        for v in scaled.evaluations.values_mut() {
            *v = *v * r;
        }
        *self += &scaled;
    }
}

impl<F: Ring> Neg for SparseMultilinearExtension<F> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for v in self.evaluations.values_mut() {
            *v = v.negate();
        }
        self
    }
}

impl<F: Ring> Index<usize> for SparseMultilinearExtension<F> {
    type Output = F;

    fn index(&self, index: usize) -> &Self::Output {
        self.evaluations.get(&index).unwrap_or(&self.zero)
    }
}

impl<F: Ring> Zero for SparseMultilinearExtension<F> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: BTreeMap::new(),
            zero: F::ZERO,
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations.is_empty()
    }
}

// Utilities
fn tuples_to_treemap<F: Ring>(tuples: &[(usize, F)]) -> BTreeMap<usize, F> {
    BTreeMap::from_iter(tuples.iter().map(|(i, v)| (*i, *v)))
}

fn treemap_to_hashmap<F: Ring>(map: &BTreeMap<usize, F>) -> HashMap<usize, F> {
    HashMap::from_iter(map.iter().map(|(i, v)| (*i, *v)))
}
fn hashmap_to_treemap<F: Ring>(map: &HashMap<usize, F>) -> BTreeMap<usize, F> {
    BTreeMap::from_iter(map.iter().map(|(i, v)| (*i, *v)))
}

fn precompute_eq<F: Ring>(g: &[F]) -> Vec<F> {
    let dim = g.len();
    let mut dp = vec![F::ZERO; 1 << dim];
    dp[0] = F::one() - g[0];
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
mod tests {
    use crypto_bigint::{U128, const_monty_params};
    use num_traits::{One, Zero};
    use proptest::prelude::*;

    use super::*;
    use crate::{field::RandomField, poly::mle::MultilinearExtension, sparse_matrix::SparseMatrix};

    const_monty_params!(ModP, U128, "0076F668F4274572E39A3EA8285319B5");
    type Fp = RandomField<ModP, { U128::LIMBS }>;

    fn any_f() -> impl Strategy<Value = Fp> {
        any::<u128>().prop_map(Fp::from)
    }

    fn any_sme() -> impl Strategy<Value = SparseMultilinearExtension<Fp>> {
        (0usize..=5).prop_flat_map(|n| {
            let len = 1usize << n;
            prop::collection::vec(any_f(), len).prop_map(move |evals| {
                // Build sparse by dropping zeros
                let tuples: Vec<(usize, Fp)> = evals
                    .iter()
                    .enumerate()
                    .filter_map(|(i, v)| if v.is_zero() { None } else { Some((i, *v)) })
                    .collect();
                SparseMultilinearExtension::from_evaluations(n, &tuples)
            })
        })
    }

    #[test]
    fn test_from_slice_and_indexing_sparse() {
        let n_vars = 3usize;
        let v = vec![Fp::from(1u64), Fp::from(0u64), Fp::from(3u64)];
        let sparse = SparseMultilinearExtension::from_slice(n_vars, &v);
        assert_eq!(sparse.num_vars(), n_vars);
        let mut expected = v.clone();
        expected.resize(1 << n_vars, Fp::zero());
        assert_eq!(sparse.to_evaluations(), expected);
        assert_eq!(sparse[0], Fp::from(1u64));
        assert_eq!(sparse[1], Fp::zero());
    }

    #[test]
    fn test_from_sparse_slice_and_matrix() {
        // from_sparse_slice
        let tuples = vec![(0usize, Fp::from(5u64)), (5usize, Fp::from(7u64))];
        let sme = SparseMultilinearExtension::from_sparse_slice(3, &tuples);
        let evals = sme.to_evaluations();
        assert_eq!(evals[0], Fp::from(5u64));
        assert_eq!(evals[5], Fp::from(7u64));
        assert!(evals.iter().enumerate().all(|(i, v)| if i == 0 || i == 5 {
            true
        } else {
            v.is_zero()
        }));

        let m: SparseMatrix<Fp> = SparseMatrix::from(vec![
            vec![Fp::from(5u64), Fp::zero()],
            vec![Fp::zero(), Fp::zero()],
            vec![Fp::zero(), Fp::from(7u64)],
        ]);
        let sm2 = SparseMultilinearExtension::from_matrix(&m);
        assert_eq!(sm2.num_vars(), 3);
        let evals2 = sm2.to_evaluations();
        assert_eq!(evals2[0], Fp::from(5u64));
        assert_eq!(evals2[5], Fp::from(7u64));
        assert!(evals2.iter().enumerate().all(|(i, v)| if i == 0 || i == 5 {
            true
        } else {
            v.is_zero()
        }));
    }

    #[test]
    fn test_fix_variables_and_evaluate_sparse() {
        let evals = [
            Fp::from(10u64),
            Fp::from(20u64),
            Fp::from(30u64),
            Fp::from(40u64),
        ];
        let mut tuples: Vec<(usize, Fp)> = vec![];
        for (i, v) in evals.iter().cloned().enumerate() {
            if !v.is_zero() {
                tuples.push((i, v));
            }
        }
        let sme = SparseMultilinearExtension::from_evaluations(2, &tuples);

        for (idx, &(x0, x1)) in [
            (Fp::zero(), Fp::zero()),
            (Fp::one(), Fp::zero()),
            (Fp::zero(), Fp::one()),
            (Fp::one(), Fp::one()),
        ]
        .iter()
        .enumerate()
        {
            let val = sme.evaluate(&[x0, x1]);
            assert_eq!(val, evals[idx]);
        }

        let mut m2 = sme.clone();
        m2.fix_variables(&[Fp::one()]);
        assert_eq!(m2.num_vars(), 1);
        assert_eq!(m2.to_evaluations(), vec![Fp::from(20u64), Fp::from(40u64)]);
    }

    #[test]
    fn test_relabel_properties_sparse() {
        let tuples: Vec<(usize, Fp)> = vec![(0, 1u64.into()), (1, 2u64.into()), (2, 3u64.into())];
        let sme = SparseMultilinearExtension::from_evaluations(3, &tuples);
        let out = sme.relabel(0, 1, 1);
        let out2 = sme.relabel(1, 0, 1);
        assert_eq!(out.to_evaluations(), out2.to_evaluations());

        let twice = sme.relabel(0, 1, 1).relabel(0, 1, 1);
        assert_eq!(twice.to_evaluations(), sme.to_evaluations());
    }

    #[test]
    fn test_zero_and_neg_and_add_assign_sparse() {
        let z: SparseMultilinearExtension<Fp> = Zero::zero();
        assert_eq!(z.num_vars, 0);
        assert!(z.is_zero());

        let a = SparseMultilinearExtension::from_evaluations(
            2,
            &vec![(0usize, 1u64.into()), (3, 4u64.into())],
        );
        let b = SparseMultilinearExtension::from_evaluations(
            2,
            &vec![(0usize, 5u64.into()), (1, 6u64.into())],
        );

        let mut c = a.clone();
        c += &b;
        let mut expected = a.to_evaluations();
        let bev = b.to_evaluations();
        for i in 0..expected.len() {
            expected[i] += &bev[i];
        }
        assert_eq!(c.to_evaluations(), expected);

        let mut d = a.clone();
        d += (Fp::from(2u64), &b);
        let mut expected2 = a.to_evaluations();
        for i in 0..expected2.len() {
            expected2[i] += Fp::from(2u64) * bev[i];
        }
        assert_eq!(d.to_evaluations(), expected2);

        let neg = -a.clone();
        let mut expected_neg = a.to_evaluations();
        for v in expected_neg.iter_mut() {
            *v = -*v;
        }
        assert_eq!(neg.to_evaluations(), expected_neg);
    }

    fn any_aligned_pair_with_point_sparse() -> impl Strategy<
        Value = (
            SparseMultilinearExtension<Fp>,
            SparseMultilinearExtension<Fp>,
            Vec<Fp>,
        ),
    > {
        (0usize..=5).prop_flat_map(|n| {
            let len = 1usize << n;
            prop::collection::vec(any_f(), len).prop_flat_map(move |e1| {
                let n2 = n;
                prop::collection::vec(any_f(), len).prop_flat_map(move |e2| {
                    let n3 = n2;
                    prop::collection::vec(any_f(), n3).prop_map({
                        let e1v = e1.clone();
                        let e2v = e2.clone();
                        move |r| {
                            let t1: Vec<(usize, Fp)> = e1v
                                .iter()
                                .cloned()
                                .enumerate()
                                .filter_map(|(i, v)| if v.is_zero() { None } else { Some((i, v)) })
                                .collect();
                            let t2: Vec<(usize, Fp)> = e2v
                                .iter()
                                .cloned()
                                .enumerate()
                                .filter_map(|(i, v)| if v.is_zero() { None } else { Some((i, v)) })
                                .collect();
                            (
                                SparseMultilinearExtension::from_evaluations(n3, &t1),
                                SparseMultilinearExtension::from_evaluations(n3, &t2),
                                r,
                            )
                        }
                    })
                })
            })
        })
    }

    proptest! {
        #[test]
        fn prop_eval_add_is_linear_sparse((p1, p2, r) in any_aligned_pair_with_point_sparse()) {
            let lhs = (p1.clone() + &p2).evaluate(&r);
            let rhs = p1.evaluate(&r) + p2.evaluate(&r);
            prop_assert_eq!(lhs, rhs);
        }

        #[test]
        fn prop_fix_vars_commutes_with_eval_sparse(p in any_sme(), r in prop::collection::vec(any_f(), 0..=6), k in 0usize..=6) {
            let n = p.num_vars();
            prop_assume!(r.len() >= n);
            let r = &r[..n];
            let k = k.min(n);
            let mut pfixed = p.clone();
            pfixed.fix_variables(&r[..k]);
            let lhs = pfixed.evaluate(&r[k..]);
            let rhs = p.evaluate(r);
            prop_assert_eq!(lhs, rhs);
        }

        #[test]
        fn prop_relabel_is_involution_sparse(p in any_sme().prop_filter("n>=3", |p| p.num_vars()>=3)) {
            let q = p.relabel(0,1,1).relabel(0,1,1);
            prop_assert_eq!(q.to_evaluations(), p.to_evaluations());
        }
    }
}
