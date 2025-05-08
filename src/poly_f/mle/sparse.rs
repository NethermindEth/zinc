use ark_ff::Zero;

use crate::biginteger::BigInt;
use crate::field::conversion::FieldMap;
use crate::field_config::FieldConfig;

use ark_std::{
    cfg_iter,
    collections::BTreeMap,
    ops::{Add, AddAssign, Index, Neg, Sub, SubAssign},
    vec::*,
};
#[cfg(feature = "parallel")]
use rayon::iter::*;

use crate::field::RandomField;

use super::MultilinearExtension;

use hashbrown::HashMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SparseMultilinearExtension<const N: usize> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: BTreeMap<usize, RandomField<N>>,
    /// Number of variables
    pub num_vars: usize,
    zero: RandomField<N>,
    /// Field in which the MLE is operating
    pub config: *const FieldConfig<N>,
}

impl<const N: usize> MultilinearExtension<N> for SparseMultilinearExtension<N> {
    fn fix_variables(&mut self, partial_point: &[RandomField<N>], config: *const FieldConfig<N>) {
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
                let dst_entry = result.entry(new_idx).or_insert(RandomField::zero());
                *dst_entry += &(gz * src_entry.1);
            }
            last = result;
        }
        let evaluations = hashmap_to_treemap(&last);

        self.evaluations = evaluations;
        self.num_vars -= dim;
        self.zero = RandomField::<N>::zero();
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
}
impl<const N: usize> Zero for SparseMultilinearExtension<N> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: tuples_to_treemap(&Vec::new()),
            zero: RandomField::<N>::zero(),
            config: std::ptr::null(),
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations.is_empty()
    }
}
impl<const N: usize> Add for SparseMultilinearExtension<N> {
    type Output = SparseMultilinearExtension<N>;

    fn add(self, other: SparseMultilinearExtension<N>) -> Self {
        &self + &other
    }
}
impl<'a, const N: usize> Add<&'a SparseMultilinearExtension<N>> for &SparseMultilinearExtension<N> {
    type Output = SparseMultilinearExtension<N>;

    fn add(self, rhs: &'a SparseMultilinearExtension<N>) -> Self::Output {
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
            *evaluations.entry(i).or_insert(RandomField::<N>::zero()) += &v;
        }
        let evaluations: Vec<_> = evaluations
            .into_iter()
            .filter(|(_, v)| !v.is_zero())
            .collect();

        Self::Output {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars: self.num_vars,
            zero: RandomField::<N>::zero(),
            config: self.config,
        }
    }
}

impl<const N: usize> AddAssign for SparseMultilinearExtension<N> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + &other;
    }
}
impl<'a, const N: usize> AddAssign<&'a SparseMultilinearExtension<N>>
    for SparseMultilinearExtension<N>
{
    fn add_assign(&mut self, rhs: &'a SparseMultilinearExtension<N>) {
        *self = &*self + rhs;
    }
}
impl<const N: usize> AddAssign<(RandomField<N>, &SparseMultilinearExtension<N>)>
    for SparseMultilinearExtension<N>
{
    fn add_assign(&mut self, (r, other): (RandomField<N>, &SparseMultilinearExtension<N>)) {
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
            zero: RandomField::<N>::zero(),
            config: self.config,
        };
        *self += &other;
    }
}
impl<const N: usize> Neg for SparseMultilinearExtension<N> {
    type Output = SparseMultilinearExtension<N>;

    fn neg(self) -> Self {
        let ev: Vec<_> = cfg_iter!(self.evaluations)
            .map(|(i, v)| (*i, -*v))
            .collect();
        Self::Output {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: RandomField::<N>::zero(),
            config: self.config,
        }
    }
}
impl<const N: usize> Sub for SparseMultilinearExtension<N> {
    type Output = SparseMultilinearExtension<N>;

    fn sub(self, other: SparseMultilinearExtension<N>) -> Self {
        &self - &other
    }
}
impl<'a, const N: usize> Sub<&'a SparseMultilinearExtension<N>> for &SparseMultilinearExtension<N> {
    type Output = SparseMultilinearExtension<N>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: &'a SparseMultilinearExtension<N>) -> Self::Output {
        self + &rhs.clone().neg()
    }
}
impl<const N: usize> SubAssign for SparseMultilinearExtension<N> {
    fn sub_assign(&mut self, other: SparseMultilinearExtension<N>) {
        *self = &*self - &other;
    }
}
impl<'a, const N: usize> SubAssign<&'a SparseMultilinearExtension<N>>
    for SparseMultilinearExtension<N>
{
    fn sub_assign(&mut self, rhs: &'a SparseMultilinearExtension<N>) {
        *self = &*self - rhs;
    }
}
impl<const N: usize> Index<usize> for SparseMultilinearExtension<N> {
    type Output = RandomField<N>;

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
fn tuples_to_treemap<const N: usize>(
    tuples: &[(usize, RandomField<N>)],
) -> BTreeMap<usize, RandomField<N>> {
    BTreeMap::from_iter(tuples.iter().map(|(i, v)| (*i, *v)))
}

fn treemap_to_hashmap<const N: usize>(
    map: &BTreeMap<usize, RandomField<N>>,
) -> HashMap<usize, RandomField<N>> {
    HashMap::from_iter(map.iter().map(|(i, v)| (*i, *v)))
}
fn hashmap_to_treemap<const N: usize>(
    map: &HashMap<usize, RandomField<N>>,
) -> BTreeMap<usize, RandomField<N>> {
    BTreeMap::from_iter(map.iter().map(|(i, v)| (*i, *v)))
}

// precompute  f(x) = eq(g,x)
fn precompute_eq<const N: usize>(
    g: &[RandomField<N>],
    config: *const FieldConfig<N>,
) -> Vec<RandomField<N>> {
    let dim = g.len();
    let mut dp = vec![RandomField::zero(); 1 << dim];
    dp[0] = BigInt::<N>::one().map_to_field(config) - g[0];
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
