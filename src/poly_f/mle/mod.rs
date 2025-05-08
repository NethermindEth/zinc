pub(crate) mod dense;
pub(crate) mod sparse;

use ark_std::{
    fmt::Debug,
    ops::{Add, AddAssign, Index, Neg, SubAssign},
    Zero,
};

/// This trait describes an interface for the multilinear extension
/// of an array.
/// The latter is a multilinear polynomial represented in terms of its
/// evaluations over the domain {0,1}^`num_vars` (i.e. the Boolean hypercube).
///
/// Index represents a point, which is a vector in {0,1}^`num_vars` in little
/// endian form. For example, `0b1011` represents `P(1,1,0,1)`
pub trait MultilinearExtension<const N: usize>:
    Sized
    + Clone
    + Debug
    + PartialEq
    + Eq
    + Add
    + Neg
    + Zero
    + for<'a> AddAssign<&'a Self>
    + for<'a> AddAssign<(RandomField<N>, &'a Self)>
    + for<'a> SubAssign<&'a Self>
    + Index<usize>
{
    /// Reduce the number of variables of `self` by fixing the
    /// `partial_point.len()` variables at `partial_point`.
    fn fix_variables(&mut self, partial_point: &[RandomField<N>], config: *const FieldConfig<N>);

    /// Creates a new object with the number of variables of `self` reduced by fixing the
    /// `partial_point.len()` variables at `partial_point`.
    fn fixed_variables(
        &self,
        partial_point: &[RandomField<N>],
        config: *const FieldConfig<N>,
    ) -> Self;
}
/// swap the bits of `x` from position `a..a+n` to `b..b+n` and from `b..b+n` to `a..a+n` in little endian order
pub(crate) fn swap_bits(x: usize, a: usize, b: usize, n: usize) -> usize {
    let a_bits = (x >> a) & ((1usize << n) - 1);
    let b_bits = (x >> b) & ((1usize << n) - 1);
    let local_xor_mask = a_bits ^ b_bits;
    let global_xor_mask = (local_xor_mask << a) | (local_xor_mask << b);
    x ^ global_xor_mask
}

use crate::field::RandomField;
use crate::field_config::FieldConfig;
