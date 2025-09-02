use std::{
    collections::BTreeSet,
    fmt::Debug,
    ops::{Add, AddAssign, Div, Mul, Sub},
};

use crypto_bigint::{
    Int, NonZero, Random, Uint, Word,
    subtle::{ConstantTimeEq, CtOption},
};
use num_traits::{ConstZero, One, Zero};

use crate::field::Monty;

pub trait Value:
    Copy + Clone + Debug + PartialEq + Eq + Send + Sync + Sized + Random + Default
{
}

impl<T> Value for T where
    T: Copy + Clone + Debug + PartialEq + Eq + Send + Sync + Sized + Random + Default
{
}

pub trait Abelian:
    Value
    + Add<Self, Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + AddAssign<Self>
    + for<'a> AddAssign<&'a Self>
    + Add<Self, Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + AddAssign<Self>
    + for<'a> AddAssign<&'a Self>
    + Sub<Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + Negate<Output = Self>
    + Zero
    + ConstZero
{
}

impl<T> Abelian for T where
    T: Value
        + Add<Self, Output = Self>
        + for<'a> Add<&'a Self, Output = Self>
        + AddAssign<Self>
        + for<'a> AddAssign<&'a Self>
        + Add<Self, Output = Self>
        + for<'a> Add<&'a Self, Output = Self>
        + AddAssign<Self>
        + for<'a> AddAssign<&'a Self>
        + Sub<Self, Output = Self>
        + for<'a> Sub<&'a Self, Output = Self>
        + Negate<Output = Self>
        + Zero
        + ConstZero
{
}

pub trait Monoid:
    Value
    + Mul<Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    // + MulAssign<Self>
    // + for<'a> MulAssign<&'a Self>
    + Div<NonZero<Self>, Output = CtOption<Self>>
    + One
    + Zero
    + ConstantTimeEq
{
}

impl<T> Monoid for T where
    T: Value
        + Mul<Self, Output = Self>
        + for<'a> Mul<&'a Self, Output = Self>
        // + MulAssign<Self>
        // + for<'a> MulAssign<&'a Self>
        + Div<NonZero<Self>, Output = CtOption<Self>>
        + One
        + Zero
        + ConstantTimeEq
{
}

pub trait Ring: Monoid + Abelian {}

impl<T> Ring for T where T: Monoid + Abelian {}

pub trait Transcript<const LIMBS: usize> {
    fn get_encoding_element(&mut self) -> Int<LIMBS>;
    fn get_word(&mut self) -> Word;
    fn sample_unique_columns(
        &mut self,
        range: ark_std::ops::Range<usize>,
        columns: &mut BTreeSet<usize>,
        count: usize,
    ) -> usize;
}

pub trait ToBytes {
    fn to_be_bytes(&self) -> Vec<u8>;
    fn to_le_bytes(&self) -> Vec<u8>;
}

pub trait FromBytes: Sized {
    /// Constructs an instance from a little-endian byte slice.
    fn from_bytes_le(bytes: &[u8]) -> Option<Self>;

    /// Constructs an instance from a big-endian byte slice.
    fn from_bytes_be(bytes: &[u8]) -> Option<Self>;
}

pub trait FromBits {
    fn from_be_bits(bits: &[bool]) -> Self;
    fn from_le_bits(bits: &[bool]) -> Self;
}

pub trait Negate {
    type Output;
    fn negate(&self) -> Self::Output;
}

pub trait Field<const N: usize>:
    Ring
    + Default
    + ToBytes
    + FromBytes
    + From<u64>
    + From<u128>
    + From<i64>
    + From<i128>
    + From<Uint<N>>
    + From<Int<N>>
    + for<'a> From<&'a Self>
    + for<'a> AddAssign<&'a Self>
    + MapIterable
{
    type Monty: Monty<N>;
}

pub trait MapIterable: Sized {
    fn map_iterable<'a, const M: usize, I: IntoIterator<Item = &'a Int<M>>>(
        iterable: I,
    ) -> Vec<Self>;
}
