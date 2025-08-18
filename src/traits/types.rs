use ark_std::{
    UniformRand,
    fmt::Debug,
    iter::Sum,
    ops::{
        Add, AddAssign, Index, IndexMut, Mul, Neg, Range, RangeTo, RemAssign, Shr, Sub, SubAssign,
    },
    str::FromStr,
    vec::Vec,
};
use crypto_bigint::Random;
use num_bigint::{ToBigInt, ToBigUint};
use num_traits::{ConstOne, ConstZero, One, Zero};

use crate::traits::{FromBytes, FromRef, MapsToField, ToBytes};

/// Trait for integer types used as field element representations.
pub trait BigInteger:
    From<u64>
    + From<u32>
    + Debug
    + FromBytes
    + FromStr
    + Clone
    + PartialEq
    + Eq
    + PartialOrd
    + UniformRand
    + Default
    + From<u8>
    + From<u16>
    + From<u32>
    + From<u64>
    + From<u128>
    + From<usize> // + From<i8>
// + From<i16>
// + From<i32>
{
    type W: Words;
    /// Converts the integer to its word representation.
    fn to_words(&self) -> Self::W;
    /// Returns the integer one.
    fn one() -> Self;
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
    fn add_with_carry(&mut self, other: &Self) -> bool;
    fn sub_with_borrow(&mut self, other: &Self) -> bool;
    fn mul(&self, other: &Self) -> (Self, Self);
    /// Constructs from big-endian bits.
    fn from_bits_be(bits: &[bool]) -> Self;
    /// Constructs from little-endian bits.
    fn from_bits_le(bits: &[bool]) -> Self;
    /// Returns the number of bits.
    fn num_bits(&self) -> u32;
    /// Constructs from words.
    fn new(words: Self::W) -> Self;
    /// Converts to big-endian bytes.
    fn to_bytes_be(self) -> Vec<u8>;
    /// Converts to little-endian bytes.
    fn to_bytes_le(self) -> Vec<u8>;
    fn demontgomery(&self, modulus: &Self, inv: u64) -> Self;
    fn last_mut(&mut self) -> &mut u64;
}

/// Trait for field configuration types.
pub trait Config: PartialEq + Eq + Clone {
    type B: BigInteger;
    /// Returns the modulus for the field.
    fn modulus(&self) -> &Self::B;
    /// Multiplies two integers in the field.
    fn mul_assign(&self, a: &mut Self::B, b: &Self::B);
    /// Returns the R^2 value for Montgomery reduction.
    fn r2(&self) -> &Self::B;
    /// Constructs a new config from a modulus.
    fn new(modulus: Self::B) -> Self;

    fn add_assign(&self, a: &mut Self::B, b: &Self::B);
    fn sub_assign(&self, a: &mut Self::B, b: &Self::B);

    fn reduce_modulus(&self, a: &mut Self::B, carry: bool);

    fn inverse(&self, a: &Self::B) -> Option<Self::B>;
    fn r(&self) -> &Self::B;
    fn inv(&self) -> u64;
}

/// Trait for references to field configuration.
pub trait ConfigReference: Copy + Clone + PartialEq + Eq + Debug + Send + Sync {
    type C: Config<B = Self::B>;
    type B: BigInteger<W = Self::W> + From<Self::I> + ToBigInt + ToBigUint + MapsToField<Self>;
    type I: Integer<W = Self::W, Uint = Self::U> + FromRef<Self::B>;
    type U: Uinteger<W = Self::W, Int = Self::I>;
    type W: Words;
    const N: usize;
    /// Returns a reference to the config, if available.
    fn reference(&self) -> &Self::C;
}

/// Trait for word-based representations of integers.
pub trait Words:
    Default
    + Index<usize, Output = Self::Word>
    + IndexMut<usize>
    + Index<Range<usize>, Output = [Self::Word]>
    + IndexMut<Range<usize>, Output = [Self::Word]>
    + Index<RangeTo<usize>, Output = [Self::Word]>
    + IndexMut<RangeTo<usize>, Output = [Self::Word]>
    + Clone
{
    type Word: PrimitiveConversions + Sized;
    /// Returns the number of words.
    fn num_words() -> usize;
}

/// Trait for cryptographic integer types.
pub trait Integer:
    Debug
    + Clone
    + PartialEq
    + Eq
    + PartialOrd
    + Default
    + Send
    + Sync
    + Zero
    + crypto_bigint::Zero
    + One
    + ConstZero
    + ConstOne
    + Neg<Output = Self>
    + Add<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + Mul<Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + AddAssign<Self>
    + for<'a> AddAssign<&'a Self>
    + RemAssign<Self>
    + Sum
    + FromRef<Self>
    + From<i64>
    + From<i32>
    + From<i8>
    + From<u8>
    + Random
    + ToBytes
{
    type W: Words;
    type Uint: Uinteger<W = Self::W>;
    type I: BigInteger<W = Self::W> + FromRef<Self>;

    fn num_bits() -> usize {
        Self::W::num_words() * <Self::W as Words>::Word::bits()
    }

    /// Constructs from words.
    fn from_words(words: Self::W) -> Self;
    fn as_words(&self) -> &[u64];
    fn from_i64(value: i64) -> Self;
    fn abs(&self) -> Self::Uint;
}

/// Trait for cryptographic unsigned integer types.
pub trait Uinteger: Clone + FromBytes + One + for<'a> SubAssign<&'a Self> {
    type W: Words;
    type Int: Integer<W = Self::W>;
    type PrimalityTest: PrimalityTest<Self>;
    /// Constructs from words.
    fn from_words(words: Self::W) -> Self;
    /// Converts to the signed integer type.
    fn as_int(&self) -> Self::Int;
    /// Converts to words.
    fn to_words(self) -> Self::W;
    fn is_even(&self) -> bool;
}

pub trait ZipTypes: Send + Sync {
    /// Width of elements in witness/polynomial evaluations on hypercube
    type N: Integer;

    /// Width of elements in the encoding matrices
    type L: Integer + FromRef<Self::N>;

    /// Width of elements in the code
    type K: Integer + FromRef<Self::N> + FromRef<Self::L>;

    /// Width of elements in linear combination of code rows
    type M: Integer + FromRef<Self::N> + FromRef<Self::L> + FromRef<Self::K>;
}

pub trait PrimalityTest<U: Uinteger> {
    type Inner;
    fn new(candidate: U) -> Self;
    fn is_probably_prime(&self) -> bool;
}

pub trait PrimitiveConversion<T> {
    fn from_primitive(value: T) -> Self;
}

pub trait InSameField {
    fn is_in_same_field(&self, other: &Self) -> bool;
}

macro_rules! impl_single_primitive_conversion {
    ($o:ty, $i:ty) => {
        impl PrimitiveConversion<$i> for $o {
            #[inline(always)]
            fn from_primitive(value: $i) -> $o {
                value as $o
            }
        }
    };
}

macro_rules! impl_primitive_conversion {
    ($o:ty) => {
        impl_single_primitive_conversion!($o, u8);
        impl_single_primitive_conversion!($o, u16);
        impl_single_primitive_conversion!($o, u32);
        impl_single_primitive_conversion!($o, u64);
        impl_single_primitive_conversion!($o, u128);
        impl_single_primitive_conversion!($o, usize);
        impl_single_primitive_conversion!($o, i8);
        impl_single_primitive_conversion!($o, i16);
        impl_single_primitive_conversion!($o, i32);
        impl_single_primitive_conversion!($o, i64);
        impl_single_primitive_conversion!($o, i128);
        impl_single_primitive_conversion!($o, isize);
    };
}

impl_primitive_conversion!(u8);
impl_primitive_conversion!(u16);
impl_primitive_conversion!(u32);
impl_primitive_conversion!(u64);
impl_primitive_conversion!(u128);
impl_primitive_conversion!(usize);
impl_primitive_conversion!(i8);
impl_primitive_conversion!(i16);
impl_primitive_conversion!(i32);
impl_primitive_conversion!(i64);
impl_primitive_conversion!(i128);
impl_primitive_conversion!(isize);

pub trait PrimitiveConversions:
    PrimitiveConversion<u8>
    + PrimitiveConversion<u16>
    + PrimitiveConversion<u32>
    + PrimitiveConversion<u64>
    + PrimitiveConversion<u128>
    + PrimitiveConversion<usize>
    + PrimitiveConversion<i8>
    + PrimitiveConversion<i16>
    + PrimitiveConversion<i32>
    + PrimitiveConversion<i64>
    + PrimitiveConversion<i128>
    + PrimitiveConversion<isize>
    + Shr<usize, Output = Self>
{
    fn bits() -> usize;
}
impl<T> PrimitiveConversions for T
where
    T: PrimitiveConversion<u8>
        + PrimitiveConversion<u16>
        + PrimitiveConversion<u32>
        + PrimitiveConversion<u64>
        + PrimitiveConversion<u128>
        + PrimitiveConversion<usize>
        + PrimitiveConversion<i8>
        + PrimitiveConversion<i16>
        + PrimitiveConversion<i32>
        + PrimitiveConversion<i64>
        + PrimitiveConversion<i128>
        + PrimitiveConversion<isize>
        + Shr<usize, Output = Self>,
{
    #[inline(always)]
    fn bits() -> usize {
        size_of::<Self>() * 8
    }
}
