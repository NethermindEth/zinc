use ark_std::{
    fmt::Debug,
    ops::{
        Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign, Neg, Range, RangeTo, RemAssign, Sub,
        SubAssign,
    },
    vec::Vec,
};
use crypto_bigint::Random;
use num_traits::{ConstOne, ConstZero, One, Zero};

use crate::{
    traits::{FieldMap, FromBytes},
    transcript::KeccakTranscript,
    zip::pcs::utils::ToBytes,
};

/// Trait for field elements, requiring arithmetic, assignment, random generation, and conversion traits.
/// Used as a bound for generic code over finite fields.
pub trait Field:
    Zero
    + Eq
    + PartialEq
    + One
    + Clone
    + Default
    + Sync
    + Send
    + for<'a> Mul<&'a Self, Output = Self>
    + Neg<Output = Self>
    + Debug
    + for<'a> AddAssign<&'a Self>
    + Sub<Self, Output = Self>
    + MulAssign
    + SubAssign
    + Random
    + AddAssign
    + Div<Self, Output = Self>
    + for<'a> core::iter::Product<&'a Self>
    + core::iter::Sum
    + for<'a> MulAssign<&'a Self>
{
    /// Integer representation type for the field element.
    type I: Integer<W = Self::W> + From<Self::CryptoInt> + FieldMap<Self, Output = Self>;
    /// Field configuration type.
    type C: Config<I = Self::I>;
    /// Reference to field configuration.
    type Cr: ConfigReference<C = Self::C>;
    /// Word representation type.
    type W: Words;
    /// Cryptographic integer type.
    type CryptoInt: CryptoInteger<W = Self::W, Uint = Self::CryptoUint> + for<'a> From<&'a Self::I>;
    /// Cryptographic unsigned integer type.
    type CryptoUint: CryptoUinteger<W = Self::W, Int = Self::CryptoInt>;
    /// Debug representation for the field.
    type DebugField: Debug + From<Self> + Send + Sync;
    /// Creates a new field element from config and value, without checking.
    fn new_unchecked(config: Self::Cr, value: Self::I) -> Self;
    /// Creates a new field element from value, without config.
    fn without_config(value: Self::I) -> Self;
    /// Generates a random field element with the given config.
    fn rand_with_config<R: ark_std::rand::Rng + ?Sized>(rng: &mut R, config: Self::Cr) -> Self;
    /// Sets the field configuration.
    fn set_config(&mut self, config: Self::Cr);
    /// Returns a reference to the integer value.
    fn value(&self) -> &Self::I;
    /// Returns a mutable reference to the integer value.
    fn value_mut(&mut self) -> &mut Self::I;
    /// Absorbs the field element into a Keccak transcript.
    fn absorb_into_transcript(&self, transcript: &mut KeccakTranscript);
}

/// Trait for integer types used as field element representations.
pub trait Integer: From<u64> + From<u32> + Debug + FromBytes + Copy {
    type W: Words;
    /// Converts the integer to its word representation.
    fn to_words(&self) -> Self::W;
    /// Returns the integer one.
    fn one() -> Self;
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
}

/// Trait for field configuration types.
pub trait Config: PartialEq + Eq {
    type I: Integer;
    /// Returns the modulus for the field.
    fn modulus(&self) -> &Self::I;
    /// Multiplies two integers in the field.
    fn mul_assign(&self, a: &mut Self::I, b: &Self::I);
    /// Returns the R^2 value for Montgomery reduction.
    fn r2(&self) -> &Self::I;
    /// Constructs a new config from a modulus.
    fn new(modulus: Self::I) -> Self;
}

/// Trait for references to field configuration.
pub trait ConfigReference: Copy + Clone + PartialEq + Eq + Debug + Send + Sync {
    type C: Config;
    /// Returns a reference to the config, if available.
    fn reference(&self) -> Option<&Self::C>;
    /// Creates a new config reference from a pointer.
    #[allow(clippy::missing_safety_doc)] // TODO Should be documented.
    unsafe fn new(config_ptr: *mut Self::C) -> Self;
    /// Returns a pointer to the config, if available.
    fn pointer(&self) -> Option<*mut Self::C>;
    /// Constant representing no config reference.
    const NONE: Self;
}

/// Trait for word-based representations of integers.
pub trait Words:
    Default
    + Index<usize, Output = u64>
    + IndexMut<usize>
    + Index<Range<usize>, Output = [u64]>
    + IndexMut<Range<usize>, Output = [u64]>
    + Index<RangeTo<usize>, Output = [u64]>
    + IndexMut<RangeTo<usize>, Output = [u64]>
    + Copy
{
    /// Returns the number of words.
    fn num_words() -> usize;
}

/// Trait for cryptographic integer types.
pub trait CryptoInteger:
    Zero
    + ConstZero
    + One
    + ConstOne
    + crypto_bigint::Zero
    + PartialOrd
    + PartialEq
    + Eq
    + RemAssign<Self>
    + Clone
    + Add<Output = Self>
    + Mul<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + From<i64>
    + From<i32>
    + From<i8>
    + From<u8>
    + Default
    + Random
    + Send
    + Sync
    + Debug
    + for<'a> From<&'a Self>
    + ToBytes
{
    type W: Words;
    type Uint: CryptoUinteger<W = Self::W>;
    type I: Integer<W = Self::W> + for<'a> From<&'a Self>;
    /// Constructs from words.
    fn from_words(words: Self::W) -> Self;
    fn as_words(&self) -> &[u64];
    fn from_i64(value: i64) -> Self;
    fn abs(&self) -> Self::Uint;
}

/// Trait for cryptographic unsigned integer types.
pub trait CryptoUinteger: Clone + FromBytes + One + for<'a> SubAssign<&'a Self> {
    type W: Words;
    type Int: CryptoInteger<W = Self::W>;
    type PrimalityTest: PrimalityTest<Self>;
    /// Constructs from words.
    fn from_words(words: Self::W) -> Self;
    /// Converts to the signed integer type.
    fn as_int(&self) -> Self::Int;
    /// Converts to words.
    fn to_words(self) -> Self::W;
    fn is_even(&self) -> bool;
}

pub trait PrimalityTest<U: CryptoUinteger> {
    type Inner;
    fn new(candidate: U) -> Self;
    fn is_probably_prime(&self) -> bool;
}
