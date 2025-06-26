use ark_std::{
    fmt::Debug,
    ops::{
        AddAssign, Div, Index, IndexMut, Mul, MulAssign, Neg, Range, RangeTo, RemAssign, Sub,
        SubAssign,
    },
    vec::Vec,
    UniformRand,
};
use crypto_bigint::{NonZero, Random, RandomMod};
use num_traits::{One, Zero};

use crate::{
    traits::{FieldMap, FromBytes},
    transcript::KeccakTranscript,
};

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
    + UniformRand
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
    type I: Integer<W = Self::W> + From<Self::CryptoInt> + FieldMap<Self, Output = Self>;
    type C: Config<I = Self::I>;
    type Cr: ConfigReference<C = Self::C>;
    type W: Words;
    type CryptoInt: CryptoInt<W = Self::W, Uint = Self::CryptoUint> + for<'a> From<&'a Self::I>;
    type CryptoUint: CryptoUint<W = Self::W, Int = Self::CryptoInt>;
    type DebugField: Debug + From<Self> + Send + Sync;
    fn new_unchecked(config: Self::Cr, value: Self::I) -> Self;
    fn without_config(value: Self::I) -> Self;
    fn rand_with_config<R: ark_std::rand::Rng + ?Sized>(rng: &mut R, config: Self::Cr) -> Self;
    fn set_config(&mut self, config: Self::Cr);
    fn value(&self) -> &Self::I;
    fn value_mut(&mut self) -> &mut Self::I;
    fn absorb_into_transcript(&self, transcript: &mut KeccakTranscript);
}

pub trait Integer: From<u64> + From<u32> + Debug + FromBytes + Copy {
    type W: Words;
    fn to_words(&self) -> Self::W;

    fn one() -> Self;
    fn from_bits_be(bits: &[bool]) -> Self;

    fn from_bits_le(bits: &[bool]) -> Self;
    fn num_bits(&self) -> u32;
    fn new(words: Self::W) -> Self;
    fn to_bytes_be(self) -> Vec<u8>;

    fn to_bytes_le(self) -> Vec<u8>;
}
pub trait Config: PartialEq + Eq {
    type I: Integer;
    fn modulus(&self) -> &Self::I;
    fn mul_assign(&self, a: &mut Self::I, b: &Self::I);

    fn r2(&self) -> &Self::I;
    fn new(modulus: Self::I) -> Self;
}
pub trait ConfigReference: Copy + Clone + PartialEq + Eq + Debug + Send + Sync {
    type C: Config;
    fn reference(&self) -> Option<&Self::C>;

    #[allow(clippy::missing_safety_doc)] // TODO Should be documented.
    unsafe fn new(config_ptr: *mut Self::C) -> Self;

    fn pointer(&self) -> Option<*mut Self::C>;
    const NONE: Self;
}

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
    fn num_words() -> usize;
}

pub trait CryptoInt: crypto_bigint::Zero + PartialOrd + RemAssign<NonZero<Self>> {
    type W: Words;
    type Uint: CryptoUint<W = Self::W>;
    type I: Integer<W = Self::W> + for<'a> From<&'a Self>;
    fn from_words(words: Self::W) -> Self;
}
pub trait CryptoUint:
    crypto_bigint::Integer + RemAssign<NonZero<Self>> + FromBytes + RandomMod + Copy
{
    type W: Words;
    type Int: CryptoInt<W = Self::W>;
    fn from_words(words: Self::W) -> Self;
    fn as_int(&self) -> Self::Int;
    fn to_words(self) -> Self::W;
}
