use ark_std::{
    fmt::Debug,
    ops::{
        AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Range, RangeTo, RemAssign, Sub, SubAssign,
    },
    UniformRand,
};
use crypto_bigint::{NonZero, Random};
use num_traits::{One, Zero};

use crate::traits::FieldMap;

pub trait Field:
    Zero
    + Eq
    + PartialEq
    + One
    + Clone
    + Copy
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
{
    type I: Integer<Self::W> + From<Self::CryptoInt> + FieldMap<Self, Output = Self>;
    type C: Config<Self::W, Self::I>;
    type Cr: ConfigReference<Self::W, Self::I, Self::C>;
    type W: Words;
    type CryptoInt: CryptoInt<Self::W, Uint = Self::CryptoUint> + for<'a> From<&'a Self::I>;
    type CryptoUint: CryptoUint<Self::W, Int = Self::CryptoInt>;
    fn new_unchecked(config: Self::Cr, value: Self::I) -> Self;
    fn without_config(value: Self::I) -> Self;
}

pub trait Integer<W: Words>: From<u64> + From<u32> + Debug {
    fn to_words(&self) -> W;

    fn one() -> Self;
}
pub trait Config<W: Words, I: Integer<W>> {
    fn modulus(&self) -> &I;
    fn mul_assign(&self, a: &mut I, b: &I);

    fn r2(&self) -> &I;
}
pub trait ConfigReference<W: Words, I: Integer<W>, C: Config<W, I>>:
    Copy + Clone + PartialEq + Eq + Debug
{
    fn reference(&self) -> Option<&C>;

    #[allow(clippy::missing_safety_doc)] // TODO Should be documented.
    unsafe fn new(config_ptr: *mut C) -> Self;

    fn pointer(&self) -> Option<*mut C>;
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

pub trait CryptoInt<W: Words>: crypto_bigint::Zero + RemAssign<NonZero<Self>> {
    type Uint: CryptoUint<W>;
    fn from_words(words: W) -> Self;
}
pub trait CryptoUint<W: Words>: crypto_bigint::Zero + RemAssign<NonZero<Self>> {
    type Int: CryptoInt<W>;
    fn from_words(words: W) -> Self;
    fn as_int(&self) -> Self::Int;
}
