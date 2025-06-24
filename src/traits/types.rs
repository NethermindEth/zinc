use ark_std::{
    fmt::Debug,
    ops::{AddAssign, Index, IndexMut, Mul, Neg, Range, RangeTo, RemAssign},
    UniformRand,
};
use crypto_bigint::NonZero;
use num_traits::{One, Zero};

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
{
    type I: Integer<Self::W> + From<Self::CryptoInt>;
    type C: Config<Self::W, Self::I>;
    type Cr: ConfigReference<Self::W, Self::I, Self::C>;
    type W: Words;
    type CryptoInt: CryptoInt<Self::W, Uint = Self::CryptoUint> + for<'a> From<&'a Self::I>;
    type CryptoUint: CryptoUint<Self::W, Int = Self::CryptoInt>;
    fn new_unchecked(config: Self::Cr, value: Self::I) -> Self;
    fn without_config(value: Self::I) -> Self;
}

pub trait Integer<W: Words>: From<u64> + Debug {
    fn to_words(&self) -> W;
}
pub trait Config<W: Words, I: Integer<W>> {
    fn modulus(&self) -> &I;
    fn mul_assign(&self, a: &mut I, b: &I);

    fn r2(&self) -> &I;
}
pub trait ConfigReference<W: Words, I: Integer<W>, C: Config<W, I>>: Copy + Clone {
    fn reference(&self) -> Option<&C>;

    #[allow(clippy::missing_safety_doc)] // TODO Should be documented.
    unsafe fn new(config_ptr: *mut C) -> Self;

    fn pointer(&self) -> Option<*mut C>;
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
