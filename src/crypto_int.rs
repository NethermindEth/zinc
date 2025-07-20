use ark_std::{
    cmp::Ordering,
    fmt::Debug,
    ops::{Add, AddAssign, Mul, RemAssign, Sub},
    rand::RngCore,
    vec::Vec,
};
use crypto_bigint::{
    subtle::{Choice, ConstantTimeEq},
    Int, NonZero, Random, Uint,
};
use num_traits::{ConstOne, ConstZero, One, Zero};

use crate::{
    biginteger::{BigInt, Words},
    traits::{CryptoInteger, CryptoUinteger, FromBytes},
    zip::pcs::utils::ToBytes,
};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct CryptoInt<const N: usize>(Int<N>);

impl<const N: usize> Zero for CryptoInt<N> {
    #[inline]
    fn zero() -> Self {
        Self(Int::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<const N: usize> Add<Self> for CryptoInt<N> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<const N: usize> ConstZero for CryptoInt<N> {
    const ZERO: Self = Self(Int::ZERO);
}

impl<const N: usize> One for CryptoInt<N> {
    #[inline]
    fn one() -> Self {
        Self(Int::one())
    }
}

impl<const N: usize> ConstOne for CryptoInt<N> {
    const ONE: Self = Self(Int::ONE);
}

impl<const N: usize> Mul<Self> for CryptoInt<N> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl<const N: usize> ConstantTimeEq for CryptoInt<N> {
    #[inline]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl<const N: usize> PartialOrd for CryptoInt<N> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<const N: usize> RemAssign<Self> for CryptoInt<N> {
    #[inline]
    fn rem_assign(&mut self, rhs: Self) {
        self.0.rem_assign(
            NonZero::new(rhs.0).expect("Cannot reduce modulo zero: field modulus is zero"),
        );
    }
}

impl<'a, const N: usize> Add<&'a Self> for CryptoInt<N> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<'a, const N: usize> Mul<&'a Self> for CryptoInt<N> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl<'a, const N: usize> AddAssign<&'a Self> for CryptoInt<N> {
    #[inline]
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0 += &rhs.0;
    }
}

impl<'a, const N: usize> Sub<&'a Self> for CryptoInt<N> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<const N: usize> From<i64> for CryptoInt<N> {
    #[inline]
    fn from(value: i64) -> Self {
        Self(Int::from(value))
    }
}
impl<const N: usize> From<i32> for CryptoInt<N> {
    #[inline]
    fn from(value: i32) -> Self {
        Self(Int::from(value))
    }
}
impl<const N: usize> From<i8> for CryptoInt<N> {
    #[inline]
    fn from(value: i8) -> Self {
        Self(Int::from(value))
    }
}
impl<const N: usize> From<u8> for CryptoInt<N> {
    #[inline]
    fn from(value: u8) -> Self {
        Self(Int::from(value as i8))
    }
}

impl<const N: usize> Default for CryptoInt<N> {
    #[inline]
    fn default() -> Self {
        Self(Int::default())
    }
}

impl<const N: usize> Random for CryptoInt<N> {
    #[inline]
    fn random(rng: &mut (impl RngCore + ?Sized)) -> Self {
        Self(Int::random(rng))
    }
}

impl<'a, const N: usize, const M: usize> From<&'a CryptoInt<M>> for CryptoInt<N> {
    #[inline]
    fn from(value: &'a CryptoInt<M>) -> Self {
        Self(Int::from(&value.0))
    }
}

impl<const N: usize> ToBytes for CryptoInt<N> {
    // Manual impl for generic type
    fn to_bytes(&self) -> Vec<u8> {
        self.0
            .to_words()
            .iter()
            .flat_map(|word| word.to_be_bytes())
            .collect()
    }
}

impl<const N: usize> CryptoInteger for CryptoInt<N> {
    type W = crate::biginteger::Words<N>;
    type Uint = Uint<N>;
    type I = BigInt<N>;

    fn from_words(words: Words<N>) -> Self {
        Self(Int::from_words(words.0))
    }

    fn as_words(&self) -> &[u64] {
        self.0.as_words()
    }

    fn from_i64(value: i64) -> Self {
        Self(Int::from_i64(value))
    }

    fn abs(&self) -> Self::Uint {
        self.0.abs()
    }
}

impl<const N: usize> CryptoUinteger for Uint<N> {
    type W = crate::biginteger::Words<N>;
    type Int = CryptoInt<N>;

    fn from_words(words: Words<N>) -> Self {
        Self::from_words(words.0)
    }

    fn as_int(&self) -> Self::Int {
        CryptoInt(self.as_int())
    }

    fn to_words(self) -> Words<N> {
        Words(self.to_words())
    }
}

impl<const N: usize> FromBytes for Uint<N> {
    fn from_bytes_le(bytes: &[u8]) -> Option<Self> {
        Some(Self::from_le_slice(bytes))
    }

    fn from_bytes_be(bytes: &[u8]) -> Option<Self> {
        Some(Self::from_be_slice(bytes))
    }
}
