use ark_std::ops::{Mul, SubAssign};
use crypto_bigint::{Integer, Odd, Uint};
use crypto_primes::hazmat::MillerRabin;
use num_traits::One;

use crate::{
    biginteger::Words,
    crypto_int::CryptoInt,
    traits::{types::PrimalityTest, CryptoUinteger, FromBytes},
};

#[derive(Clone)]
pub struct CryptoUint<const N: usize>(pub(crate) Uint<N>);

impl<const N: usize> Mul<Self> for CryptoUint<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl<const N: usize> One for CryptoUint<N> {
    fn one() -> Self {
        Self(Uint::ONE)
    }
}

impl<'a, const N: usize> SubAssign<&'a Self> for CryptoUint<N> {
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0 -= rhs.0;
    }
}

impl<const N: usize> CryptoUinteger for CryptoUint<N> {
    type W = Words<N>;
    type Int = CryptoInt<N>;
    type PrimalityTest = MillerRabin<Uint<N>>;

    fn from_words(words: Words<N>) -> Self {
        Self(Uint::from_words(words.0))
    }

    fn as_int(&self) -> Self::Int {
        CryptoInt(self.0.as_int())
    }

    fn to_words(self) -> Words<N> {
        Words(self.0.to_words())
    }

    fn is_even(&self) -> bool {
        self.0.is_even().unwrap_u8() == 1
    }
}

impl<const N: usize> FromBytes for CryptoUint<N> {
    fn from_bytes_le(bytes: &[u8]) -> Option<Self> {
        Some(Self(Uint::from_le_slice(bytes)))
    }

    fn from_bytes_be(bytes: &[u8]) -> Option<Self> {
        Some(Self(Uint::from_be_slice(bytes)))
    }
}

impl<const N: usize> PrimalityTest<CryptoUint<N>> for MillerRabin<Uint<N>> {
    type Inner = Uint<N>;

    fn new(candidate: CryptoUint<N>) -> Self {
        Self::new(Odd::new(candidate.0).unwrap())
    }

    fn is_probably_prime(&self) -> bool {
        self.test_base_two().is_probably_prime()
    }
}

// pub(crate) struct MillerRabin<U: CryptoUinteger>(<Self as PrimalityTest<U>>::Inner) where MillerRabin<U>: PrimalityTest<U>;
//
// impl<const N: usize> PrimalityTest for MillerRabin<U> {
//     type Inner = crypto_primes::hazmat::MillerRabin<U::GenericType>;
//     fn new(candidate: U) -> Self {
//         Self(
//             crypto_primes::hazmat::MillerRabin::new(Odd::new(candidate.into_generic_type()).unwrap()),
//         )
//     }
//
//     fn is_probably_prime(&self) -> bool {
//         self.0.test_base_two().is_probably_prime()
//     }
// }
