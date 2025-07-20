use crypto_bigint::Uint;

use crate::{
    biginteger::Words,
    crypto_int::CryptoInt,
    traits::{CryptoUinteger, FromBytes},
};

impl<const N: usize> CryptoUinteger for Uint<N> {
    type W = Words<N>;
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
