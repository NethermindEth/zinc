use crypto_bigint::{Int, Uint};

use crate::{
    biginteger::Words,
    traits::{CryptoInt, CryptoUint, FromBytes},
};

impl<const N: usize> CryptoInt for Int<N> {
    type W = crate::biginteger::Words<N>;
    type Uint = Uint<N>;

    fn from_words(words: Words<N>) -> Self {
        Self::from_words(words.0)
    }
}
impl<const N: usize> CryptoUint for Uint<N> {
    type W = crate::biginteger::Words<N>;
    type Int = Int<N>;

    fn from_words(words: Words<N>) -> Self {
        Self::from_words(words.0)
    }

    fn as_int(&self) -> Self::Int {
        self.as_int()
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
