use crypto_bigint::{Int, Uint};

use crate::{
    biginteger::Words,
    traits::{CryptoInt, CryptoUint},
};

impl<const N: usize> CryptoInt<crate::biginteger::Words<N>> for Int<N> {
    type Uint = Uint<N>;

    fn from_words(words: Words<N>) -> Self {
        Self::from_words(words.0)
    }
}
impl<const N: usize> CryptoUint<crate::biginteger::Words<N>> for Uint<N> {
    type Int = Int<N>;

    fn from_words(words: Words<N>) -> Self {
        Self::from_words(words.0)
    }

    fn as_int(&self) -> Self::Int {
        self.as_int()
    }
}
