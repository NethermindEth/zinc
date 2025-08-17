use num_traits::One;

use crate::{
    traits::{BigInteger, ConfigReference, FromBytes, Uinteger, types::PrimalityTest},
    transcript::KeccakTranscript,
};

fn hash_int<C: ConfigReference>(hasher: &mut KeccakTranscript) -> C::U {
    let n_bytes = C::N * 8;
    let bytes = hasher.get_random_bytes(n_bytes);
    hasher.absorb(&bytes);
    C::U::from_bytes_be(&bytes).expect("Failed to convert bytes to CryptoUint")
}

pub fn get_prime<C: ConfigReference>(hasher: &mut KeccakTranscript) -> C::B {
    let prime = loop {
        let mut prime_candidate = hash_int::<C>(hasher);
        if prime_candidate.is_even() {
            prime_candidate -= &C::U::one();
        }
        let mr = <C::U as Uinteger>::PrimalityTest::new(prime_candidate.clone());

        if mr.is_probably_prime() {
            break prime_candidate;
        }
    };
    C::B::new(prime.to_words())
}

#[cfg(test)]
mod test {
    use crate::{field::ConfigRef, prime_gen::get_prime, transcript::KeccakTranscript};

    #[test]
    fn test_prime_generator() {
        let mut hasher = KeccakTranscript::new();
        const N: usize = 3;
        get_prime::<ConfigRef<N>>(&mut hasher);
    }
}
