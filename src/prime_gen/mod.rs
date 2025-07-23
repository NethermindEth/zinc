use crypto_bigint::{Integer as CryptoInteger, Odd};
use crypto_primes::hazmat::MillerRabin;

use crate::{
    traits::{CryptoUint, Field, FromBytes, Integer, Words},
    transcript::KeccakTranscript,
};

fn hash_int<F: Field>(hasher: &mut KeccakTranscript) -> F::CryptoUint {
    let n_bytes = F::W::num_words() * 8;
    let bytes = hasher.get_random_bytes(n_bytes);
    hasher.absorb(&bytes);
    F::CryptoUint::from_bytes_be(&bytes).expect("Failed to convert bytes to CryptoUint")
}

pub fn get_prime<F: Field>(hasher: &mut KeccakTranscript) -> F::I {
    let prime = loop {
        let mut prime_candidate: F::CryptoUint = hash_int::<F>(hasher);
        if prime_candidate.is_even().unwrap_u8() == 1 {
            prime_candidate -= F::CryptoUint::one();
        }
        let mr = MillerRabin::new(Odd::new(prime_candidate).unwrap());

        if mr.test_base_two().is_probably_prime() {
            break prime_candidate;
        }
    };
    F::I::new(prime.to_words())
}

#[cfg(test)]
mod test {
    use crate::{field::RandomField, prime_gen::get_prime, transcript::KeccakTranscript};

    #[test]
    fn test_prime_generator() {
        let mut hasher = KeccakTranscript::new();
        const N: usize = 3;
        get_prime::<RandomField<N>>(&mut hasher);
    }
}
