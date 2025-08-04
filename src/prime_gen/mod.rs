use num_traits::One;

use crate::{
    traits::{types::PrimalityTest, BigInteger, Field, FromBytes, Uinteger, Words},
    transcript::KeccakTranscript,
};

fn hash_int<F: Field>(hasher: &mut KeccakTranscript) -> F::U {
    let n_bytes = F::W::num_words() * 8;
    let bytes = hasher.get_random_bytes(n_bytes);
    hasher.absorb(&bytes);
    F::U::from_bytes_be(&bytes).expect("Failed to convert bytes to CryptoUint")
}

pub fn get_prime<F: Field>(hasher: &mut KeccakTranscript) -> F::B {
    let prime = loop {
        let mut prime_candidate: F::U = hash_int::<F>(hasher);
        if prime_candidate.is_even() {
            prime_candidate -= &F::U::one();
        }
        let mr = <<F as Field>::U as Uinteger>::PrimalityTest::new(prime_candidate.clone());

        if mr.is_probably_prime() {
            break prime_candidate;
        }
    };
    F::B::new(prime.to_words())
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
