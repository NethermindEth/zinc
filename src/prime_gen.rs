use crypto_bigint::{Integer, Odd, Uint, Word};
use crypto_primes::hazmat::MillerRabin;

use crate::transcript::KeccakTranscript;

fn hash_int<const LIMBS: usize>(hasher: &mut KeccakTranscript) -> Uint<LIMBS> {
    let n_bytes = LIMBS * size_of::<Word>();
    let bytes = hasher.get_random_bytes(n_bytes);
    hasher.absorb(&bytes);

    Uint::from_be_slice(&bytes)
}

pub fn get_prime<const LIMBS: usize>(hasher: &mut KeccakTranscript) -> Uint<LIMBS> {
    loop {
        let mut prime_candidate = hash_int::<LIMBS>(hasher);
        if prime_candidate.is_even().into() {
            prime_candidate -= Uint::ONE;
        }

        let prime_test = MillerRabin::new(Odd::new(prime_candidate).unwrap());

        let p = prime_test.test_base_two();
        if p.is_probably_prime() {
            break prime_candidate;
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{field::WORD_FACTOR, prime_gen::get_prime, transcript::KeccakTranscript};

    #[test]
    fn test_prime_generator() {
        let mut hasher = KeccakTranscript::new();
        const N: usize = 3 * WORD_FACTOR;
        get_prime::<N>(&mut hasher);
    }
}
