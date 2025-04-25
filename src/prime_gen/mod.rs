use crate::{biginteger::BigInt, transcript::KeccakTranscript};
use crypto_bigint::{Integer, Odd, Uint};
use crypto_primes::hazmat::MillerRabin;

fn hash_int<const N: usize>(hasher: &mut KeccakTranscript) -> Uint<N> {
    let n_bytes = N * 8;
    let bytes = hasher.get_random_bytes(n_bytes);
    hasher.absorb(&bytes);
    Uint::<N>::from_be_slice(&bytes)
}

pub fn get_prime<const N: usize>(hasher: &mut KeccakTranscript) -> BigInt<N> {
    let prime = loop {
        let mut prime_candidate: Uint<N> = hash_int(hasher);
        if prime_candidate.is_even().unwrap_u8() == 1 {
            prime_candidate -= Uint::<N>::one();
        }
        let mr = MillerRabin::new(Odd::new(prime_candidate).unwrap());

        if mr.test_base_two().is_probably_prime() {
            break prime_candidate;
        }
    };
    BigInt::new(prime.to_words())
}

#[cfg(test)]
mod test {
    use crate::{prime_gen::get_prime, transcript::KeccakTranscript};

    #[test]
    fn test_prime_generator() {
        let mut hasher = KeccakTranscript::new();
        const N: usize = 3;

        println!("{:?}", get_prime::<N>(&mut hasher))
    }
}
