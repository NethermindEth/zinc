use ark_std::vec::Vec;
use crypto_bigint::{Int, Uint, Word, modular::ConstMontyParams};
use num_traits::Zero;
use sha3::{Digest, Keccak256};

use crate::traits::{Field, FromBits, Transcript};

/// A cryptographic transcript implementation using the Keccak-256 hash function.
/// Used for Fiat-Shamir transformations in zero-knowledge proof systems.
#[derive(Clone)]
pub struct KeccakTranscript {
    /// The underlying Keccak-256 hasher that maintains the transcript state.
    hasher: Keccak256,
}

impl Default for KeccakTranscript {
    fn default() -> Self {
        Self::new()
    }
}

impl KeccakTranscript {
    pub fn new() -> Self {
        Self {
            hasher: Keccak256::new(),
        }
    }

    /// Absorbs arbitrary bytes into the transcript.
    /// This updates the internal state of the hasher with the provided data.
    pub fn absorb(&mut self, v: &[u8]) {
        self.hasher.update(v);
    }

    /// Generates a specified number of pseudorandom bytes based on the current transcript state.
    /// Uses a counter-based approach to generate enough bytes from the hasher.
    pub fn get_random_bytes(&mut self, length: usize) -> Vec<u8> {
        let mut result = Vec::with_capacity(length);
        let mut counter = 0;
        while result.len() < length {
            let mut temp_hasher = self.hasher.clone();
            temp_hasher.update(i32::to_be_bytes(counter));
            let hash = temp_hasher.finalize();
            result.extend_from_slice(&hash);

            counter += 1;
        }

        result.truncate(length);
        result
    }

    /// Absorbs a field element into the transcript.
    /// Delegates to the field element's implementation of absorb_into_transcript.
    pub fn absorb_random_field<F: Field<LIMBS>, const LIMBS: usize>(&mut self, v: &F) {
        //       self.absorb(&[0x3]);
        //       self.absorb(&ToBytes::to_be_bytes(MOD::PARAMS.modulus().as_ref()));
        //       self.absorb(&[0x5]);

        self.absorb(&[0x1]);
        self.absorb(&v.to_be_bytes());
        self.absorb(&[0x3])
    }

    /// Absorbs a slice of field elements into the transcript.
    /// Processes each field element in the slice sequentially.
    pub fn absorb_slice<F: Field<LIMBS>, const LIMBS: usize>(&mut self, slice: &[F]) {
        for field_element in slice.iter() {
            self.absorb_random_field(field_element);
        }
    }

    /// Internal helper that generates two 128-bit limbs from the current transcript state.
    /// Updates the transcript state.
    fn get_challenge_limbs(&mut self) -> (u128, u128) {
        let challenge = self.hasher.clone().finalize();

        // Interpret the digest as big-endian halves but keep the original ordering used by this protocol:
        // lo = first 16 bytes, hi = last 16 bytes.
        let lo = u128::from_be_bytes(challenge[0..16].try_into().unwrap());
        let hi = u128::from_be_bytes(challenge[16..32].try_into().unwrap());

        self.hasher.update([0x00]);
        self.hasher.update(challenge);
        self.hasher.update([0x01]);

        (lo, hi)
    }

    /// Generates a pseudorandom field element as a challenge based on the current transcript state.
    #[allow(clippy::not_unsafe_ptr_arg_deref)]
    pub fn get_challenge<F: Field<LIMBS>, const LIMBS: usize>(&mut self) -> F {
        let (lo, hi) = self.get_challenge_limbs();
        let modulus = F::Monty::PARAMS.modulus().as_ref();
        let challenge_num_bits = modulus.bits() - 1;
        if LIMBS == 1 {
            let lo_mask = (1u64 << challenge_num_bits) - 1;

            let truncated_lo = lo as u64 & lo_mask;

            let challenge = truncated_lo.into();
            return challenge;
        }
        if challenge_num_bits < 128 {
            let lo_mask = (1u128 << challenge_num_bits) - 1;

            let truncated_lo = lo & lo_mask;

            truncated_lo.into()
        } else if challenge_num_bits >= 256 {
            let two_to_128 = F::from(Uint::<LIMBS>::from_le_bits(
                &(0..196).map(|i| i == 128).collect::<Vec<bool>>(),
            ));

            F::from(lo) + two_to_128 * F::from(hi)
        } else {
            let hi_bits_to_keep = challenge_num_bits - 128;
            let hi_mask = (1u128 << hi_bits_to_keep) - 1;

            let truncated_hi = hi & hi_mask;

            let two_to_128 = F::from(Uint::<LIMBS>::from_le_bits(
                &(0..196).map(|i| i == 128).collect::<Vec<bool>>(),
            ));

            F::from(lo) + two_to_128 * F::from(truncated_hi)
        }
    }

    /// Generates pseudorandom field elements as challenges based on the current transcript state.
    pub fn get_challenges<F: Field<LIMBS>, const LIMBS: usize>(&mut self, n: usize) -> Vec<F> {
        let mut challenges = Vec::with_capacity(n);
        challenges.extend((0..n).map(|_| self.get_challenge::<F, LIMBS>()));
        challenges
    }

    /// Generates a pseudorandom [Integer] as a challenge based on the current transcript state.
    pub fn get_integer_challenge<const LIMBS: usize>(&mut self) -> Int<LIMBS> {
        let mut words: [Word; LIMBS] = [Word::zero(); LIMBS];
        for word in words.iter_mut().take(LIMBS) {
            let mut challenge = [0u8; size_of::<Word>()];
            let rand_bytes = self.get_random_bytes(size_of::<Word>());
            challenge.copy_from_slice(&rand_bytes);
            self.hasher.update([0x12]);
            self.hasher.update(challenge);
            self.hasher.update([0x34]);
            *word = Word::from_le_bytes(challenge);
        }

        Int::from_words(words)
    }

    /// Generates pseudorandom [CryptoInt]s as challenges based on the current transcript state.
    pub fn get_integer_challenges<const LIMBS: usize>(&mut self, n: usize) -> Vec<Int<LIMBS>> {
        (0..n).map(|_| self.get_integer_challenge()).collect()
    }

    /// Generates a pseudorandom `usize` within the given range bounds based on the current transcript state.
    fn get_usize_in_range(&mut self, range: &ark_std::ops::Range<usize>) -> usize {
        let challenge = self.hasher.clone().finalize();

        self.hasher.update([0x88]);
        self.hasher.update(challenge);
        self.hasher.update([0x11]);

        let num = usize::from_le_bytes(challenge[..size_of::<usize>()].try_into().unwrap());
        range.start + (num % (range.end - range.start))
    }
}

impl<const LIMBS: usize> Transcript<LIMBS> for KeccakTranscript {
    fn get_encoding_element(&mut self) -> Int<LIMBS> {
        let byte = self.get_random_bytes(1)[0];
        // cancels all bits and depends only on whether the random byte LSB is 0 or 1
        let bit = byte & 1;
        Int::from(bit as i8)
    }

    fn get_word(&mut self) -> Word {
        self.get_integer_challenge::<1>().as_words()[0]
    }

    fn sample_unique_columns(
        &mut self,
        range: ark_std::ops::Range<usize>,
        columns: &mut ark_std::collections::BTreeSet<usize>,
        count: usize,
    ) -> usize {
        let mut added = 0;
        while added < count {
            let candidate = self.get_usize_in_range(&range);
            if columns.insert(candidate) {
                added += 1;
            }
        }
        added
    }
}
#[cfg(test)]
mod tests {
    use crypto_bigint::{U256, Uint, const_monty_params};

    use super::KeccakTranscript;
    use crate::field::{F256, RandomField, WORD_FACTOR};

    const_monty_params!(
        ModP,
        U256,
        "0800000000000011000000000000000000000000000000000000000000000001"
    );

    type F = F256<ModP>;

    #[test]
    fn test_keccak_transcript() {
        let mut transcript = KeccakTranscript::new();

        transcript.absorb(b"This is a test string!");
        let challenge = transcript.get_challenge::<F, { 4 * WORD_FACTOR }>();

        let expected =
            Uint::from_be_hex("018841C8CCF58149C2077504AEFE70BB2935E37F8168922CCFEF7DEE720FEFBE");
        let expected = RandomField::from(expected);

        assert_eq!(challenge, expected);
    }
}
