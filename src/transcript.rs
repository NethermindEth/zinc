use ark_std::vec::Vec;
use sha3::{Digest, Keccak256};

use crate::{
    traits::{
        BigInteger, Config, ConfigReference, Field, FieldMap, Integer, PrimitiveConversion, Words,
    },
    zip::pcs::structs::ZipTranscript,
};

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
    pub fn absorb_random_field<F: Field>(&mut self, v: &F) {
        v.absorb_into_transcript(self)
    }

    /// Absorbs a slice of field elements into the transcript.
    /// Processes each field element in the slice sequentially.
    pub fn absorb_slice<F: Field>(&mut self, slice: &[F]) {
        for field_element in slice.iter() {
            self.absorb_random_field(field_element);
        }
    }

    /// Internal helper that generates two 128-bit limbs from the current transcript state.
    /// Updates the transcript state.
    fn get_challenge_limbs(&mut self) -> (u128, u128) {
        let challenge = self.hasher.clone().finalize();

        let lo = u128::from_be_bytes(challenge[0..16].try_into().unwrap());
        let hi = u128::from_be_bytes(challenge[16..32].try_into().unwrap());

        self.hasher.update([0x00]);
        self.hasher.update(challenge);
        self.hasher.update([0x01]);

        (lo, hi)
    }

    /// Generates a pseudorandom field element as a challenge based on the current transcript state.
    #[allow(clippy::not_unsafe_ptr_arg_deref)]
    pub fn get_challenge<F: Field>(&mut self, config_ref: F::R) -> F {
        let (lo, hi) = self.get_challenge_limbs();
        let config = config_ref.reference().expect("Field config cannot be none");
        let modulus = config.modulus();
        let challenge_num_bits = modulus.num_bits() - 1;
        if F::W::num_words() == 1 {
            let lo_mask = (1u64 << challenge_num_bits) - 1;

            let truncated_lo = lo as u64 & lo_mask;

            let mut challenge: F = truncated_lo.map_to_field(config_ref);
            challenge.set_config(config_ref);
            return challenge;
        }
        if challenge_num_bits < 128 {
            let lo_mask = (1u128 << challenge_num_bits) - 1;

            let truncated_lo = lo & lo_mask;

            let mut challenge: F = truncated_lo.map_to_field(config_ref);
            challenge.set_config(config_ref);
            challenge
        } else if challenge_num_bits >= 256 {
            let two_to_128 = F::B::from_bits_le(&(0..196).map(|i| i == 128).collect::<Vec<bool>>())
                .map_to_field(config_ref);

            let mut challenge: F = <u128 as FieldMap<F>>::map_to_field(&lo, config_ref)
                + two_to_128 * <u128 as FieldMap<F>>::map_to_field(&hi, config_ref);
            challenge.set_config(config_ref);
            challenge
        } else {
            let hi_bits_to_keep = challenge_num_bits - 128;
            let hi_mask = (1u128 << hi_bits_to_keep) - 1;

            let truncated_hi = hi & hi_mask;

            let two_to_128 = F::B::from_bits_le(&(0..196).map(|i| i == 128).collect::<Vec<bool>>())
                .map_to_field(config_ref);

            let mut ret: F = <u128 as FieldMap<F>>::map_to_field(&lo, config_ref)
                + two_to_128 * <u128 as FieldMap<F>>::map_to_field(&truncated_hi, config_ref);
            ret.set_config(config_ref);
            ret
        }
    }

    pub fn get_challenges<F: Field>(&mut self, n: usize, config: F::R) -> Vec<F> {
        let mut challenges = Vec::with_capacity(n);
        challenges.extend((0..n).map(|_| self.get_challenge::<F>(config)));
        challenges
    }

    pub fn get_integer_challenge<I: Integer>(&mut self) -> I {
        let mut words = I::W::default();

        for i in 0..I::W::num_words() {
            let mut challenge = [0u8; 8];
            challenge.copy_from_slice(self.get_random_bytes(8).as_slice());
            self.hasher.update([0x12]);
            self.hasher.update(challenge);
            self.hasher.update([0x34]);
            words[i] = PrimitiveConversion::from_primitive(u64::from_le_bytes(challenge));
        }

        I::from_words(words)
    }

    pub fn get_integer_challenges<I: Integer>(&mut self, n: usize) -> Vec<I> {
        (0..n).map(|_| self.get_integer_challenge()).collect()
    }

    fn get_usize_in_range(&mut self, range: &ark_std::ops::Range<usize>) -> usize {
        let challenge = self.hasher.clone().finalize();

        self.hasher.update([0x88]);
        self.hasher.update(challenge);
        self.hasher.update([0x11]);

        let num = usize::from_le_bytes(challenge[..8].try_into().unwrap());
        range.start + (num % (range.end - range.start))
    }
}

impl<I: Integer> ZipTranscript<I> for KeccakTranscript {
    fn get_encoding_element(&mut self) -> I {
        let byte = self.get_random_bytes(1)[0];
        // cancels all bits and depends only on whether the random byte LSB is 0 or 1
        let bit = byte & 1;
        I::from(bit as i8)
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
    use ark_std::str::FromStr;

    use super::KeccakTranscript;
    use crate::{
        field::{BigInt, ConfigRef, FieldConfig, RandomField},
        traits::{Config, FieldMap},
    };

    #[test]
    fn test_keccak_transcript() {
        let mut transcript = KeccakTranscript::new();
        let config = FieldConfig::new(
            BigInt::<32>::from_str(
                "3618502788666131213697322783095070105623107215331596699973092056135872020481",
            )
            .unwrap(),
        );
        let field_config = ConfigRef::from(&config);

        transcript.absorb(b"This is a test string!");
        let challenge: RandomField<32> = transcript.get_challenge(field_config);

        let expected_bigint = BigInt::<32>::from_str(
            "693058076479703886486101269644733982722902192016595549603371045888466087870",
        )
        .unwrap();
        let expected_field = expected_bigint.map_to_field(field_config);

        assert_eq!(challenge, expected_field);
    }
}
