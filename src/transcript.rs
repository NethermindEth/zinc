use ark_std::vec::Vec;
use crypto_bigint::Int;
use sha3::{Digest, Keccak256};

use crate::{
    biginteger::BigInt, field::RandomField, field_config::ConfigRef, traits::FieldMap,
    zip::pcs::structs::ZipTranscript,
};

#[derive(Clone)]
pub struct KeccakTranscript {
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

    pub fn absorb(&mut self, v: &[u8]) {
        self.hasher.update(v);
    }

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

    pub fn absorb_random_field<const N: usize>(&mut self, v: &RandomField<N>) {
        match v {
            RandomField::Raw { value } => {
                self.absorb(&[0x1]);
                self.absorb(&value.to_bytes_be());
                self.absorb(&[0x3])
            }
            RandomField::Initialized { config, value } => {
                let config = config.reference().expect("Field config cannot be none");

                self.absorb(&[0x3]);
                self.absorb(&config.modulus().to_bytes_be());
                self.absorb(&[0x5]);

                self.absorb(&[0x1]);
                self.absorb(&value.to_bytes_be());
                self.absorb(&[0x3])
            }
        }
    }

    pub fn absorb_slice<const N: usize>(&mut self, slice: &[RandomField<N>]) {
        for field_element in slice.iter() {
            self.absorb_random_field(field_element);
        }
    }

    fn get_challenge_limbs(&mut self) -> (u128, u128) {
        let challenge = self.hasher.clone().finalize();

        let lo = u128::from_be_bytes(challenge[0..16].try_into().unwrap());
        let hi = u128::from_be_bytes(challenge[16..32].try_into().unwrap());

        self.hasher.update([0x00]);
        self.hasher.update(challenge);
        self.hasher.update([0x01]);

        (lo, hi)
    }

    #[allow(clippy::not_unsafe_ptr_arg_deref)]
    pub fn get_challenge<'cfg, const N: usize>(
        &mut self,
        config: ConfigRef<'cfg, N>,
    ) -> RandomField<'cfg, N> {
        let (lo, hi) = self.get_challenge_limbs();
        let config = config.reference().expect("Field config cannot be none");
        let modulus = config.modulus();
        let challenge_num_bits = modulus.num_bits() - 1;
        if N == 1 {
            let lo_mask = (1u64 << challenge_num_bits) - 1;

            let truncated_lo = lo as u64 & lo_mask;

            let mut challenge = truncated_lo.map_to_field(ConfigRef::from(config));
            challenge.set_config(ConfigRef::from(config));
            return challenge;
        }
        if challenge_num_bits < 128 {
            let lo_mask = (1u128 << challenge_num_bits) - 1;

            let truncated_lo = lo & lo_mask;

            let mut challenge = truncated_lo.map_to_field(ConfigRef::from(config));
            challenge.set_config(ConfigRef::from(config));
            challenge
        } else if challenge_num_bits >= 256 {
            let two_to_128 =
                BigInt::<N>::from_bits_le(&(0..196).map(|i| i == 128).collect::<Vec<bool>>())
                    .map_to_field(ConfigRef::from(config));

            let mut challenge = lo.map_to_field(ConfigRef::from(config))
                + two_to_128 * hi.map_to_field(ConfigRef::from(config));
            challenge.set_config(ConfigRef::from(config));
            challenge
        } else {
            let hi_bits_to_keep = challenge_num_bits - 128;
            let hi_mask = (1u128 << hi_bits_to_keep) - 1;

            let truncated_hi = hi & hi_mask;

            let two_to_128 =
                BigInt::<N>::from_bits_le(&(0..196).map(|i| i == 128).collect::<Vec<bool>>())
                    .map_to_field(ConfigRef::from(config));

            let mut ret = lo.map_to_field(ConfigRef::from(config))
                + two_to_128 * truncated_hi.map_to_field(ConfigRef::from(config));
            ret.set_config(ConfigRef::from(config));
            ret
        }
    }
    pub fn get_challenges<'cfg, const N: usize>(
        &mut self,
        n: usize,
        config: ConfigRef<'cfg, N>,
    ) -> Vec<RandomField<'cfg, N>> {
        let mut challenges = Vec::with_capacity(n);
        challenges.extend((0..n).map(|_| self.get_challenge(config)));
        challenges
    }

    pub fn get_integer_challenge<const N: usize>(&mut self) -> Int<N> {
        let mut words = [0u64; N];

        for word in &mut words {
            let mut challenge = [0u8; 8];
            challenge.copy_from_slice(self.get_random_bytes(8).as_slice());
            self.hasher.update([0x12]);
            self.hasher.update(challenge);
            self.hasher.update([0x34]);
            *word = u64::from_le_bytes(challenge);
        }

        Int::<N>::from_words(words)
    }

    pub fn get_integer_challenges<const N: usize>(&mut self, n: usize) -> Vec<Int<N>> {
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
impl<const L: usize> ZipTranscript<L> for KeccakTranscript {
    fn get_encoding_element(&mut self) -> Int<L> {
        self.get_integer_challenge::<L>()
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
        biginteger::BigInt,
        field_config::{ConfigRef, FieldConfig},
        traits::FieldMap,
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
        let challenge = transcript.get_challenge(field_config);

        let expected_bigint = BigInt::<32>::from_str(
            "693058076479703886486101269644733982722902192016595549603371045888466087870",
        )
        .unwrap();
        let expected_field = expected_bigint.map_to_field(field_config);

        assert_eq!(challenge, expected_field);
    }
}
