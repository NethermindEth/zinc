use sha3::{Digest, Keccak256};

use crate::{
    biginteger::BigInt,
    field::{conversion::FieldMap, RandomField},
    field_config::FieldConfig,
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

    pub fn absorb_random_field<const N: usize>(&mut self, v: &RandomField<N>) {
        match v {
            RandomField::Raw { value } => {
                self.absorb(&[0x1]);
                self.absorb(&value.to_bytes_be());
                self.absorb(&[0x3])
            }
            RandomField::Initialized { config, value } => {
                unsafe {
                    let config = config.as_ref().unwrap();
                    self.absorb(&[0x3]);
                    self.absorb(&config.modulus.to_bytes_be());
                    self.absorb(&[0x5])
                }
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
    pub fn get_challenge<const N: usize>(
        &mut self,
        config: *const FieldConfig<N>,
    ) -> RandomField<N> {
        let (lo, hi) = self.get_challenge_limbs();
        let modulus = unsafe { (*config).modulus };
        let challenge_num_bits = modulus.num_bits() - 1;
        if N == 1 {
            let lo_mask = (1u64 << challenge_num_bits) - 1;

            let truncated_lo = lo as u64 & lo_mask;

            let mut challenge = truncated_lo.map_to_field(config);
            challenge.set_config(config);
            return challenge;
        }
        if challenge_num_bits < 128 {
            let lo_mask = (1u128 << challenge_num_bits) - 1;

            let truncated_lo = lo & lo_mask;

            let mut challenge = truncated_lo.map_to_field(config);
            challenge.set_config(config);
            challenge
        } else if challenge_num_bits >= 256 {
            let two_to_128 =
                BigInt::<N>::from_bits_le(&(0..196).map(|i| i == 128).collect::<Vec<bool>>())
                    .map_to_field(config);

            let mut challenge = lo.map_to_field(config) + two_to_128 * hi.map_to_field(config);
            challenge.set_config(config);
            challenge
        } else {
            let hi_bits_to_keep = challenge_num_bits - 128;
            let hi_mask = (1u128 << hi_bits_to_keep) - 1;

            let truncated_hi = hi & hi_mask;

            let two_to_128 =
                BigInt::<N>::from_bits_le(&(0..196).map(|i| i == 128).collect::<Vec<bool>>())
                    .map_to_field(config);

            let mut ret = lo.map_to_field(config) + two_to_128 * truncated_hi.map_to_field(config);
            ret.set_config(config);
            ret
        }
    }
    pub fn get_challenges<const N: usize>(
        &mut self,
        n: usize,
        config: *const FieldConfig<N>,
    ) -> Vec<RandomField<N>> {
        let mut challenges = Vec::with_capacity(n);
        challenges.extend((0..n).map(|_| self.get_challenge(config)));
        challenges
    }

    pub fn get_integer_challenge(&mut self) -> i64 {
        let challenge = self.hasher.clone().finalize();

        let int = i64::from_be_bytes(challenge[0..8].try_into().unwrap());

        self.hasher.update([0x12]);
        self.hasher.update(challenge);
        self.hasher.update([0x34]);

        int
    }

    pub fn get_integer_challenges<const N: usize>(&mut self, n: usize) -> Vec<Int<N>> {
        (0..n).map(|_| self.get_integer_challenge()).collect()
    }
    fn get_usize_in_range(&mut self, range: &std::ops::Range<usize>) -> usize {
        let challenge = self.hasher.clone().finalize();

        self.hasher.update([0x88]);
        self.hasher.update(challenge);
        self.hasher.update([0x11]);

        let num = usize::from_le_bytes(challenge[..8].try_into().unwrap());
        range.start + (num % (range.end - range.start))
    }
}
impl ZipTranscript for KeccakTranscript {
    fn get_encoding_element(&mut self) -> i128 {
        let challenge = self.hasher.clone().finalize();

        let int = i128::from_be_bytes(challenge[0..16].try_into().unwrap());

        self.hasher.update([0x89]);
        self.hasher.update(challenge);
        self.hasher.update([0x11]);

        int
    }

    fn sample_unique_columns(
        &mut self,
        range: std::ops::Range<usize>,
        columns: &mut std::collections::BTreeSet<usize>,
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
    use std::str::FromStr;

    use crate::{biginteger::BigInt, field::conversion::FieldMap, field_config::FieldConfig};

    use super::KeccakTranscript;

    #[test]
    fn test_keccak_transcript() {
        let mut transcript = KeccakTranscript::new();
        let field_config = FieldConfig::new(
            BigInt::<32>::from_str(
                "3618502788666131213697322783095070105623107215331596699973092056135872020481",
            )
            .unwrap(),
        );

        transcript.absorb(b"This is a test string!");
        let challenge = transcript.get_challenge(&field_config);

        let expected_bigint = BigInt::<32>::from_str(
            "693058076479703886486101269644733982722902192016595549603371045888466087870",
        )
        .unwrap();
        let expected_field = expected_bigint.map_to_field(&field_config);

        assert_eq!(challenge, expected_field);
    }
}
