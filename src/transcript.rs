use sha3::{Digest, Keccak256};

use crate::{biginteger::BigInt, field::RandomField, field_config::FieldConfig};

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

    pub fn get_challenge_limbs(&mut self) -> (u128, u128) {
        let challenge = self.hasher.clone().finalize();

        let lo = u128::from_be_bytes(challenge[0..16].try_into().unwrap());
        let hi = u128::from_be_bytes(challenge[16..32].try_into().unwrap());

        self.hasher.update([0x00]);
        self.hasher.update(challenge);
        self.hasher.update([0x01]);

        (lo, hi)
    }

    pub fn get_challenge<const N: usize>(&mut self, config: &FieldConfig<N>) -> RandomField<N> {
        let (lo, hi) = self.get_challenge_limbs();
        let (lo, hi) = (RandomField::from(lo), RandomField::from(hi));

        let two_to_128 = RandomField::from_bigint(
            config,
            BigInt::from_bits_le(&(0..196).map(|i| i == 128).collect::<Vec<bool>>()),
        )
        .unwrap();

        lo + two_to_128 * hi
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::Zero;

    use crate::{biginteger::BigInt, field::RandomField, field_config::FieldConfig};

    use super::KeccakTranscript;

    #[test]
    fn test_keccak_transcript() {
        let mut transcript = KeccakTranscript::new();
        let field_cofig = FieldConfig::new(
            BigInt::<32>::from_str(
                "3618502788666131213697322783095070105623107215331596699973092056135872020481",
            )
            .unwrap(),
        );

        transcript.absorb(b"This is a test string!");
        let challenge = transcript.get_challenge(&field_cofig);

        // TODO: fill in the appropriate value once From<u128> is implemented
        // for RandomField
        assert_eq!(challenge, RandomField::zero());
    }
}
