mod arithmetic;
mod comparison;
mod conversion;
mod numeric;

use crypto_bigint::{
    Random, Uint, Word,
    modular::{ConstMontyForm, ConstMontyParams},
    rand_core::TryRngCore,
};

use crate::traits::{Field, FromBytes, ToBytes};

pub trait Monty<const LIMBS: usize>: ConstMontyParams<LIMBS> {}

pub(crate) type Form<MOD, const LIMBS: usize> = ConstMontyForm<MOD, LIMBS>;

impl<T: ConstMontyParams<LIMBS>, const LIMBS: usize> Monty<LIMBS> for T {}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct RandomField<MOD: Monty<LIMBS>, const LIMBS: usize>(pub(crate) Form<MOD, LIMBS>);

impl<MOD: Monty<LIMBS>, const LIMBS: usize> Field<LIMBS> for RandomField<MOD, LIMBS> {
    type Monty = MOD;
}

impl<MOD: Monty<LIMBS>, const LIMBS: usize> Random for RandomField<MOD, LIMBS> {
    fn try_random<R: TryRngCore + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
        Ok(Self(Form::try_random(rng)?))
    }
}

impl<const LIMBS: usize> ToBytes for Uint<LIMBS> {
    fn to_be_bytes(&self) -> Vec<u8> {
        // crypto_bigint::Uint exposes words in little-endian limb order (least significant first).
        // For a big-endian byte array of the whole integer, we must output the most significant
        // word first, each serialized in big-endian.
        self.as_words()
            .iter()
            .rev()
            .flat_map(|w| w.to_be_bytes())
            .collect()
    }

    fn to_le_bytes(&self) -> Vec<u8> {
        // Little-endian byte array: least significant word first, each in little-endian.
        self.as_words()
            .iter()
            .flat_map(|w| w.to_le_bytes())
            .collect()
    }
}

impl<MOD: Monty<LIMBS>, const LIMBS: usize> ToBytes for RandomField<MOD, LIMBS> {
    fn to_be_bytes(&self) -> Vec<u8> {
        ToBytes::to_be_bytes(self.0.as_montgomery())
    }

    fn to_le_bytes(&self) -> Vec<u8> {
        ToBytes::to_le_bytes(self.0.as_montgomery())
    }
}

impl<MOD: Monty<LIMBS>, const LIMBS: usize> FromBytes for RandomField<MOD, LIMBS> {
    fn from_bytes_le(bytes: &[u8]) -> Option<Self> {
        // Interpret input as Montgomery representation, mirroring ToBytes
        let u = Uint::<LIMBS>::from_bytes_le(bytes)?;
        Some(Self(Form::from_montgomery(u)))
    }

    fn from_bytes_be(bytes: &[u8]) -> Option<Self> {
        // Interpret input as Montgomery representation, mirroring ToBytes
        let u = Uint::<LIMBS>::from_bytes_be(bytes)?;
        Some(Self(Form::from_montgomery(u)))
    }
}

impl<const LIMBS: usize> FromBytes for Uint<LIMBS> {
    fn from_bytes_le(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != LIMBS * size_of::<Word>() {
            None
        } else {
            Some(Uint::from_le_slice(bytes))
        }
    }

    fn from_bytes_be(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != LIMBS * size_of::<Word>() {
            None
        } else {
            Some(Uint::from_be_slice(bytes))
        }
    }
}

#[cfg(target_pointer_width = "64")]
pub const WORD_FACTOR: usize = 1;
#[cfg(target_pointer_width = "32")]
pub const WORD_FACTOR: usize = 2;

pub type F64<MOD> = RandomField<MOD, { WORD_FACTOR }>;
pub type F128<MOD> = RandomField<MOD, { 2 * WORD_FACTOR }>;
pub type F192<MOD> = RandomField<MOD, { 3 * WORD_FACTOR }>;
pub type F256<MOD> = RandomField<MOD, { 4 * WORD_FACTOR }>;
pub type F320<MOD> = RandomField<MOD, { 5 * WORD_FACTOR }>;
pub type F384<MOD> = RandomField<MOD, { 6 * WORD_FACTOR }>;
pub type F448<MOD> = RandomField<MOD, { 7 * WORD_FACTOR }>;
pub type F512<MOD> = RandomField<MOD, { 8 * WORD_FACTOR }>;
pub type F576<MOD> = RandomField<MOD, { 9 * WORD_FACTOR }>;
pub type F640<MOD> = RandomField<MOD, { 10 * WORD_FACTOR }>;
pub type F704<MOD> = RandomField<MOD, { 11 * WORD_FACTOR }>;
pub type F768<MOD> = RandomField<MOD, { 12 * WORD_FACTOR }>;
pub type F832<MOD> = RandomField<MOD, { 13 * WORD_FACTOR }>;
pub type F896<MOD> = RandomField<MOD, { 14 * WORD_FACTOR }>;
pub type F960<MOD> = RandomField<MOD, { 15 * WORD_FACTOR }>;
pub type F1024<MOD> = RandomField<MOD, { 16 * WORD_FACTOR }>;
pub type F1280<MOD> = RandomField<MOD, { 20 * WORD_FACTOR }>;
pub type F1536<MOD> = RandomField<MOD, { 24 * WORD_FACTOR }>;
pub type F1792<MOD> = RandomField<MOD, { 28 * WORD_FACTOR }>;
pub type F2048<MOD> = RandomField<MOD, { 32 * WORD_FACTOR }>;
pub type F3072<MOD> = RandomField<MOD, { 48 * WORD_FACTOR }>;
pub type F3584<MOD> = RandomField<MOD, { 56 * WORD_FACTOR }>;
pub type F4096<MOD> = RandomField<MOD, { 64 * WORD_FACTOR }>;
pub type F4224<MOD> = RandomField<MOD, { 66 * WORD_FACTOR }>;
pub type F4352<MOD> = RandomField<MOD, { 68 * WORD_FACTOR }>;
pub type F6144<MOD> = RandomField<MOD, { 96 * WORD_FACTOR }>;
pub type F8192<MOD> = RandomField<MOD, { 128 * WORD_FACTOR }>;
pub type F16384<MOD> = RandomField<MOD, { 256 * WORD_FACTOR }>;
pub type F32768<MOD> = RandomField<MOD, { 512 * WORD_FACTOR }>;

#[cfg(test)]
mod tests {
    use crypto_bigint::{U128, const_monty_params};

    use super::*;

    const_monty_params!(ModP, U128, "7fffffffffffffffffffffffffffffff");

    type F = RandomField<ModP, { U128::LIMBS }>;

    #[test]
    fn basic_add_smoke() {
        let a: F = 123u64.into();
        let b: F = 456u64.into();
        assert_eq!(a + b, F::from(579u64));
    }
}
