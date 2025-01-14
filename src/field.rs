use std::ops::Add;

use crate::field_config::{self, FieldConfig};

use ark_ff::{BigInt, BigInteger};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Valid};

pub struct RandomField<'config, const N: usize> {
    pub config: &'config FieldConfig<N>,
    pub value: BigInt<N>,
}

impl<'config, const N: usize> RandomField<'config, N> {
    pub fn new_unchecked(config: &'config FieldConfig<N>, value: BigInt<N>) -> Self {
        RandomField { config, value }
    }
    /// Convert from `BigInteger` to `RandomField`
    ///
    /// If `BigInteger` is greater then field modulus return `None`
    pub fn from_bigint(config: &'config FieldConfig<N>, value: BigInt<N>) -> Option<Self> {
        if value.is_zero() {
            Some(Self::new_unchecked(config, value))
        } else if value >= config.modulus {
            None
        } else {
            let mut r = value;
            config.mul_assign(&mut r, &config.r2);
            Some(Self::new_unchecked(config, r))
        }
    }

    pub fn into_bigint(&self) -> BigInt<N> {
        let mut r = self.value.0;
        // Montgomery Reduction
        for i in 0..N {
            let k = r[i].wrapping_mul(self.config.inv);
            let mut carry = 0;

            field_config::mac_with_carry(r[i], k, self.config.modulus.0[0], &mut carry);
            for j in 1..N {
                r[(j + i) % N] = field_config::mac_with_carry(
                    r[(j + i) % N],
                    k,
                    self.config.modulus.0[j],
                    &mut carry,
                );
            }
            r[i % N] = carry;
        }

        BigInt::new(r)
    }
}

impl<'a, 'config, const N: usize> Add<&'a RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        &self + rhs
    }
}

impl<'a, 'config, const N: usize> Add<&'a RandomField<'config, N>> for &RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.

        let config_ptr_lhs: *const FieldConfig<N> = self.config;
        let config_ptr_rhs: *const FieldConfig<N> = rhs.config;

        if config_ptr_lhs != config_ptr_rhs {
            panic!("cannot add field elements of different fields");
        }

        todo!()
    }
}

impl<const N: usize> std::fmt::Debug for RandomField<'_, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl<const N: usize> CanonicalSerialize for RandomField<'_, N> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        todo!()
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        todo!()
    }
}

impl<const N: usize> Valid for RandomField<'_, N> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        todo!()
    }
}
impl<const N: usize> CanonicalDeserialize for RandomField<'_, N> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        todo!()
    }
}
#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::BigInteger256;

    use crate::field_config::FieldConfig;

    use super::RandomField;

    #[test]
    fn test_bigint_conversion() {
        let field_config = FieldConfig::new(
            BigInteger256::from_str("695962179703626800597079116051991347").unwrap(),
            BigInteger256::from_str("2").unwrap(),
        );

        let bigint = BigInteger256::from_str("695962179703").unwrap();

        let field_elem = RandomField::from_bigint(&field_config, bigint).unwrap();
        assert_eq!(bigint, field_elem.into_bigint());
        let bigint = BigInteger256::from_str("695962179703626800597079116051991346").unwrap();

        let field_elem = RandomField::from_bigint(&field_config, bigint).unwrap();
        assert_eq!(bigint, field_elem.into_bigint())
    }
}
