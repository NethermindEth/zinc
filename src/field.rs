use std::ops::{Add, Mul};

use ark_ff::{BigInt, BigInteger, One, Zero};

use crate::field_config::{self, FieldConfig};

use ark_serialize::{buffer_byte_size, CanonicalDeserialize, CanonicalSerialize, Valid};

#[derive(Copy, Clone)]
pub struct RandomField<'config, const N: usize> {
    pub config: Option<&'config FieldConfig<N>>,
    pub value: BigInt<N>,
}

impl<'config, const N: usize> RandomField<'config, N> {
    fn new_unchecked(config: Option<&'config FieldConfig<N>>, value: BigInt<N>) -> Self {
        RandomField { config, value }
    }
    /// Convert from `BigInteger` to `RandomField`
    ///
    /// If `BigInteger` is greater then field modulus return `None`
    pub fn from_bigint(config: &'config FieldConfig<N>, value: BigInt<N>) -> Option<Self> {
        if value.is_zero() {
            Some(Self::new_unchecked(Some(config), value))
        } else if value >= config.modulus {
            None
        } else {
            let mut r = value;
            config.mul_assign(&mut r, &config.r2);
            Some(Self::new_unchecked(Some(config), r))
        }
    }

    pub fn into_bigint(&self) -> BigInt<N> {
        if self.is_zero() {
            return BigInt::zero();
        }

        if self.value == BigInt::one() && self.config.is_none() {
            return BigInt::one();
        }

        let config = self
            .config
            .expect("This field element has no associated field");
        let mut r = self.value.0;
        // Montgomery Reduction
        for i in 0..N {
            let k = r[i].wrapping_mul(config.inv);
            let mut carry = 0;

            field_config::mac_with_carry(r[i], k, config.modulus.0[0], &mut carry);
            for j in 1..N {
                r[(j + i) % N] = field_config::mac_with_carry(
                    r[(j + i) % N],
                    k,
                    config.modulus.0[j],
                    &mut carry,
                );
            }
            r[i % N] = carry;
        }

        BigInt::new(r)
    }
}

impl<'config, const N: usize> Add<RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: RandomField<'config, N>) -> RandomField<'config, N> {
        &self + &rhs
    }
}

impl<'a, 'config, const N: usize> Add<&'a RandomField<'config, N>> for &RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        if rhs.is_zero() {
            return *self;
        }
        if self.is_zero() {
            return *rhs;
        }
        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.
        let lconfig = self
            .config
            .expect("This field element has no associated field");
        let rconfig = rhs
            .config
            .expect("This field element has no associated field");
        let config_ptr_lhs: *const FieldConfig<N> = lconfig;
        let config_ptr_rhs: *const FieldConfig<N> = rconfig;

        if config_ptr_lhs != config_ptr_rhs {
            panic!("cannot add field elements of different fields");
        }

        todo!()
    }
}

impl<'config, const N: usize> Mul<RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn mul(self, _: RandomField<'config, N>) -> RandomField<'config, N> {
        todo!()
    }
}

impl<const N: usize> std::fmt::Debug for RandomField<'_, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.config {
            Some(config) => write!(f, "{} in the field Z_{}", self.value, config.modulus),
            None => write!(f, "{}", self.value),
        }
    }
}

impl<const N: usize> Zero for RandomField<'_, N> {
    fn zero() -> Self {
        Self::new_unchecked(None, BigInt::zero())
    }

    fn is_zero(&self) -> bool {
        self.value == BigInt::zero()
    }

    fn set_zero(&mut self) {
        self.value = BigInt::zero()
    }
}

impl<const N: usize> One for RandomField<'_, N> {
    fn one() -> Self {
        Self::new_unchecked(None, BigInt::one())
    }

    fn set_one(&mut self) {
        self.value = BigInt::one()
    }

    fn is_one(&self) -> bool {
        match self.config {
            Some(conf) => self.value == conf.r,
            None => self.value == BigInt::one(),
        }
    }
}

impl<const N: usize> PartialEq for RandomField<'_, N> {
    fn eq(&self, other: &Self) -> bool {
        if self.config.is_none() && other.config.is_none() {
            return self.is_one() && other.is_one() || self.is_zero() && other.is_zero();
        }

        let config_ptr_lhs: *const FieldConfig<N> = self.config.unwrap();
        let config_ptr_rhs: *const FieldConfig<N> = other.config.unwrap();

        self.value == other.value && config_ptr_lhs == config_ptr_rhs
    }
}

impl<const N: usize> Eq for RandomField<'_, N> {}

unsafe impl<const N: usize> Send for RandomField<'_, N> {}

unsafe impl<const N: usize> Sync for RandomField<'_, N> {}

impl<const N: usize> CanonicalSerialize for RandomField<'_, N> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        _: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        let mut bytes = vec![0u8; N * 8 + 9];

        if self.config.is_none() {
            bytes[N * 8 + 8] = 1;

            if self.is_one() {
                bytes[0] = 1;
            }
            writer.write_all(&bytes)?;
            return Ok(());
        } else {
            let config_ptr: *const FieldConfig<N> = self.config.unwrap();

            // Convert the pointer to a u64 and write the bytes
            let ptr_bytes = (config_ptr as u64).to_be_bytes();
            bytes[N * 8..N * 8 + 8].copy_from_slice(&ptr_bytes);
        }

        // Fill in the bytes for the limbs
        self.value.0.iter().enumerate().for_each(|(i, limb)| {
            let limb_bytes = limb.to_be_bytes();
            bytes[i * 8..(i + 1) * 8].copy_from_slice(&limb_bytes);
        });

        // Write the entire byte slice to the writer
        writer.write_all(&bytes)?;

        Ok(())
    }

    fn serialized_size(&self, _: ark_serialize::Compress) -> usize {
        buffer_byte_size(N * 64)
    }
}

impl<const N: usize> Valid for RandomField<'_, N> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

impl<const N: usize> CanonicalDeserialize for RandomField<'_, N> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        mut reader: R,
        _compress: ark_serialize::Compress,
        _validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        let mut bytes = Vec::with_capacity(N * 8 + 9);

        reader.read_to_end(&mut bytes)?;

        if bytes[N * 8 + 8] == 1 {
            if bytes[0] == 1 {
                return Ok(Self::one());
            } else {
                return Ok(Self::zero());
            }
        }

        let mut value = BigInt::<N>::from(0u64);
        value
            .0
            .iter_mut()
            .zip(bytes[0..N * 8].chunks(8))
            .for_each(|(other, this)| {
                *other = u64::from_be_bytes(this.try_into().expect("Slice has incorrect length"));
            });

        let ptr_bytes = &bytes[N * 8..N * 8 + 8];

        let address = u64::from_be_bytes(ptr_bytes.try_into().expect("Invalid slice length"))
            as *const FieldConfig<N>;

        Ok(Self::new_unchecked(unsafe { address.as_ref() }, value))
    }
}
#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::{BigInt, BigInteger256, BigInteger64, One, Zero};
    use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};

    use super::RandomField;
    use crate::field_config::FieldConfig;

    use ark_std::io::Cursor;

    #[test]
    fn test_bigint_conversion() {
        let field_config = FieldConfig::new(
            BigInteger256::from_str("695962179703626800597079116051991347").unwrap(),
        );

        let bigint = BigInteger256::from_str("695962179703").unwrap();

        let field_elem = RandomField::from_bigint(&field_config, bigint).unwrap();
        assert_eq!(bigint, field_elem.into_bigint());
        let bigint = BigInteger256::from_str("695962179703626800597079116051991346").unwrap();

        let field_elem = RandomField::from_bigint(&field_config, bigint).unwrap();
        assert_eq!(bigint, field_elem.into_bigint())
    }

    #[test]
    fn test_serialization() {
        let field_config = FieldConfig::new(BigInteger64::from_str("17").unwrap());

        let mut serialized = Vec::new();
        let field_elem =
            RandomField::<'_, 1>::new_unchecked(Some(&field_config), BigInt::from(13u32));

        let _ = field_elem.serialize_with_mode(&mut serialized, Compress::No);

        let cursor = Cursor::new(&serialized);
        let deser_field_elem =
            RandomField::<'_, 1>::deserialize_with_mode(cursor, Compress::No, Validate::No)
                .unwrap();

        assert_eq!(field_elem, deser_field_elem);
    }

    #[test]
    fn test_zero_serialization() {
        let mut serialized = Vec::new();
        let field_elem = RandomField::<'_, 1>::zero();

        let _ = field_elem.serialize_with_mode(&mut serialized, Compress::No);

        let cursor = Cursor::new(&serialized);
        let deser_field_elem =
            RandomField::<'_, 1>::deserialize_with_mode(cursor, Compress::No, Validate::No)
                .unwrap();

        assert!(deser_field_elem.is_zero());
    }

    #[test]
    fn test_one_serialization() {
        let mut serialized = Vec::new();
        let field_elem = RandomField::<'_, 1>::one();

        let _ = field_elem.serialize_with_mode(&mut serialized, Compress::No);

        let cursor = Cursor::new(&serialized);
        let deser_field_elem =
            RandomField::<'_, 1>::deserialize_with_mode(cursor, Compress::No, Validate::No)
                .unwrap();

        assert!(deser_field_elem.is_one());
    }
}
