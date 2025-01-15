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
        writer: W,
        _: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        let output_byte_size = buffer_byte_size(N as usize * 64);

        // Write out `self` to a temporary buffer.
        // The size of the buffer is $byte_size + 1 because `F::BIT_SIZE`
        // is at most 8 bits.
        let mut bytes = crate::const_helpers::SerBuffer::zeroed();
        bytes.copy_from_u64_slice(&self.into_bigint().0);

        bytes.write_up_to(writer, output_byte_size)?;

        if self.config.is_none() {
            bytes[output_byte_size - 1] |= 1;
        } else {
            let config_ptr: *const FieldConfig<N> = self.config.unwrap();
            bytes.buffers[N] = (config_ptr as u64).to_le_bytes()
        }
        Ok(())
    }

    fn serialized_size(&self, _: ark_serialize::Compress) -> usize {
        buffer_byte_size(N as usize * 64)
    }
}

impl<const N: usize> Valid for RandomField<'_, N> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}
impl<const N: usize> CanonicalDeserialize for RandomField<'_, N> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        // Calculate the number of bytes required to represent a field element
        // serialized with `flags`.
        let output_byte_size = buffer_byte_size(N as usize * 64 + 64);

        let mut masked_bytes = crate::const_helpers::SerBuffer::zeroed();
        masked_bytes.read_exact_up_to(reader, output_byte_size)?;
        if masked_bytes.last == 1 {
            masked_bytes.last = 0;
            return Ok(Self::new_unchecked(None, masked_bytes.to_bigint()));
        }

        let ptr: *const FieldConfig<N> = masked_bytes.buffers[N].as_ptr() as *const FieldConfig<N>;

        masked_bytes.buffers[N] = [0u8; 8];
        let value = masked_bytes.to_bigint();
        Ok(Self::new_unchecked(unsafe { ptr.as_ref() }, value))
    }
}
#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::{BigInt, BigInteger256, One, Zero};
    use ark_serialize::CanonicalSerialize;

    use crate::field_config::FieldConfig;

    use super::RandomField;

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
        use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
        use ark_std::io::Cursor;

        let field_config = FieldConfig::new(
            BigInteger256::from_str("695962179703626800597079116051991347").unwrap(),
        );

        let mut serialized = Vec::new();
        let bigint = RandomField::<'_, 4>::new_unchecked(Some(&field_config), BigInt::zero());

        let mut cursor = Cursor::new(&serialized);
        bigint.serialize_with_mode(&mut serialized, Compress::No);
        println!("{:?}", serialized)
    }
}
