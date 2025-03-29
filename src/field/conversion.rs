use i256::{I256, I512};

use crate::biginteger::BigInt;
use crate::field::RandomField;
use crate::field::RandomField::Raw;
use crate::field_config::FieldConfig;
use crate::traits::FromBytes;

impl<const N: usize> From<u128> for RandomField<N> {
    fn from(other: u128) -> Self {
        let mut value = BigInt::default();
        if N == 1 {
            panic!("Integer is 128 bits but field is 64 bits")
        } else {
            value.0[0] = ((other << 64) >> 64) as u64;
            value.0[1] = (other >> 64) as u64;
        }
        Raw { value }
    }
}

macro_rules! impl_from_uint {
    ($type:ty) => {
        impl<const N: usize> From<$type> for RandomField<N> {
            fn from(value: $type) -> Self {
                let value = BigInt::from(value);
                Raw { value }
            }
        }
    };
}

impl_from_uint!(u64);
impl_from_uint!(u32);
impl_from_uint!(u16);
impl_from_uint!(u8);

impl<const N: usize> From<bool> for RandomField<N> {
    fn from(value: bool) -> Self {
        let value = BigInt::from(value as u8);
        Raw { value }
    }
}

impl<const N: usize> FromBytes for RandomField<N> {
    fn from_bytes_le(bytes: &[u8]) -> Option<Self> {
        Some(Raw {
            value: BigInt::<N>::from_bytes_le(bytes)?,
        })
    }

    fn from_bytes_be(bytes: &[u8]) -> Option<Self> {
        Some(Raw {
            value: BigInt::<N>::from_bytes_be(bytes)?,
        })
    }
}

impl<const N: usize> RandomField<N> {
    pub fn from_bytes_le_with_config(config: *const FieldConfig<N>, bytes: &[u8]) -> Option<Self> {
        let value = BigInt::<N>::from_bytes_le(bytes);

        Self::from_bigint(config, value?)
    }

    pub fn from_bytes_be_with_config(config: *const FieldConfig<N>, bytes: &[u8]) -> Option<Self> {
        let value = BigInt::<N>::from_bytes_be(bytes);

        Self::from_bigint(config, value?)
    }
}

pub trait FieldMap {
    type Output<const N: usize>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N>;
}

macro_rules! impl_field_map_for_int {
    ($type:ty, $bits:expr) => {
        impl FieldMap for $type {
            type Output<const N: usize> = RandomField<N>;
            fn map_to_field<const N: usize>(
                &self,
                config: *const FieldConfig<N>,
            ) -> Self::Output<N> {
                if config.is_null() {
                    panic!("Cannot convert signed integer to prime field element without a modulus")
                }
                unsafe {
                    let modulus: [u64; N] = (*config).modulus.0;
                    let abs_val = self.unsigned_abs();

                    // Calculate how many u64 limbs we need based on bits
                    const LIMBS: usize = ($bits + 63) / 64;
                    let mut val = [0u64; LIMBS];

                    // Fill val array based on size
                    if LIMBS == 1 {
                        val[0] = abs_val as u64;
                    } else {
                        for i in 0..LIMBS {
                            val[i] = (abs_val >> (i * 64)) as u64;
                        }
                    }

                    let mut r: BigInt<N> = match N {
                        n if n < LIMBS => {
                            let mut wider_modulus = [0u64; LIMBS];
                            wider_modulus[..N].copy_from_slice(&modulus);
                            let mut value = crypto_bigint::Uint::<LIMBS>::from_words(val);
                            let modu = crypto_bigint::Uint::<LIMBS>::from_words(wider_modulus);

                            value %= crypto_bigint::NonZero::from_uint(modu);
                            let mut result = [0u64; N];
                            result.copy_from_slice(&value.to_words()[..N]);

                            BigInt(result)
                        }
                        n if n == LIMBS => {
                            let mut value_N = [0u64; N];
                            value_N.copy_from_slice(&val);

                            let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                            let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                            value %= crypto_bigint::NonZero::from_uint(modu);
                            BigInt(value.to_words())
                        }
                        _ => {
                            let mut wider_value = [0u64; N];
                            wider_value[..LIMBS].copy_from_slice(&val);
                            let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                            let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                            wider %= crypto_bigint::NonZero::from_uint(modu);
                            BigInt(wider.to_words())
                        }
                    };

                    (*config).mul_assign(&mut r, &(*config).r2);

                    let mut elem = RandomField::<N>::new_unchecked(config, r);
                    if *self < 0 {
                        elem = -elem;
                    }

                    elem
                }
            }
        }

        impl FieldMap for &$type {
            type Output<const N: usize> = RandomField<N>;
            fn map_to_field<const N: usize>(
                &self,
                config: *const FieldConfig<N>,
            ) -> Self::Output<N> {
                (*self).map_to_field(config)
            }
        }
    };
}

// Usage with bit sizes
impl_field_map_for_int!(i8, 8);
impl_field_map_for_int!(i16, 16);
impl_field_map_for_int!(i32, 32);
impl_field_map_for_int!(i64, 64);
impl_field_map_for_int!(i128, 128);

// Separate implementation for I256
impl FieldMap for I256 {
    type Output<const N: usize> = RandomField<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        if config.is_null() {
            panic!("Cannot convert signed integer to prime field element without a modulus")
        }
        unsafe {
            let modulus: [u64; N] = (*config).modulus.0;
            let val: [u64; 4] = self.abs().to_le_u64();

            let mut r: BigInt<N> = match N {
                n if n < 4 => {
                    let mut wider_modulus: [u64; 4] = [0; 4];
                    wider_modulus[..N].copy_from_slice(&modulus);
                    let mut value = crypto_bigint::Uint::<4>::from_words(val);
                    let modu = crypto_bigint::Uint::<4>::from_words(wider_modulus);

                    value %= crypto_bigint::NonZero::from_uint(modu);
                    let mut result = [0u64; N];
                    result.copy_from_slice(&value.to_words()[..N]);

                    BigInt(result)
                }
                4 => {
                    let mut value_N: [u64; N] = [0; N];
                    value_N.copy_from_slice(&val);

                    let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    value %= crypto_bigint::NonZero::from_uint(modu);
                    BigInt(value.to_words())
                }
                _ => {
                    let mut wider_value: [u64; N] = [0; N];
                    wider_value[..4].copy_from_slice(&val);
                    let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    wider %= crypto_bigint::NonZero::from_uint(modu);
                    BigInt(wider.to_words())
                }
            };

            (*config).mul_assign(&mut r, &(*config).r2);

            let mut elem = RandomField::<N>::new_unchecked(config, r);
            if self.is_negative() {
                elem = -elem;
            }

            elem
        }
    }
}

// Reference implementation for I256
impl FieldMap for &I256 {
    type Output<const N: usize> = RandomField<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        (*self).map_to_field(config)
    }
}

// Implementation for I512
impl FieldMap for I512 {
    type Output<const N: usize> = RandomField<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        if config.is_null() {
            panic!("Cannot convert signed integer to prime field element without a modulus")
        }
        unsafe {
            let modulus: [u64; N] = (*config).modulus.0;
            let val: [u64; 8] = self.abs().to_le_u64();

            let mut r: BigInt<N> = match N {
                n if n < 8 => {
                    let mut wider_modulus: [u64; 8] = [0; 8];
                    wider_modulus[..N].copy_from_slice(&modulus);
                    let mut value = crypto_bigint::Uint::<8>::from_words(val);
                    let modu = crypto_bigint::Uint::<8>::from_words(wider_modulus);

                    value %= crypto_bigint::NonZero::from_uint(modu);
                    let mut result = [0u64; N];
                    result.copy_from_slice(&value.to_words()[..N]);

                    BigInt(result)
                }
                8 => {
                    let mut value_N: [u64; N] = [0; N];
                    value_N.copy_from_slice(&val);

                    let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    value %= crypto_bigint::NonZero::from_uint(modu);
                    BigInt(value.to_words())
                }
                _ => {
                    let mut wider_value: [u64; N] = [0; N];
                    wider_value[..8].copy_from_slice(&val);
                    let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    wider %= crypto_bigint::NonZero::from_uint(modu);
                    BigInt(wider.to_words())
                }
            };

            (*config).mul_assign(&mut r, &(*config).r2);

            let mut elem = RandomField::<N>::new_unchecked(config, r);
            if self.is_negative() {
                elem = -elem;
            }

            elem
        }
    }
}

impl FieldMap for &I512 {
    type Output<const N: usize> = RandomField<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        (*self).map_to_field(config)
    }
}

macro_rules! impl_field_map_for_uint {
    ($type:ty, $bits:expr) => {
        impl FieldMap for $type {
            type Output<const N: usize> = RandomField<N>;
            fn map_to_field<const N: usize>(
                &self,
                config: *const FieldConfig<N>,
            ) -> Self::Output<N> {
                if config.is_null() {
                    panic!(
                        "Cannot convert unsigned integer to prime field element without a modulus"
                    )
                }
                unsafe {
                    let modulus: [u64; N] = (*config).modulus.0;

                    // Calculate how many u64 limbs we need based on bits
                    const LIMBS: usize = ($bits + 63) / 64;
                    let mut val = [0u64; LIMBS];

                    // Fill val array based on size
                    if LIMBS == 1 {
                        val[0] = *self as u64;
                    } else {
                        for i in 0..LIMBS {
                            val[i] = (*self >> (i * 64)) as u64;
                        }
                    }

                    let mut r: BigInt<N> = match N {
                        n if n < LIMBS => {
                            let mut wider_modulus = [0u64; LIMBS];
                            wider_modulus[..N].copy_from_slice(&modulus);
                            let mut value = crypto_bigint::Uint::<LIMBS>::from_words(val);
                            let modu = crypto_bigint::Uint::<LIMBS>::from_words(wider_modulus);

                            value %= crypto_bigint::NonZero::from_uint(modu);
                            let mut result = [0u64; N];
                            result.copy_from_slice(&value.to_words()[..N]);

                            BigInt(result)
                        }
                        n if n == LIMBS => {
                            let mut value_N = [0u64; N];
                            value_N.copy_from_slice(&val);

                            let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                            let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                            value %= crypto_bigint::NonZero::from_uint(modu);
                            BigInt(value.to_words())
                        }
                        _ => {
                            let mut wider_value = [0u64; N];
                            wider_value[..LIMBS].copy_from_slice(&val);
                            let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                            let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                            wider %= crypto_bigint::NonZero::from_uint(modu);
                            BigInt(wider.to_words())
                        }
                    };

                    (*config).mul_assign(&mut r, &(*config).r2);
                    RandomField::<N>::new_unchecked(config, r)
                }
            }
        }

        impl FieldMap for &$type {
            type Output<const N: usize> = RandomField<N>;
            fn map_to_field<const N: usize>(
                &self,
                config: *const FieldConfig<N>,
            ) -> Self::Output<N> {
                (*self).map_to_field(config)
            }
        }
    };
}

// Usage with bit sizes
impl_field_map_for_uint!(u32, 32);
impl_field_map_for_uint!(u64, 64);

// Implementation for bool
impl FieldMap for bool {
    type Output<const N: usize> = RandomField<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        if config.is_null() {
            panic!("Cannot convert boolean to prime field element without a modulus")
        }
        unsafe {
            let mut r = BigInt::from(*self as u64);
            (*config).mul_assign(&mut r, &(*config).r2);
            RandomField::<N>::new_unchecked(config, r)
        }
    }
}

impl FieldMap for &bool {
    type Output<const N: usize> = RandomField<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        (*self).map_to_field(config)
    }
}

#[cfg(test)]
mod tests {
    use i256::{I256, I512};

    use crate::field::conversion::FieldMap;
    use crate::field_config::FieldConfig;
    use crate::traits::FromBytes;
    use crate::{biginteger::BigInt, create_field_config, field::RandomField};
    use std::str::FromStr;

    fn test_from<T: Clone, const N: usize>(value: T, value_str: &str)
    where
        RandomField<N>: From<T>,
    {
        let raw_element = RandomField::<N>::from(value);
        assert_eq!(
            raw_element,
            RandomField::Raw {
                value: BigInt::from_str(value_str).unwrap()
            }
        )
    }

    #[test]
    fn converts_u128_to_random_field() {
        test_from::<u128, 2>(
            243043087159742188419721163456177516,
            "243043087159742188419721163456177516",
        );
    }

    #[test]
    #[should_panic(expected = "Integer is 128 bits but field is 64 bits")]
    fn panics_when_u128_does_not_fit_in_n1() {
        test_from::<u128, 1>(243043087159742188419721163456177516, "");
    }

    #[test]
    fn converts_u64_to_random_field() {
        test_from::<u64, 1>(23, "23");
    }

    #[test]
    fn converts_u32_to_random_field() {
        test_from::<u32, 1>(23, "23");
    }

    #[test]
    fn converts_u16_to_random_field() {
        test_from::<u16, 1>(23, "23");
    }

    #[test]
    fn converts_u8_to_random_field() {
        test_from::<u8, 1>(23, "23");
    }

    #[test]
    fn converts_false_to_zero() {
        test_from::<bool, 1>(false, "0");
    }

    #[test]
    fn converts_true_to_one() {
        test_from::<bool, 1>(true, "1");
    }

    #[test]
    fn converts_from_bytes_le_with_config_valid() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const _;

        let bytes = [0x05, 0, 0, 0, 0, 0, 0, 0];
        let expected = BigInt::from_str("5").unwrap();

        let result = RandomField::<1>::from_bytes_le_with_config(config_ptr, &bytes).unwrap();
        assert_eq!(result.into_bigint(), expected);
    }

    #[test]
    fn converts_from_bytes_be_with_config_valid() {
        let config = FieldConfig::new(
            BigInt::<32>::from_str(
                "3618502788666131213697322783095070105623107215331596699973092056135872020481",
            )
            .unwrap(),
        );
        let config_ptr = &config as *const FieldConfig<32>;

        let bytes = [
            0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00,
        ]; // Value: 5 (big-endian)
        let expected = BigInt::<32>::from_bytes_be(&bytes).unwrap();

        let result = RandomField::<32>::from_bytes_be_with_config(config_ptr, &bytes);
        assert_eq!(result.map(|x| x.into_bigint()), Some(expected));
    }

    #[test]
    fn converts_from_bytes_le_with_config_zero() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const FieldConfig<1>;

        let bytes = [0x00; 8]; // All zeros
        let expected = RandomField::Initialized {
            config: config_ptr,
            value: BigInt::<1>::zero(),
        };

        let result = RandomField::<1>::from_bytes_le_with_config(config_ptr, &bytes);
        assert_eq!(result, Some(expected));
    }

    #[test]
    fn converts_from_bytes_be_with_config_zero() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const FieldConfig<1>;

        let bytes = [0x00; 8]; // All zeros
        let expected = RandomField::Initialized {
            config: config_ptr,
            value: BigInt::<1>::zero(),
        };

        let result = RandomField::<1>::from_bytes_be_with_config(config_ptr, &bytes);
        assert_eq!(result, Some(expected));
    }

    #[test]
    fn converts_from_bytes_le_with_config_out_of_range() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const FieldConfig<1>;

        let bytes = [0x65, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // Value: 101 (modulus is 23)
        let result = RandomField::<1>::from_bytes_le_with_config(config_ptr, &bytes);
        assert!(result.is_none());
    }

    #[test]
    fn converts_from_bytes_be_with_config_out_of_range() {
        let config = FieldConfig::new(BigInt::<32>::from_str("37129241769965749").unwrap());
        let config_ptr = &config as *const FieldConfig<32>;

        let bytes = [0x65, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // Value: 101
        let result = RandomField::<32>::from_bytes_be_with_config(config_ptr, &bytes);
        assert!(result.is_none());
    }

    #[test]
    fn converts_from_bytes_le_with_config_exact_modulus() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const FieldConfig<1>;

        let bytes = [0x17, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // Value: 23 (modulus)
        let result = RandomField::<1>::from_bytes_le_with_config(config_ptr, &bytes);
        assert!(result.is_none()); // Must be strictly less than modulus
    }

    #[test]
    fn converts_from_bytes_be_with_config_exact_modulus() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const FieldConfig<1>;

        let bytes = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x17]; // Value: 23 (big-endian)
        let result = RandomField::<1>::from_bytes_be_with_config(config_ptr, &bytes);
        assert!(result.is_none());
    }

    #[test]
    fn converts_from_bytes_le_with_config_leading_zeros() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const FieldConfig<1>;

        let bytes = [0b0000_0001]; // Value: 1 with leading zeros
        let expected =
            RandomField::from_bigint(config_ptr, BigInt::<1>::from_bytes_le(&bytes).unwrap());

        let result = RandomField::<1>::from_bytes_le_with_config(config_ptr, &bytes);
        assert_eq!(result, expected);
    }

    #[test]
    fn converts_from_bytes_be_with_config_leading_zeros() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const FieldConfig<1>;

        let bytes = [0x01]; //1 with leading zeros (big-endian);

        let result = RandomField::<1>::from_bytes_be_with_config(config_ptr, &bytes).unwrap();
        assert_eq!(result.into_bigint(), BigInt::one());
    }

    #[test]
    fn test_signed_integers_field_map() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const _;

        // Test i8
        let i8_val: i8 = 5;
        let i8_result = i8_val.map_to_field(config_ptr);
        assert_eq!(i8_result.into_bigint(), BigInt::from(5u64));

        let i8_neg: i8 = -5;
        let i8_neg_result = i8_neg.map_to_field(config_ptr);
        assert_eq!(i8_neg_result.into_bigint(), BigInt::from(18u64)); // 23 - 5 = 18

        // Test i16
        let i16_val: i16 = 10;
        let i16_result = i16_val.map_to_field(config_ptr);
        assert_eq!(i16_result.into_bigint(), BigInt::from(10u64));

        // Test i32
        let i32_val: i32 = 15;
        let i32_result = i32_val.map_to_field(config_ptr);
        assert_eq!(i32_result.into_bigint(), BigInt::from(15u64));

        // Test i64
        let i64_val: i64 = 20;
        let i64_result = i64_val.map_to_field(config_ptr);
        assert_eq!(i64_result.into_bigint(), BigInt::from(20u64));

        // Test i128
        let i128_val: i128 = 22;
        let i128_result = i128_val.map_to_field(config_ptr);
        assert_eq!(i128_result.into_bigint(), BigInt::from(22u64));

        // Test reference implementations
        let i8_ref_result = (&i8_val).map_to_field(config_ptr);
        assert_eq!(i8_ref_result.into_bigint(), BigInt::from(5u64));

        let i32_ref_result = (&i32_val).map_to_field(config_ptr);
        assert_eq!(i32_ref_result.into_bigint(), BigInt::from(15u64));
    }

    #[test]
    fn test_unsigned_integers_field_map() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const _;

        // Test u32
        let u32_val: u32 = 5;
        let u32_result = u32_val.map_to_field(config_ptr);
        assert_eq!(u32_result.into_bigint(), BigInt::from(5u64));

        // Test u64
        let u64_val: u64 = 10;
        let u64_result = u64_val.map_to_field(config_ptr);
        assert_eq!(u64_result.into_bigint(), BigInt::from(10u64));

        // Test reference implementations
        let u32_ref_result = (&u32_val).map_to_field(config_ptr);
        assert_eq!(u32_ref_result.into_bigint(), BigInt::from(5u64));

        let u64_ref_result = (&u64_val).map_to_field(config_ptr);
        assert_eq!(u64_ref_result.into_bigint(), BigInt::from(10u64));
    }

    #[test]
    fn test_i256_field_map() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const _;

        // Test positive I256
        let i256_val = I256::from_str("5").unwrap();
        let i256_result = i256_val.map_to_field(config_ptr);
        assert_eq!(i256_result.into_bigint(), BigInt::from(5u64));

        // Test negative I256
        let i256_neg = I256::from_str("-5").unwrap();
        let i256_neg_result = i256_neg.map_to_field(config_ptr);
        assert_eq!(i256_neg_result.into_bigint(), BigInt::from(18u64)); // 23 - 5 = 18

        // Test reference implementation
        let i256_ref_result = (&i256_val).map_to_field(config_ptr);
        assert_eq!(i256_ref_result.into_bigint(), BigInt::from(5u64));
    }

    #[test]
    fn test_i512_field_map() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const _;

        // Test positive I512
        let i512_val = I512::from_str("5").unwrap();
        let i512_result = i512_val.map_to_field(config_ptr);
        assert_eq!(i512_result.into_bigint(), BigInt::from(5u64));

        // Test negative I512
        let i512_neg = I512::from_str("-5").unwrap();
        let i512_neg_result = i512_neg.map_to_field(config_ptr);
        assert_eq!(i512_neg_result.into_bigint(), BigInt::from(18u64)); // 23 - 5 = 18

        // Test reference implementation
        let i512_ref_result = (&i512_val).map_to_field(config_ptr);
        assert_eq!(i512_ref_result.into_bigint(), BigInt::from(5u64));
    }

    #[test]
    #[should_panic(expected = "Cannot convert signed integer to prime field element without a modulus")]
    fn test_signed_field_map_null_config() {
        let i32_val: i32 = 5;
        i32_val.map_to_field::<1>(std::ptr::null());
    }

    #[test]
    #[should_panic(expected = "Cannot convert unsigned integer to prime field element without a modulus")]
    fn test_unsigned_field_map_null_config() {
        let u32_val: u32 = 5;
        u32_val.map_to_field::<1>(std::ptr::null());
    }
}
