use ark_std::{vec::Vec, Zero};

use crate::{
    biginteger::BigInt,
    crypto_int::CryptoInt,
    field::{RandomField, RandomField::Raw},
    field_config::ConfigRef,
    traits::{
        Config, ConfigReference, CryptoInteger, CryptoUinteger, Field, FieldMap, FromBytes,
        Integer, Words,
    },
};

impl<const N: usize> From<u128> for RandomField<'_, N> {
    fn from(value: u128) -> Self {
        let value = BigInt::from(value);

        Raw { value }
    }
}

macro_rules! impl_from_uint {
    ($type:ty) => {
        impl<const N: usize> From<$type> for RandomField<'_, N> {
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

impl<const N: usize> From<bool> for RandomField<'_, N> {
    fn from(value: bool) -> Self {
        let value = BigInt::from(value as u8);
        Raw { value }
    }
}

impl<const N: usize> FromBytes for RandomField<'_, N> {
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

impl<'cfg, const N: usize> RandomField<'cfg, N> {
    pub fn from_bytes_le_with_config(config: ConfigRef<'cfg, N>, bytes: &[u8]) -> Option<Self> {
        let value = BigInt::<N>::from_bytes_le(bytes);

        Self::from_bigint(config, value?)
    }

    pub fn from_bytes_be_with_config(config: ConfigRef<'cfg, N>, bytes: &[u8]) -> Option<Self> {
        let value = BigInt::<N>::from_bytes_be(bytes);

        Self::from_bigint(config, value?)
    }
}

// Implementation of FieldMap for signed integers
macro_rules! impl_field_map_for_int {
    ($t:ty) => {
        impl<F: Field> FieldMap<F> for $t {
            type Output = F;

            fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
                let config = match config_ref.reference() {
                    Some(config) => config,
                    None => {
                        panic!("Cannot convert integer to prime field element without a modulus")
                    }
                };
                let value = self.abs_diff(0);
                let mut words = F::W::default();

                words[0] = value as u64;

                if (ark_std::mem::size_of::<$t>() + 7) / 8 > 1 && F::W::num_words() > 1 {
                    words[1] = ((value as u128) >> 64) as u64;
                }
                let mut value = F::CryptoUint::from_words(words).as_int();
                let modulus = F::CryptoInt::from_words(config.modulus().to_words());

                value %= modulus;
                let mut value = F::I::from(value);
                config.mul_assign(&mut value, config.r2());

                let mut r = F::new_unchecked(config_ref, value);

                if self < &<$t>::zero() {
                    r = -r;
                }

                r
            }
        }
    };
}

impl_field_map_for_int!(i8);
impl_field_map_for_int!(u8);
impl_field_map_for_int!(i16);
impl_field_map_for_int!(u16);
impl_field_map_for_int!(i32);
impl_field_map_for_int!(u32);
impl_field_map_for_int!(i64);
impl_field_map_for_int!(u64);
impl_field_map_for_int!(i128);
impl_field_map_for_int!(u128);
impl_field_map_for_int!(isize);
impl_field_map_for_int!(usize);

// Implementation for bool
impl<F: Field> FieldMap<F> for bool {
    type Output = F;

    fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
        let config = match config_ref.reference() {
            Some(config) => config,
            None => panic!("Cannot convert boolean to prime field element without a modulus"),
        };

        let mut r = F::I::from(*self as u64);
        config.mul_assign(&mut r, config.r2());
        F::new_unchecked(config_ref, r)
    }
}

impl<F: Field> FieldMap<F> for &bool {
    type Output = F;
    fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
        (*self).map_to_field(config_ref)
    }
}

// Implementation for Int<N>
impl<F: Field, T: CryptoInteger> FieldMap<F> for T
where
    T::I: FieldMap<F, Output = F>,
{
    type Output = F;

    fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
        let local_type_bigint = T::I::from(self);
        let res = local_type_bigint.map_to_field(config_ref);
        if self < &<T as Zero>::zero() {
            return -res;
        }
        res
    }
}

// Implementation of FieldMap for BigInt<N>
impl<F: Field, const M: usize> FieldMap<F> for BigInt<M>
where
    for<'a> CryptoInt<M>: From<&'a F::I>,
    F::I: From<CryptoInt<M>>,
    for<'a> F::CryptoInt: From<&'a BigInt<M>>,
{
    type Output = F;

    fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
        let config = match config_ref.reference() {
            Some(config) => config,
            None => panic!("Cannot convert BigInt to prime field element without a modulus"),
        };

        let mut value = if M > F::W::num_words() {
            let modulus: CryptoInt<M> = config.modulus().into();
            let mut value: CryptoInt<M> = self.into();
            value %= modulus;

            F::I::from(value)
        } else {
            let modulus: F::CryptoInt = config.modulus().into();
            let mut value: F::CryptoInt = self.into();
            value %= modulus;

            value.into()
        };

        config.mul_assign(&mut value, config.r2());

        F::new_unchecked(config_ref, value)
    }
}

// Implementation of FieldMap for reference to BigInt<N>
impl<F: Field, const M: usize> FieldMap<F> for &BigInt<M>
where
    BigInt<M>: FieldMap<F, Output = F>,
{
    type Output = F;
    fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
        (*self).map_to_field(config_ref)
    }
}

// Implementation of FieldMap for Vec<T>

impl<F: Field, T: FieldMap<F>> FieldMap<F> for Vec<T> {
    type Output = Vec<T::Output>;

    fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config_ref)).collect()
    }
}

impl<F: Field, T: FieldMap<F>> FieldMap<F> for &Vec<T> {
    type Output = Vec<T::Output>;

    fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config_ref)).collect()
    }
}

impl<F: Field, T: FieldMap<F>> FieldMap<F> for &[T] {
    type Output = Vec<T::Output>;

    fn map_to_field(&self, config_ref: F::Cr) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config_ref)).collect()
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{fmt::Debug, format, str::FromStr};

    use crate::{
        biginteger::BigInt,
        field::RandomField,
        field_config::{ConfigRef, FieldConfig},
        traits::{Config, ConfigReference, Field, FieldMap, FromBytes},
    };

    fn test_from<F: Field + From<T>, T: Clone>(value: T, value_str: &str)
    where
        F::I: FromStr,
        <F::I as FromStr>::Err: Debug,
    {
        let raw_element = F::from(value);
        assert_eq!(
            raw_element,
            F::without_config(F::I::from_str(value_str).unwrap())
        )
    }

    #[test]
    fn converts_u128_to_random_field() {
        test_from::<RandomField<2>, u128>(
            243043087159742188419721163456177516,
            "243043087159742188419721163456177516",
        );
    }

    #[test]
    #[should_panic(expected = "Integer is 128 bits but field is 64 bits")]
    fn panics_when_u128_does_not_fit_in_n1() {
        test_from::<RandomField<1>, u128>(243043087159742188419721163456177516, "");
    }

    #[test]
    fn converts_u64_to_random_field() {
        test_from::<RandomField<1>, u64>(23, "23");
    }

    #[test]
    fn converts_u32_to_random_field() {
        test_from::<RandomField<1>, u32>(23, "23");
    }

    #[test]
    fn converts_u16_to_random_field() {
        test_from::<RandomField<1>, u16>(23, "23");
    }

    #[test]
    fn converts_u8_to_random_field() {
        test_from::<RandomField<1>, u8>(23, "23");
    }

    #[test]
    fn converts_false_to_zero() {
        test_from::<RandomField<1>, bool>(false, "0");
    }

    #[test]
    fn converts_true_to_one() {
        test_from::<RandomField<1>, bool>(true, "1");
    }

    #[test]
    fn converts_from_bytes_le_with_config_valid() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0x05, 0, 0, 0, 0, 0, 0, 0];
        let expected = BigInt::from_str("5").unwrap();

        let result = RandomField::<1>::from_bytes_le_with_config(config, &bytes).unwrap();
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
        let config_ptr = ConfigRef::from(&config);

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
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0x00; 8]; // All zeros
        let expected = RandomField::Initialized {
            config,
            value: BigInt::<1>::zero(),
        };

        let result = RandomField::<1>::from_bytes_le_with_config(config, &bytes);
        assert_eq!(result, Some(expected));
    }

    #[test]
    fn converts_from_bytes_be_with_config_zero() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0x00; 8]; // All zeros
        let expected = RandomField::Initialized {
            config,
            value: BigInt::<1>::zero(),
        };

        let result = RandomField::<1>::from_bytes_be_with_config(config, &bytes);
        assert_eq!(result, Some(expected));
    }

    #[test]
    fn converts_from_bytes_le_with_config_out_of_range() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0x65, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // Value: 101 (modulus is 23)
        let result = RandomField::<1>::from_bytes_le_with_config(config, &bytes);
        assert!(result.is_none());
    }

    #[test]
    fn converts_from_bytes_be_with_config_out_of_range() {
        let config = FieldConfig::new(BigInt::<32>::from_str("37129241769965749").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0x65, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // Value: 101
        let result = RandomField::<32>::from_bytes_be_with_config(config, &bytes);
        assert!(result.is_none());
    }

    #[test]
    fn converts_from_bytes_le_with_config_exact_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0x17, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // Value: 23 (modulus)
        let result = RandomField::<1>::from_bytes_le_with_config(config, &bytes);
        assert!(result.is_none()); // Must be strictly less than modulus
    }

    #[test]
    fn converts_from_bytes_be_with_config_exact_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x17]; // Value: 23 (big-endian)
        let result = RandomField::<1>::from_bytes_be_with_config(config, &bytes);
        assert!(result.is_none());
    }

    #[test]
    fn converts_from_bytes_le_with_config_leading_zeros() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0b0000_0001]; // Value: 1 with leading zeros
        let expected = BigInt::<1>::from_bytes_le(&bytes)
            .unwrap()
            .map_to_field(config);

        let result = RandomField::<1>::from_bytes_le_with_config(config, &bytes);
        assert_eq!(result, Some(expected));
    }

    #[test]
    fn converts_from_bytes_be_with_config_leading_zeros() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let bytes = [0x01]; //1 with leading zeros (big-endian);

        let result = RandomField::<1>::from_bytes_be_with_config(config, &bytes).unwrap();
        assert_eq!(result.into_bigint(), BigInt::one());
    }

    macro_rules! test_signed_type_full_range {
        ($type:ty, $field:expr, $config:expr, $N:expr) => {{
            // Test full range for primitive types
            for x in <$type>::MIN..=<$type>::MAX {
                let result: RandomField<$N> = x.map_to_field($config);
                let ref_result: RandomField<$N> = (&x).map_to_field($config);
                let expected = if x < 0 {
                    BigInt::<$N>::from(($field as i64 + x as i64) as u64)
                } else {
                    BigInt::<$N>::from(x as u64)
                };
                assert_eq!(
                    result.into_bigint(),
                    expected,
                    "conversion failed for value: {}",
                    x
                );
                assert_eq!(
                    ref_result.into_bigint(),
                    expected,
                    "reference conversion failed for value: {}",
                    x
                );
            }
        }};
    }

    macro_rules! test_signed_type_edge_cases {
        ($type:ty, $field:expr, $config:expr, $N:expr) => {{
            // Test zero
            let zero = <$type>::from_str("0").unwrap();
            let zero_result: RandomField<$N> = zero.map_to_field($config);
            assert_eq!(
                zero_result.into_bigint(),
                BigInt::<$N>::zero(),
                "Zero value should map to field zero"
            );

            // Test maximum value
            let max = <$type>::from_str(&format!("{}", <$type>::MAX)).unwrap();
            let max_result: RandomField<$N> = max.map_to_field($config);
            assert!(
                max_result.into_bigint() < BigInt::<$N>::from($field),
                "Maximum value should be less than field modulus"
            );

            // Test minimum value
            let min = -<$type>::from_str(&format!("{}", <$type>::MAX)).unwrap();
            let min_f = <$type as FieldMap<RandomField<$N>>>::map_to_field(&min, $config);
            assert!(
                min_f.into_bigint() < BigInt::<$N>::from($field),
                "Minimum value should wrap to valid field element"
            );

            // Test positive boundary
            let pos = <$type>::from_str("5").unwrap();
            let pos_result: RandomField<$N> = pos.map_to_field($config);
            assert_eq!(
                pos_result.into_bigint(),
                BigInt::<$N>::from(5u64),
                "Positive value should map directly to field"
            );

            // Test negative boundary
            let neg = <$type>::from_str("-5").unwrap();
            let neg_result: RandomField<$N> = neg.map_to_field($config);
            assert_eq!(
                neg_result.into_bigint(),
                BigInt::<$N>::from(($field as i64 - 5) as u64),
                "Negative value should wrap around field modulus"
            );

            // Test reference conversions
            let ref_zero: RandomField<$N> = (&zero).map_to_field($config);
            assert_eq!(
                ref_zero.into_bigint(),
                BigInt::<$N>::zero(),
                "Reference to zero should map to field zero"
            );

            let ref_max: RandomField<$N> = (&max).map_to_field($config);
            assert!(
                ref_max.into_bigint() < BigInt::<$N>::from($field),
                "Reference to maximum value should be less than field modulus"
            );

            let ref_min: RandomField<$N> = (&min).map_to_field($config);
            assert!(
                ref_min.into_bigint() < BigInt::<$N>::from($field),
                "Reference to minimum value should wrap to valid field element"
            );
        }};
    }

    #[test]
    fn test_signed_integers_field_map() {
        let field = 18446744069414584321_u64;
        let config = FieldConfig::new(BigInt::from_str("18446744069414584321").unwrap());
        let config: ConfigRef<1> = ConfigRef::from(&config);

        // Test primitive types with full range
        test_signed_type_full_range!(i8, field, config, 1);
        test_signed_type_full_range!(i16, field, config, 1);

        // Test larger primitive types with edge cases only
        test_signed_type_edge_cases!(i32, field, config, 1);
        test_signed_type_edge_cases!(i64, field, config, 1);
        test_signed_type_edge_cases!(i128, field, config, 1);
    }

    macro_rules! test_unsigned_type_full_range {
        ($type:ty, $field:expr, $config:expr, $N:expr) => {{
            // Test full range for small unsigned types
            for x in <$type>::MIN..=<$type>::MAX {
                let result: RandomField<$N> = x.map_to_field($config);
                let ref_result: RandomField<$N> = (&x).map_to_field($config);
                let expected = BigInt::<$N>::from(x as u64);
                assert_eq!(
                    result.into_bigint(),
                    expected,
                    "conversion failed for value: {}",
                    x
                );
                assert_eq!(
                    ref_result.into_bigint(),
                    expected,
                    "reference conversion failed for value: {}",
                    x
                );
            }
        }};
    }

    macro_rules! test_unsigned_type_edge_cases {
        ($type:ty, $field:expr, $config:expr, $N:expr) => {{
            // Test zero
            let zero = <$type>::MIN;
            let zero_result: RandomField<$N> = zero.map_to_field($config);
            assert_eq!(
                zero_result.into_bigint(),
                BigInt::<$N>::zero(),
                "Zero value should map to field zero"
            );

            // Test maximum value
            let max = <$type>::MAX;
            let max_result: RandomField<$N> = max.map_to_field($config);
            assert!(
                max_result.into_bigint() < BigInt::<$N>::from($field),
                "Maximum value should be less than field modulus"
            );

            // Test boundary value - using literal instead of From
            let boundary: $type = 5;
            let boundary_result: RandomField<$N> = boundary.map_to_field($config);
            assert_eq!(
                boundary_result.into_bigint(),
                BigInt::<$N>::from(5u64),
                "Boundary value should map directly to field"
            );

            // Test reference conversions
            let ref_zero: RandomField<$N> = (&zero).map_to_field($config);
            assert_eq!(
                ref_zero.into_bigint(),
                BigInt::<$N>::zero(),
                "Reference to zero should map to field zero"
            );

            let ref_max: RandomField<$N> = (&max).map_to_field($config);
            assert!(
                ref_max.into_bigint() < BigInt::<$N>::from($field),
                "Reference to maximum value should be less than field modulus"
            );
        }};
    }

    #[test]
    fn test_unsigned_integers_field_map() {
        let field_1 = 18446744069414584321_u64;
        let config_1 = FieldConfig::new(BigInt::from_str("18446744069414584321").unwrap());
        let config = ConfigRef::from(&config_1);
        // Test small types with full range
        test_unsigned_type_full_range!(u8, field_1, config, 1);
        test_unsigned_type_full_range!(u16, field_1, config, 1);

        // Test larger types with edge cases only
        test_unsigned_type_edge_cases!(u32, field_1, config, 1);
        test_unsigned_type_edge_cases!(u64, field_1, config, 1);
        test_unsigned_type_edge_cases!(u128, field_1, config, 1);
    }

    #[test]
    #[should_panic(expected = "Cannot convert integer to prime field element without a modulus")]
    fn test_signed_field_map_null_config() {
        let i32_val: i32 = 5;
        <i32 as FieldMap<RandomField<1>>>::map_to_field(&i32_val, ConfigRef::<1>::NONE);
    }

    #[test]
    #[should_panic(expected = "Cannot convert integer to prime field element without a modulus")]
    fn test_unsigned_field_map_null_config() {
        let u32_val: u32 = 5;
        <u32 as FieldMap<RandomField<1>>>::map_to_field(&u32_val, ConfigRef::<1>::NONE);
    }
}

#[cfg(test)]
mod bigint_field_map_tests {
    use ark_std::str::FromStr;

    use super::*;
    use crate::field_config::FieldConfig;

    #[test]
    fn test_bigint_smaller_than_field() {
        // Using a 2-limb field config with 1-limb BigInt
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = ConfigRef::from(&config);

        let small_bigint = BigInt::<1>::from(12345u64);
        let result: RandomField<2> = small_bigint.map_to_field(config_ptr);

        assert_eq!(
            result.into_bigint().first(),
            12345u64,
            "Small BigInt should be preserved in larger field"
        );
    }

    #[test]
    fn test_bigint_equal_size() {
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = ConfigRef::from(&config);

        let value = BigInt::<2>::from_str("12345678901234567890").unwrap();
        let result: RandomField<2> = value.map_to_field(config_ptr);

        // The result should be the value modulo the field modulus
        let expected = BigInt::<2>::from_str("12345678901234567890").unwrap();
        assert_eq!(
            result.into_bigint(),
            expected,
            "Equal size BigInt should be correctly converted"
        );
    }

    #[test]
    fn test_bigint_larger_than_field() {
        // Using a 1-limb field config with 2-limb BigInt
        let modulus = BigInt::<1>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = ConfigRef::from(&config);

        let large_value = BigInt::<2>::from_str("123456789012345678901").unwrap();
        let result: RandomField<1> = large_value.map_to_field(config_ptr);

        let expected = BigInt::<1>::from(12776324595858172975u64);
        assert_eq!(
            result.into_bigint(),
            expected,
            "Larger BigInt should be correctly reduced modulo field modulus"
        );
    }

    #[test]
    fn test_bigint_zero() {
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = ConfigRef::from(&config);

        let zero = BigInt::<2>::zero();
        let result: RandomField<2> = zero.map_to_field(config_ptr);

        assert!(
            result.into_bigint().is_zero(),
            "Zero BigInt should map to zero field element"
        );
    }

    #[test]
    fn test_bigint_reference() {
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = ConfigRef::from(&config);

        let value = BigInt::<2>::from_str("12345").unwrap();
        let result: RandomField<2> = value.map_to_field(config_ptr);
        let direct_result: RandomField<2> = value.map_to_field(config_ptr);

        assert_eq!(
            result.into_bigint(),
            direct_result.into_bigint(),
            "Reference implementation should match direct implementation"
        );
    }

    #[test]
    #[should_panic(expected = "Cannot convert BigInt to prime field element without a modulus")]
    fn test_null_config() {
        let value = BigInt::<2>::from(123u64);
        let _result: RandomField<2> = value.map_to_field(ConfigRef::NONE);
    }

    #[test]
    fn test_bigint_max_value() {
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = ConfigRef::from(&config);

        // Create a BigInt with all bits set to 1
        let max_value = BigInt::from([u64::MAX, u64::MAX]);

        let result: RandomField<2> = max_value.map_to_field(config_ptr);

        assert!(
            result.into_bigint() < modulus,
            "Result should be properly reduced modulo field modulus"
        );
    }
}
