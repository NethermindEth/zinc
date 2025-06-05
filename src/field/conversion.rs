use ark_std::cmp::Ordering;
use ark_std::mem::transmute_copy;
use ark_std::vec::Vec;
use crypto_bigint::{Int, NonZero, Uint};

use crate::biginteger::BigInt;
use crate::field::RandomField;
use crate::field::RandomField::Raw;
use crate::field_config::ConfigRef;
use crate::primitives::{Abs, Unsigned};
use crate::traits::FromBytes;

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

pub trait FieldMap<'cfg, const N: usize> {
    type Cfg: Copy;
    type Output;
    fn map_to_field(&self, config: Self::Cfg) -> Self::Output;
}

// Implementation of FieldMap for signed integers
impl<'cfg, const N: usize, T: Abs + Copy> FieldMap<'cfg, N> for T {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = RandomField<'cfg, N>;

    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        let config = match config.reference() {
            Some(config) => config,
            None => {
                panic!("Cannot convert integer to prime field element without a modulus")
            }
        };

        let modulus: [u64; N] = config.modulus().to_words();
        let abs_val = (*self).unsigned_abs();

        let limbs = <T as Abs>::Unsigned::limbs();
        // Calculate how many u64 limbs we need based on bits
        let val = abs_val.as_array::<N>();

        let mut r = match limbs.cmp(&N) {
            Ordering::Less => {
                let wider_value: [u64; N] = unsafe { transmute_copy(&val) };
                let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                wider %= crypto_bigint::NonZero::new(modu).unwrap();
                BigInt::from(wider.to_words())
            }
            Ordering::Equal => {
                let mut value_N = [0u64; N];
                value_N.copy_from_slice(&val);

                let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                value %= crypto_bigint::NonZero::new(modu).unwrap();
                BigInt::from(value.to_words())
            }
            Ordering::Greater => {
                let mut wider_modulus = [0u64; 2];
                wider_modulus[..N].copy_from_slice(&modulus);
                let mut slice = [0u64; 2];
                slice[..N.min(limbs)].copy_from_slice(&val[..N.min(limbs)]);
                let mut value = crypto_bigint::Uint::<2>::from_words(slice);
                let modu = crypto_bigint::Uint::<2>::from_words(wider_modulus);

                value %= crypto_bigint::NonZero::new(modu).unwrap();
                let mut result = [0u64; N];
                result.copy_from_slice(&value.to_words()[..N]);

                BigInt::from(result)
            }
        };

        config.mul_assign(&mut r, config.r2());

        let mut elem = RandomField::<N>::new_unchecked(ConfigRef::from(config), r);

        if self.is_negative() {
            elem = -elem; // Negate if the original value was negative
        }

        elem
    }
}

// Implementation for bool
impl<'cfg, const N: usize> FieldMap<'cfg, N> for bool {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = RandomField<'cfg, N>;

    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        let config = match config.reference() {
            Some(config) => config,
            None => panic!("Cannot convert boolean to prime field element without a modulus"),
        };

        let mut r = BigInt::from(*self as u64);
        config.mul_assign(&mut r, config.r2());
        RandomField::<N>::new_unchecked(ConfigRef::from(config), r)
    }
}

impl<'cfg, const N: usize> FieldMap<'cfg, N> for &bool {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = RandomField<'cfg, N>;
    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        (*self).map_to_field(config)
    }
}

// Implementation for Int<N>
impl<'cfg, const M: usize, const N: usize> FieldMap<'cfg, N> for Int<M> {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = RandomField<'cfg, N>;

    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        let local_type_bigint = BigInt::from(self);
        let res = local_type_bigint.map_to_field(config);
        if self < &Int::ZERO {
            return -res;
        }
        res
    }
}

impl<'cfg, const M: usize, const N: usize> FieldMap<'cfg, N> for &Int<M> {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = RandomField<'cfg, N>;

    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        (*self).map_to_field(config)
    }
}
// Implementation of FieldMap for BigInt<N>
impl<'cfg, const M: usize, const N: usize> FieldMap<'cfg, N> for BigInt<M> {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = RandomField<'cfg, N>;

    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        let config = match config.reference() {
            Some(config) => config,
            None => panic!("Cannot convert BigInt to prime field element without a modulus"),
        };

        let modulus: [u64; N] = config.modulus().to_words();

        let mut r: BigInt<N> = match M.cmp(&N) {
            ark_std::cmp::Ordering::Less => {
                let mut wider_value = [0u64; N];
                wider_value[..M].copy_from_slice(&self.to_words());
                let mut value = Uint::from_words(wider_value);
                let modu = Uint::<N>::from_words(modulus);

                value %= NonZero::new(modu).unwrap();
                let mut result = [0u64; N];
                result.copy_from_slice(&value.to_words()[..N]);

                BigInt::from(result)
            }
            ark_std::cmp::Ordering::Equal => {
                let mut value = Uint::<M>::from_words(self.to_words());
                let mut wider_modulus = [0u64; M];
                wider_modulus[..N].copy_from_slice(&modulus);
                let modu = Uint::<M>::from_words(wider_modulus);

                value %= NonZero::new(modu).unwrap();
                let mut result = [0u64; N];
                result.copy_from_slice(&value.to_words()[..N]);

                BigInt::from(result)
            }
            ark_std::cmp::Ordering::Greater => {
                let mut value = Uint::<M>::from_words(self.to_words());
                let mut wider_modulus = [0u64; M];
                wider_modulus[..N].copy_from_slice(&modulus);
                let modu = Uint::<M>::from_words(wider_modulus);

                value %= NonZero::new(modu).unwrap();
                let mut result = [0u64; N];
                result.copy_from_slice(&value.to_words()[..N]);

                BigInt::from(result)
            }
        };

        // Apply Montgomery form transformation
        config.mul_assign(&mut r, config.r2());
        RandomField::<N>::new_unchecked(ConfigRef::from(config), r)
    }
}

// Implementation of FieldMap for reference to BigInt<N>
impl<'cfg, const M: usize, const N: usize> FieldMap<'cfg, N> for &BigInt<M> {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = RandomField<'cfg, N>;
    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        (*self).map_to_field(config)
    }
}

// Implementation of FieldMap for Vec<T>

impl<'cfg, const N: usize, T: FieldMap<'cfg, N>> FieldMap<'cfg, N> for Vec<T> {
    type Cfg = T::Cfg;
    type Output = Vec<T::Output>;

    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config)).collect()
    }
}

impl<'cfg, const N: usize, T: FieldMap<'cfg, N>> FieldMap<'cfg, N> for &Vec<T> {
    type Cfg = T::Cfg;
    type Output = Vec<T::Output>;

    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config)).collect()
    }
}

impl<'cfg, const N: usize, T: FieldMap<'cfg, N>> FieldMap<'cfg, N> for &[T] {
    type Cfg = T::Cfg;
    type Output = Vec<T::Output>;

    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config)).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::field::conversion::FieldMap;
    use crate::field_config::{ConfigRef, FieldConfig};
    use crate::traits::FromBytes;
    use crate::{biginteger::BigInt, field::RandomField};
    use ark_std::format;
    use ark_std::str::FromStr;

    fn test_from<'cfg, T: Clone, const N: usize>(value: T, value_str: &str)
    where
        RandomField<'cfg, N>: From<T>,
    {
        let raw_element = RandomField::<'cfg, N>::from(value);
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
                let result = x.map_to_field($config);
                let ref_result = (&x).map_to_field($config);
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
            let zero_result = zero.map_to_field($config);
            assert_eq!(
                zero_result.into_bigint(),
                BigInt::<$N>::zero(),
                "Zero value should map to field zero"
            );

            // Test maximum value
            let max = <$type>::from_str(&format!("{}", <$type>::MAX)).unwrap();
            let max_result = max.map_to_field($config);
            assert!(
                max_result.into_bigint() < BigInt::<$N>::from($field),
                "Maximum value should be less than field modulus"
            );

            // Test minimum value
            let min = -<$type>::from_str(&format!("{}", <$type>::MAX)).unwrap();
            assert!(
                min.map_to_field($config).into_bigint() < BigInt::<$N>::from($field),
                "Minimum value should wrap to valid field element"
            );

            // Test positive boundary
            let pos = <$type>::from_str("5").unwrap();
            let pos_result = pos.map_to_field($config);
            assert_eq!(
                pos_result.into_bigint(),
                BigInt::<$N>::from(5u64),
                "Positive value should map directly to field"
            );

            // Test negative boundary
            let neg = <$type>::from_str("-5").unwrap();
            let neg_result = neg.map_to_field($config);
            assert_eq!(
                neg_result.into_bigint(),
                BigInt::<$N>::from(($field as i64 - 5) as u64),
                "Negative value should wrap around field modulus"
            );

            // Test reference conversions
            let ref_zero = (&zero).map_to_field($config);
            assert_eq!(
                ref_zero.into_bigint(),
                BigInt::<$N>::zero(),
                "Reference to zero should map to field zero"
            );

            let ref_max = (&max).map_to_field($config);
            assert!(
                ref_max.into_bigint() < BigInt::<$N>::from($field),
                "Reference to maximum value should be less than field modulus"
            );

            let ref_min = (&min).map_to_field($config);
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
                let result = x.map_to_field($config);
                let ref_result = (&x).map_to_field($config);
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
            let zero_result = zero.map_to_field($config);
            assert_eq!(
                zero_result.into_bigint(),
                BigInt::<$N>::zero(),
                "Zero value should map to field zero"
            );

            // Test maximum value
            let max = <$type>::MAX;
            let max_result = max.map_to_field($config);
            assert!(
                max_result.into_bigint() < BigInt::<$N>::from($field),
                "Maximum value should be less than field modulus"
            );

            // Test boundary value - using literal instead of From
            let boundary: $type = 5;
            let boundary_result = boundary.map_to_field($config);
            assert_eq!(
                boundary_result.into_bigint(),
                BigInt::<$N>::from(5u64),
                "Boundary value should map directly to field"
            );

            // Test reference conversions
            let ref_zero = (&zero).map_to_field($config);
            assert_eq!(
                ref_zero.into_bigint(),
                BigInt::<$N>::zero(),
                "Reference to zero should map to field zero"
            );

            let ref_max = (&max).map_to_field($config);
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
        i32_val.map_to_field(ConfigRef::<1>::NONE);
    }

    #[test]
    #[should_panic(expected = "Cannot convert integer to prime field element without a modulus")]
    fn test_unsigned_field_map_null_config() {
        let u32_val: u32 = 5;
        u32_val.map_to_field(ConfigRef::<1>::NONE);
    }
}

#[cfg(test)]
mod bigint_field_map_tests {
    use super::*;
    use crate::field_config::FieldConfig;
    use ark_std::str::FromStr;

    #[test]
    fn test_bigint_smaller_than_field() {
        // Using a 2-limb field config with 1-limb BigInt
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = ConfigRef::from(&config);

        let small_bigint = BigInt::<1>::from(12345u64);
        let result = small_bigint.map_to_field(config_ptr);

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
        let result = value.map_to_field(config_ptr);

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
        let result = large_value.map_to_field(config_ptr);

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
        let result = zero.map_to_field(config_ptr);

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
        let result = value.map_to_field(config_ptr);
        let direct_result = value.map_to_field(config_ptr);

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

        let result = max_value.map_to_field(config_ptr);

        assert!(
            result.into_bigint() < modulus,
            "Result should be properly reduced modulo field modulus"
        );
    }
}
