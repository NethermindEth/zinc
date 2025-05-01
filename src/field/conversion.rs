use crypto_bigint::Int;
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

pub trait FieldMap<const N: usize> {
    type Output;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output;
}

macro_rules! impl_field_map_for_int {
    ($type:ty, $bits:expr) => {
        impl<const N: usize> FieldMap<N> for $type {
            type Output = RandomField<N>;
            fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
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

                            value %= crypto_bigint::NonZero::new(modu).unwrap();
                            let mut result = [0u64; N];
                            result.copy_from_slice(&value.to_words()[..N]);

                            BigInt(result)
                        }
                        n if n == LIMBS => {
                            let mut value_N = [0u64; N];
                            value_N.copy_from_slice(&val);

                            let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                            let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                            value %= crypto_bigint::NonZero::new(modu).unwrap();
                            value %= crypto_bigint::NonZero::new(modu).unwrap();
                            BigInt(value.to_words())
                        }
                        _ => {
                            let mut wider_value = [0u64; N];
                            wider_value[..LIMBS].copy_from_slice(&val);
                            let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                            let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                            wider %= crypto_bigint::NonZero::new(modu).unwrap();
                            wider %= crypto_bigint::NonZero::new(modu).unwrap();
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

        impl<const N: usize> FieldMap<N> for &$type {
            type Output = RandomField<N>;
            fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
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
impl<const N: usize> FieldMap<N> for I256 {
    type Output = RandomField<N>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        if config.is_null() {
            panic!("Cannot convert signed integer to prime field element without a modulus")
        }
        unsafe {
            let modulus: [u64; N] = (*config).modulus.0;
            if *self == I256::MIN {
                panic!("Cannot convert I256::MIN to field element due to overflow in abs()");
            }
            let val: [u64; 4] = self.abs().to_le_u64();

            let mut r: BigInt<N> = match N {
                n if n < 4 => {
                    let mut wider_modulus: [u64; 4] = [0; 4];
                    wider_modulus[..N].copy_from_slice(&modulus);
                    let value = crypto_bigint::Uint::<4>::from_words(val);
                    let modu = crypto_bigint::Uint::<4>::from_words(wider_modulus);

                    let value = if value >= modu {
                        value % crypto_bigint::NonZero::new(modu).unwrap()
                    } else {
                        value
                    };
                    let mut result = [0u64; N];
                    result.copy_from_slice(&value.to_words()[..N]);

                    BigInt(result)
                }
                4 => {
                    let mut value_N: [u64; N] = [0; N];
                    value_N.copy_from_slice(&val);

                    let value = crypto_bigint::Uint::<N>::from_words(value_N);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    let value = if value >= modu {
                        value % crypto_bigint::NonZero::new(modu).unwrap()
                    } else {
                        value
                    };
                    BigInt(value.to_words())
                }
                _ => {
                    let mut wider_value: [u64; N] = [0; N];
                    wider_value[..4].copy_from_slice(&val);
                    let wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    let wider = if wider >= modu {
                        wider % crypto_bigint::NonZero::new(modu).unwrap()
                    } else {
                        wider
                    };
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
impl<const N: usize> FieldMap<N> for &I256 {
    type Output = RandomField<N>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        (*self).map_to_field(config)
    }
}

// Implementation for I512
impl<const N: usize> FieldMap<N> for I512 {
    type Output = RandomField<N>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
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

                    value %= crypto_bigint::NonZero::new(modu).unwrap();
                    let mut result = [0u64; N];
                    result.copy_from_slice(&value.to_words()[..N]);

                    BigInt(result)
                }
                8 => {
                    let mut value_N: [u64; N] = [0; N];
                    value_N.copy_from_slice(&val);

                    let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    value %= crypto_bigint::NonZero::new(modu).unwrap();
                    BigInt(value.to_words())
                }
                _ => {
                    let mut wider_value: [u64; N] = [0; N];
                    wider_value[..8].copy_from_slice(&val);
                    let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    wider %= crypto_bigint::NonZero::new(modu).unwrap();
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

impl<const N: usize> FieldMap<N> for &I512 {
    type Output = RandomField<N>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        (*self).map_to_field(config)
    }
}

macro_rules! impl_field_map_for_vec {
    ($type:ty) => {
        impl<const N: usize> FieldMap<N> for Vec<$type> {
            type Output = Vec<RandomField<N>>;
            fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
                self.iter().map(|x| x.map_to_field(config)).collect()
            }
        }

        impl<const N: usize> FieldMap<N> for &Vec<$type> {
            type Output = Vec<RandomField<N>>;
            fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
                self.iter().map(|x| x.map_to_field(config)).collect()
            }
        }

        impl<const M: usize> FieldMap<M> for &[$type] {
            type Output = Vec<RandomField<M>>;
            fn map_to_field(&self, config: *const FieldConfig<M>) -> Self::Output {
                self.iter().map(|x| x.map_to_field(config)).collect()
            }
        }
    };
}

impl_field_map_for_vec!(i8);
impl_field_map_for_vec!(i16);
impl_field_map_for_vec!(i32);
impl_field_map_for_vec!(i64);
impl_field_map_for_vec!(i128);
impl_field_map_for_vec!(I256);
impl_field_map_for_vec!(I512);

macro_rules! impl_field_map_for_uint {
    ($type:ty, $bits:expr) => {
        impl<const N: usize> FieldMap<N> for $type {
            type Output = RandomField<N>;
            fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
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
                            value %= crypto_bigint::NonZero::new(modu).unwrap();
                            let mut result = [0u64; N];
                            result.copy_from_slice(&value.to_words()[..N]);

                            BigInt(result)
                        }
                        n if n == LIMBS => {
                            let mut value_N = [0u64; N];
                            value_N.copy_from_slice(&val);

                            let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                            let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                            value %= crypto_bigint::NonZero::new(modu).unwrap();
                            BigInt(value.to_words())
                        }
                        _ => {
                            let mut wider_value = [0u64; N];
                            wider_value[..LIMBS].copy_from_slice(&val);
                            let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                            let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                            wider %= crypto_bigint::NonZero::new(modu).unwrap();
                            BigInt(wider.to_words())
                        }
                    };

                    (*config).mul_assign(&mut r, &(*config).r2);
                    RandomField::<N>::new_unchecked(config, r)
                }
            }
        }

        impl<const N: usize> FieldMap<N> for &$type {
            type Output = RandomField<N>;
            fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
                (*self).map_to_field(config)
            }
        }
    };
}

// Usage with bit sizes
impl_field_map_for_uint!(u8, 8);
impl_field_map_for_uint!(u16, 16);
impl_field_map_for_uint!(u32, 32);
impl_field_map_for_uint!(u64, 64);
impl_field_map_for_uint!(u128, 128);

// Implementation for bool
impl<const N: usize> FieldMap<N> for bool {
    type Output = RandomField<N>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
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

impl<const N: usize> FieldMap<N> for &bool {
    type Output = RandomField<N>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        (*self).map_to_field(config)
    }
}

// Implementation for Int<N>
impl<const M: usize, const N: usize> FieldMap<N> for Int<M> {
    type Output = RandomField<N>;

    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        let local_type_bigint = BigInt::new(*self.as_words());
        local_type_bigint.map_to_field(config)
    }
}

impl<const M: usize, const N: usize> FieldMap<N> for &Int<M> {
    type Output = RandomField<N>;

    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        (*self).map_to_field(config)
    }
}
impl<const N: usize, const M: usize> FieldMap<N> for Vec<Int<M>> {
    type Output = Vec<RandomField<N>>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config)).collect()
    }
}

impl<const N: usize, const M: usize> FieldMap<N> for &Vec<Int<M>> {
    type Output = Vec<RandomField<N>>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config)).collect()
    }
}

impl<const N: usize, const M: usize> FieldMap<N> for &[Int<M>] {
    type Output = Vec<RandomField<N>>;
    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        self.iter().map(|x| x.map_to_field(config)).collect()
    }
}

// Implementation for BigInt<N>
impl<const M: usize, const N: usize> FieldMap<N> for BigInt<M> {
    type Output = RandomField<N>;

    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
        if config.is_null() {
            panic!("Cannot convert BigInt to prime field element without a modulus")
        }

        unsafe {
            let modulus: [u64; N] = (*config).modulus.0;

            let mut r: BigInt<N> = match M.cmp(&N) {
                std::cmp::Ordering::Less => {
                    let mut wider_value = [0u64; N];
                    wider_value[..N].copy_from_slice(&self.0);
                    BigInt(wider_value)
                }
                std::cmp::Ordering::Equal => {
                    let mut value = [0u64; N];
                    value.copy_from_slice(&self.0);
                    BigInt(value)
                }
                std::cmp::Ordering::Greater => {
                    let mut value = crypto_bigint::Uint::<M>::from_words(self.0);
                    let mut wider_modulus = [0u64; M];
                    wider_modulus[..M].copy_from_slice(&modulus);
                    let modu = crypto_bigint::Uint::<M>::from_words(wider_modulus);

                    value %= crypto_bigint::NonZero::new(modu).unwrap();
                    let mut result = [0u64; N];
                    result.copy_from_slice(&value.to_words()[..N]);

                    BigInt(result)
                }
            };

            // Apply Montgomery form transformation
            (*config).mul_assign(&mut r, &(*config).r2);
            RandomField::<N>::new_unchecked(config, r)
        }
    }
}

// Implementation for reference to BigInt<N>
impl<const M: usize, const N: usize> FieldMap<N> for &BigInt<M> {
    type Output = RandomField<N>;

    fn map_to_field(&self, config: *const FieldConfig<N>) -> Self::Output {
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
        let config: *const FieldConfig<1> = &create_field_config!(23);

        let bytes = [0b0000_0001]; // Value: 1 with leading zeros
        let expected = BigInt::<1>::from_bytes_le(&bytes)
            .unwrap()
            .map_to_field(config);

        let result = RandomField::<1>::from_bytes_le_with_config(config, &bytes);
        assert_eq!(result, Some(expected));
    }

    #[test]
    fn converts_from_bytes_be_with_config_leading_zeros() {
        let config = create_field_config!(23);
        let config_ptr = &config as *const FieldConfig<1>;

        let bytes = [0x01]; //1 with leading zeros (big-endian);

        let result = RandomField::<1>::from_bytes_be_with_config(config_ptr, &bytes).unwrap();
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
            let max = <$type>::from_str(&<$type>::MAX.to_string()).unwrap();
            let max_result = max.map_to_field($config);
            assert!(
                max_result.into_bigint() < BigInt::<$N>::from($field),
                "Maximum value should be less than field modulus"
            );

            // Test minimum value
            let min = -<$type>::from_str(&<$type>::MAX.to_string()).unwrap();
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
        let field_1 = 18446744069414584321_u64;
        let config_1: *const FieldConfig<1> = &create_field_config!(field_1);

        // Test primitive types with full range
        test_signed_type_full_range!(i8, field_1, config_1, 1);
        test_signed_type_full_range!(i16, field_1, config_1, 1);

        // Test larger primitive types with edge cases only
        test_signed_type_edge_cases!(i32, field_1, config_1, 1);
        test_signed_type_edge_cases!(i64, field_1, config_1, 1);
        test_signed_type_edge_cases!(i128, field_1, config_1, 1);

        // Test big integer types with edge cases
        test_signed_type_edge_cases!(I256, field_1, config_1, 1);
        test_signed_type_edge_cases!(I512, field_1, config_1, 1);
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
        let config_1: *const FieldConfig<1> = &create_field_config!(field_1);

        // Test small types with full range
        test_unsigned_type_full_range!(u8, field_1, config_1, 1);
        test_unsigned_type_full_range!(u16, field_1, config_1, 1);

        // Test larger types with edge cases only
        test_unsigned_type_edge_cases!(u32, field_1, config_1, 1);
        test_unsigned_type_edge_cases!(u64, field_1, config_1, 1);
        test_unsigned_type_edge_cases!(u128, field_1, config_1, 1);
    }

    #[test]
    #[should_panic(
        expected = "Cannot convert signed integer to prime field element without a modulus"
    )]
    fn test_signed_field_map_null_config() {
        let i32_val: i32 = 5;
        i32_val.map_to_field(std::ptr::null::<FieldConfig<1>>());
    }

    #[test]
    #[should_panic(
        expected = "Cannot convert unsigned integer to prime field element without a modulus"
    )]
    fn test_unsigned_field_map_null_config() {
        let u32_val: u32 = 5;
        u32_val.map_to_field(std::ptr::null::<FieldConfig<1>>());
    }
}

#[cfg(test)]
mod bigint_field_map_tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_bigint_smaller_than_field() {
        // Using a 2-limb field config with 1-limb BigInt
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = &config as *const _;

        let small_bigint = BigInt::<1>::from(12345u64);
        let result = small_bigint.map_to_field(config_ptr);

        assert_eq!(
            result.into_bigint().0[0],
            12345u64,
            "Small BigInt should be preserved in larger field"
        );
    }

    #[test]
    fn test_bigint_equal_size() {
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = &config as *const _;

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
        let config_ptr = &config as *const _;

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
        let config_ptr = &config as *const _;

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
        let config_ptr = &config as *const _;

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
        let _result: RandomField<2> = value.map_to_field(std::ptr::null());
    }

    #[test]
    fn test_bigint_max_value() {
        let modulus = BigInt::<2>::from_str("18446744069414584321").unwrap();
        let config = FieldConfig::new(modulus);
        let config_ptr = &config as *const _;

        // Create a BigInt with all bits set to 1
        let mut max_value = BigInt::<2>::zero();
        max_value.0.iter_mut().for_each(|x| *x = u64::MAX);

        let result = max_value.map_to_field(config_ptr);

        assert!(
            result.into_bigint() < modulus,
            "Result should be properly reduced modulo field modulus"
        );
    }
}
