use crate::biginteger::BigInt;
use crate::field::RandomField;
use crate::field::RandomField::Raw;

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

impl<const N: usize> From<u64> for RandomField<N> {
    fn from(value: u64) -> Self {
        let value = BigInt::from(value);
        Raw { value }
    }
}

impl<const N: usize> From<u32> for RandomField<N> {
    fn from(value: u32) -> Self {
        let value = BigInt::from(value);
        Raw { value }
    }
}

impl<const N: usize> From<u16> for RandomField<N> {
    fn from(value: u16) -> Self {
        let value = BigInt::from(value);
        Raw { value }
    }
}
impl<const N: usize> From<u8> for RandomField<N> {
    fn from(value: u8) -> Self {
        let value = BigInt::from(value);
        Raw { value }
    }
}

impl<const N: usize> From<bool> for RandomField<N> {
    fn from(value: bool) -> Self {
        let value = BigInt::from(value as u8);
        Raw { value }
    }
}

#[cfg(test)]
mod tests {
    use crate::{biginteger::BigInt, field::RandomField};
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
}
