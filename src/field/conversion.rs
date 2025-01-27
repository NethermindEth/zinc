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
    use crate::{
        biginteger::{BigInt, BigInteger64},
        create_bigint, create_field_config,
        field::RandomField,
        field_config::FieldConfig,
    };
    use std::str::FromStr;

    #[test]
    fn test_from_u128() {
        let int = 243043087159742188419721163456177516u128;
        let raw_elem = RandomField::<2>::from(int);
        assert_eq!(
            raw_elem,
            RandomField::Raw {
                value: create_bigint!(2, 243043087159742188419721163456177516),
            }
        )
    }

    #[test]
    fn test_from_u32() {
        let int = 23u32;
        let raw_elem = RandomField::<1>::from(int);
        assert_eq!(
            raw_elem,
            RandomField::Raw {
                value: BigInteger64::from(23u32)
            }
        )
    }

    #[should_panic]
    #[test]
    fn test_failing_from_u128() {
        let int = 243043087159742188419721163456177516u128;
        let _ = RandomField::<1>::from(int);
    }

    #[test]
    fn test_bigint_conversion() {
        let config = create_field_config!(4, 695962179703626800597079116051991347);

        let bigint = create_bigint!(695962179703);

        let field_elem = RandomField::from_bigint(&config, bigint).unwrap();
        assert_eq!(bigint, field_elem.into_bigint());
        let bigint = create_bigint!(695962179703626800597079116051991346);

        let field_elem = RandomField::from_bigint(&config, bigint).unwrap();
        assert_eq!(bigint, field_elem.into_bigint())
    }
}
