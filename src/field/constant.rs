use ark_ff::{One, Zero};
use zeroize::Zeroize;

use crate::{
    biginteger::BigInt,
    field::{RandomField, RandomField::Raw},
};

impl<const N: usize> Zero for RandomField<'_, N> {
    fn zero() -> Self {
        Raw {
            value: BigInt::zero(),
        }
    }

    fn set_zero(&mut self) {
        *self.value_mut() = BigInt::zero()
    }

    fn is_zero(&self) -> bool {
        self.value().is_zero()
    }
}

impl<const N: usize> One for RandomField<'_, N> {
    fn one() -> Self {
        Raw {
            value: BigInt::one(),
        }
    }

    fn set_one(&mut self) {
        self.with_either_mut(
            |value| {
                *value = BigInt::one();
            },
            |config, value| {
                *value = *config.r();
            },
        );
    }

    fn is_one(&self) -> bool {
        self.with_either(
            |value| *value == BigInt::one(),
            |config, value| *value == *config.r(),
        )
    }
}

impl<const N: usize> Zeroize for RandomField<'_, N> {
    fn zeroize(&mut self) {
        unsafe { *self = ark_std::mem::zeroed() }
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::{One, Zero};
    use ark_std::str::FromStr;
    use zeroize::Zeroize;

    use crate::{
        biginteger::BigInt,
        create_bigint, create_random_field,
        field::RandomField,
        field_config::{ConfigRef, FieldConfig},
    };

    #[test]
    fn test_zero_creation() {
        let zero_elem = RandomField::<1>::zero();
        assert!(zero_elem.is_zero());
    }

    #[test]
    fn test_set_zero() {
        let config = FieldConfig::<1>::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);
        let mut elem = create_random_field!(config, 7);
        elem.set_zero();
        assert!(elem.is_zero());
    }

    #[test]
    fn test_one_creation() {
        let one_elem = RandomField::<1>::one();
        assert!(one_elem.is_one());
    }

    #[test]
    fn test_set_one() {
        let config = FieldConfig::<1>::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);
        let mut elem = create_random_field!(config, 5);
        elem.set_one();
        assert!(elem.is_one());
    }

    #[test]
    fn test_set_one_for_raw() {
        let mut raw_field = RandomField::<1>::zero();
        assert!(raw_field.is_zero());

        raw_field.set_one();

        assert!(raw_field.is_one());
    }

    #[test]
    fn test_is_zero_true() {
        let zero_elem = RandomField::<1>::zero();
        assert!(zero_elem.is_zero());
    }

    #[test]
    fn test_is_zero_false() {
        let non_zero_elem = RandomField::<1>::one();
        assert!(!non_zero_elem.is_zero());
    }

    #[test]
    fn test_is_one_true() {
        let one_elem = RandomField::<1>::one();
        assert!(one_elem.is_one());
    }

    #[test]
    fn test_is_one_false() {
        let non_one_elem = RandomField::<1>::from(3u32);
        assert!(!non_one_elem.is_one());
    }

    #[test]
    fn test_zeroize() {
        let config = FieldConfig::<1>::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);
        let mut elem = create_random_field!(config, 12);
        elem.zeroize();
        assert!(elem.is_zero());
    }

    #[test]
    fn test_zero_not_equal_one() {
        let zero_elem = RandomField::<1>::zero();
        let one_elem = RandomField::<1>::one();
        assert_ne!(zero_elem, one_elem);
    }
}
