use ark_ff::{One, Zero};
use zeroize::Zeroize;

use crate::{
    field::{RandomField, RandomField::Raw},
    traits::{BigInteger, Config, ConfigReference},
};

impl<C: ConfigReference> Zero for RandomField<C> {
    fn zero() -> Self {
        Raw {
            value: C::B::zero(),
        }
    }

    fn set_zero(&mut self) {
        *self.value_mut() = C::B::zero()
    }

    fn is_zero(&self) -> bool {
        self.value().is_zero()
    }
}

impl<C: ConfigReference> One for RandomField<C> {
    fn one() -> Self {
        Raw { value: C::B::one() }
    }

    fn set_one(&mut self) {
        self.with_either_mut(
            |value| {
                *value = C::B::one();
            },
            |config, value| {
                *value = config.r().clone();
            },
        );
    }

    fn is_one(&self) -> bool {
        self.with_either(
            |value| *value == C::B::one(),
            |config, value| *value == *config.r(),
        )
    }
}

impl<C: ConfigReference> Zeroize for RandomField<C> {
    fn zeroize(&mut self) {
        unsafe { *self = ark_std::mem::zeroed() }
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::{One, Zero};
    use zeroize::Zeroize;

    use crate::{
        big_int,
        field::{RandomField, config::ConfigRef},
        field_config, random_field,
    };

    #[test]
    fn test_zero_creation() {
        let zero_elem = RandomField::<ConfigRef<1>>::zero();
        assert!(zero_elem.is_zero());
    }

    #[test]
    fn test_set_zero() {
        let config = field_config!(23, 1);
        let config = ConfigRef::<1>::from(&config);
        let mut elem = random_field!(7u32, config);
        elem.set_zero();
        assert!(elem.is_zero());
    }

    #[test]
    fn test_one_creation() {
        let one_elem = RandomField::<ConfigRef<1>>::one();
        assert!(one_elem.is_one());
    }

    #[test]
    fn test_set_one() {
        let config = field_config!(23, 1);
        let config = ConfigRef::<1>::from(&config);
        let mut elem = random_field!(5u32, config);
        elem.set_one();
        assert!(elem.is_one());
    }

    #[test]
    fn test_set_one_for_raw() {
        let mut raw_field = RandomField::<ConfigRef<1>>::zero();
        assert!(raw_field.is_zero());

        raw_field.set_one();

        assert!(raw_field.is_one());
    }

    #[test]
    fn test_is_zero_true() {
        let zero_elem = RandomField::<ConfigRef<1>>::zero();
        assert!(zero_elem.is_zero());
    }

    #[test]
    fn test_is_zero_false() {
        let non_zero_elem = RandomField::<ConfigRef<1>>::one();
        assert!(!non_zero_elem.is_zero());
    }

    #[test]
    fn test_is_one_true() {
        let one_elem = RandomField::<ConfigRef<1>>::one();
        assert!(one_elem.is_one());
    }

    #[test]
    fn test_is_one_false() {
        let non_one_elem = RandomField::<ConfigRef<1>>::from(3u32);
        assert!(!non_one_elem.is_one());
    }

    #[test]
    fn test_zeroize() {
        let config = field_config!(23, 1);
        let config = ConfigRef::<1>::from(&config);
        let mut elem: RandomField<_> = random_field!(12, config);
        elem.zeroize();
        assert!(elem.is_zero());
    }

    #[test]
    fn test_zero_not_equal_one() {
        let zero_elem = RandomField::<ConfigRef<1>>::zero();
        let one_elem = RandomField::one();
        assert_ne!(zero_elem, one_elem);
    }
}
