use crate::biginteger::BigInt;
use crate::field::RandomField;
use crate::field::RandomField::Raw;
use ark_ff::{One, Zero};
use zeroize::Zeroize;

impl<const N: usize> Zero for RandomField<N> {
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

impl<const N: usize> One for RandomField<N> {
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
                *value = config.r;
            },
        );
    }

    fn is_one(&self) -> bool {
        self.with_either(
            |value| *value == BigInt::one(),
            |config, value| *value == config.r,
        )
    }
}

impl<const N: usize> Zeroize for RandomField<N> {
    fn zeroize(&mut self) {
        unsafe { *self = std::mem::zeroed() }
    }
}
