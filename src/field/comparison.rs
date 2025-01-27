use crate::field::RandomField;
use ark_ff::{One, Zero};

impl<const N: usize> PartialEq for RandomField<N> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_one() & other.is_one() {
            return true;
        }
        if self.is_zero() && other.is_zero() {
            return true;
        }

        if self.is_raw() && other.is_initialized() || self.is_initialized() && other.is_raw() {
            return false;
        }
        if self.is_raw() {
            self.value() == other.value()
        } else {
            self.value() == other.value()
                && self.config_ref().unwrap().modulus == other.config_ref().unwrap().modulus
        }
    }
}

impl<const N: usize> Eq for RandomField<N> {} // Eq requires PartialEq and ensures reflexivity.
