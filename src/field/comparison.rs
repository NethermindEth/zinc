use crate::field::RandomField;
use crate::field::RandomField::{Initialized, Raw};
use ark_ff::{One, Zero};

impl<const N: usize> PartialEq for RandomField<N> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_one() & other.is_one() {
            return true;
        }
        if self.is_zero() && other.is_zero() {
            return true;
        }

        match (self, other) {
            (Initialized { .. }, Raw { .. }) | (Raw { .. }, Initialized { .. }) => false,
            (Raw { .. }, Raw { .. }) => self.value() == other.value(),
            (Initialized { .. }, Initialized { .. }) => {
                self.value() == other.value() && self.config_ptr() == other.config_ptr()
            }
        }
    }
}

impl<const N: usize> Eq for RandomField<N> {} // Eq requires PartialEq and ensures reflexivity.
