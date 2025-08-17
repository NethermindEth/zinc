use ark_ff::{One, Zero};

use crate::{
    field::{
        RandomField,
        RandomField::{Initialized, Raw},
    },
    traits::ConfigReference,
};

impl<C: ConfigReference> PartialEq for RandomField<C> {
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

impl<C: ConfigReference> Eq for RandomField<C> {} // Eq requires PartialEq and ensures reflexivity.
