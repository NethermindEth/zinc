use num_traits::{ConstOne, ConstZero, One, Zero};

use crate::field::{Form, Monty, RandomField};

impl<MOD: Monty<LIMBS>, const LIMBS: usize> Zero for RandomField<MOD, LIMBS> {
    fn zero() -> Self {
        Self(Form::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<MOD: Monty<LIMBS>, const LIMBS: usize> ConstZero for RandomField<MOD, LIMBS> {
    const ZERO: Self = Self(Form::ZERO);
}

impl<MOD: Monty<LIMBS>, const LIMBS: usize> One for RandomField<MOD, LIMBS> {
    fn one() -> Self {
        Self(Form::ONE)
    }
}

impl<MOD: Monty<LIMBS>, const LIMBS: usize> ConstOne for RandomField<MOD, LIMBS> {
    const ONE: Self = Self(Form::ONE);
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{U128, const_monty_params};
    use num_traits::{One, Zero};

    use super::*;

    // Define a small prime modulus for tests (copying style from Zinc tests where modulus ~ 128-bit)
    const_monty_params!(ModP, U128, "7fffffffffffffffffffffffffffffff");

    type F = RandomField<ModP, { U128::LIMBS }>;

    #[test]
    fn zero_one_basics() {
        let z = F::zero();
        assert!(z.is_zero());
        let o = F::one();
        assert!(!o.is_zero());
        assert_ne!(z, o);
    }

    #[test]
    fn from_bool_matches_one_zero() {
        let t: F = true.into();
        let f: F = false.into();
        assert_eq!(t, F::one());
        assert_eq!(f, F::zero());
    }
}
