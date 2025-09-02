use std::cmp::Ordering;

use crypto_bigint::subtle::{Choice, ConstantTimeEq};

use crate::field::{Monty, RandomField};

impl<MOD: Monty<LIMBS>, const LIMBS: usize> ConstantTimeEq for RandomField<MOD, LIMBS> {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl<MOD: Monty<LIMBS>, const LIMBS: usize> PartialOrd for RandomField<MOD, LIMBS> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.as_montgomery().partial_cmp(other.0.as_montgomery())
    }
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{U128, const_monty_params};

    use super::*;

    const_monty_params!(ModP, U128, "7fffffffffffffffffffffffffffffff");
    type F = RandomField<ModP, { U128::LIMBS }>;

    #[test]
    fn const_time_eq_and_order() {
        let a: F = 10u64.into();
        let b: F = 10u64.into();
        let c: F = 11u64.into();
        assert_eq!(a.ct_eq(&b).unwrap_u8(), 1);
        assert_eq!(a.ct_eq(&c).unwrap_u8(), 0);
        assert!(a.partial_cmp(&c).is_some());
    }
}
