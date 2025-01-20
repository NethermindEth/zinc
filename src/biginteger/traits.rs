use std::ops::{Div, Rem, Shl, Shr};

use crypto_bigint::{
    subtle::{ConditionallySelectable, ConstantTimeEq, ConstantTimeGreater, ConstantTimeLess},
    CheckedAdd, CheckedMul, CheckedSub, Integer, Limb, NonZero, Zero,
};

use super::BigInt;

impl<const N: usize> Integer for BigInt<N> {
    const ONE: Self = BigInt::<N>::ONE;

    const MAX: Self = BigInt::<N>::MAX;

    const BITS: usize = BigInt::<N>::BITS;

    const BYTES: usize = BigInt::<N>::BYTES;

    const LIMBS: usize = N;

    fn is_odd(&self) -> crypto_bigint::subtle::Choice {
        todo!()
    }

    fn is_even(&self) -> crypto_bigint::subtle::Choice {
        !self.is_odd()
    }
}

impl<const N: usize> Zero for BigInt<N> {
    const ZERO: Self = BigInt::<N>::ZERO;

    fn is_zero(&self) -> crypto_bigint::subtle::Choice {
        self.ct_eq(&Self::ZERO)
    }
}

impl<const N: usize> Shr<usize> for BigInt<N> {
    type Output = Self;

    fn shr(self, rhs: usize) -> Self::Output {
        todo!()
    }
}

impl<const N: usize> Shl<usize> for BigInt<N> {
    type Output = Self;

    fn shl(self, rhs: usize) -> Self::Output {
        todo!()
    }
}
impl<const N: usize> ConstantTimeLess for BigInt<N> {
    fn ct_lt(&self, other: &Self) -> crypto_bigint::subtle::Choice {
        todo!()
    }
}

impl<const N: usize> ConstantTimeEq for BigInt<N> {
    fn ct_eq(&self, other: &Self) -> crypto_bigint::subtle::Choice {
        todo!()
    }
}

impl<const N: usize> ConstantTimeGreater for BigInt<N> {
    fn ct_gt(&self, other: &Self) -> crypto_bigint::subtle::Choice {
        todo!()
    }
}

impl<'a, const N: usize> CheckedAdd<&'a BigInt<N>> for BigInt<N> {
    type Output = Self;

    fn checked_add(&self, rhs: &'a BigInt<N>) -> crypto_bigint::subtle::CtOption<Self> {
        todo!()
    }
}

impl<'a, const N: usize> CheckedMul<&'a BigInt<N>> for BigInt<N> {
    type Output = Self;

    fn checked_mul(&self, rhs: &'a BigInt<N>) -> crypto_bigint::subtle::CtOption<Self> {
        todo!()
    }
}

impl<'a, const N: usize> CheckedSub<&'a BigInt<N>> for BigInt<N> {
    type Output = Self;

    fn checked_sub(&self, rhs: &'a BigInt<N>) -> crypto_bigint::subtle::CtOption<Self> {
        todo!()
    }
}

impl<const N: usize> ConditionallySelectable for BigInt<N> {
    fn conditional_select(a: &Self, b: &Self, choice: crypto_bigint::subtle::Choice) -> Self {
        todo!()
    }
}

impl<const N: usize> Rem<NonZero<BigInt<N>>> for BigInt<N> {
    type Output = Self;

    fn rem(self, rhs: NonZero<BigInt<N>>) -> Self::Output {
        todo!()
    }
}

impl<const N: usize> Div<NonZero<BigInt<N>>> for BigInt<N> {
    type Output = Self;

    fn div(self, rhs: NonZero<BigInt<N>>) -> Self::Output {
        todo!()
    }
}

impl<const N: usize> AsRef<[Limb]> for BigInt<N> {
    fn as_ref(&self) -> &[Limb] {
        todo!()
    }
}
