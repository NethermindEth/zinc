use std::ops::{Div, Rem, Shl, Shr};

use crypto_bigint::{
    subtle::{
        Choice, ConditionallySelectable, ConstantTimeEq, ConstantTimeGreater, ConstantTimeLess,
        CtOption,
    },
    CheckedAdd, CheckedMul, CheckedSub, Integer, Limb, NonZero, Zero,
};

use super::{BigInt, BigInteger};

impl<const N: usize> Integer for BigInt<N> {
    const ONE: Self = BigInt::<N>::ONE;

    const MAX: Self = BigInt::<N>::MAX;

    const BITS: usize = BigInt::<N>::BITS;

    const BYTES: usize = BigInt::<N>::BYTES;

    const LIMBS: usize = N;

    fn is_odd(&self) -> Choice {
        Choice::from(self.const_is_odd() as u8)
    }

    fn is_even(&self) -> Choice {
        !Integer::is_odd(self)
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

    fn shr(mut self, rhs: usize) -> Self::Output {
        self >>= rhs;
        self
    }
}

impl<const N: usize> Shl<usize> for BigInt<N> {
    type Output = Self;

    fn shl(mut self, rhs: usize) -> Self::Output {
        self <<= rhs;
        self
    }
}

impl<const N: usize> ConstantTimeLess for BigInt<N> {
    fn ct_lt(&self, other: &Self) -> crypto_bigint::subtle::Choice {
        ((self < other) as u8).into()
    }
}

impl<const N: usize> ConstantTimeEq for BigInt<N> {
    fn ct_eq(&self, other: &Self) -> crypto_bigint::subtle::Choice {
        ((self == other) as u8).into()
    }
}

impl<const N: usize> ConstantTimeGreater for BigInt<N> {
    fn ct_gt(&self, other: &Self) -> crypto_bigint::subtle::Choice {
        ((self > other) as u8).into()
    }
}

impl<'a, const N: usize> CheckedAdd<&'a BigInt<N>> for BigInt<N> {
    type Output = Self;

    fn checked_add(&self, rhs: &'a BigInt<N>) -> crypto_bigint::subtle::CtOption<Self> {
        let mut res = *self;
        let carry = res.add_with_carry(rhs);
        CtOption::new(res, (carry as u8).into())
    }
}

impl<'a, const N: usize> CheckedMul<&'a BigInt<N>> for BigInt<N> {
    type Output = Self;

    fn checked_mul(&self, rhs: &'a BigInt<N>) -> crypto_bigint::subtle::CtOption<Self> {
        let (lo, hi) = self.mul(rhs);
        CtOption::new(lo, crypto_bigint::Zero::is_zero(&hi))
    }
}

impl<'a, const N: usize> CheckedSub<&'a BigInt<N>> for BigInt<N> {
    type Output = Self;

    fn checked_sub(&self, rhs: &'a BigInt<N>) -> crypto_bigint::subtle::CtOption<Self> {
        let mut res = *self;
        let carry = res.sub_with_borrow(rhs);
        CtOption::new(res, (carry as u8).into())
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
