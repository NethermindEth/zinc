#![allow(clippy::not_unsafe_ptr_arg_deref)]

use std::{
    iter::Sum,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use ark_ff::{One, UniformRand, Zero};
use zeroize::Zeroize;

use crate::{
    biginteger::BigInt,
    field_config::{self, FieldConfig},
};

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub enum RandomField<const N: usize> {
    Raw {
        value: BigInt<N>,
    },
    Initialized {
        config: *const FieldConfig<N>,
        value: BigInt<N>,
    },
}

use RandomField::*;

impl<const N: usize> RandomField<N> {
    pub fn is_raw(&self) -> bool {
        matches!(self, Raw { .. })
    }

    pub fn is_initialized(&self) -> bool {
        matches!(self, Initialized { .. })
    }

    pub fn with_raw_value_or<F, A>(&self, f: F, default: A) -> A
    where
        F: Fn(&BigInt<N>) -> A,
    {
        match self {
            Raw { value } => f(value),
            _ => default,
        }
    }

    pub fn with_raw_value_mut_or<F, A>(&mut self, f: F, default: A) -> A
    where
        F: Fn(&mut BigInt<N>) -> A,
    {
        match self {
            Raw { value } => f(value),
            _ => default,
        }
    }

    pub fn with_init_value<'a, F, A>(&'a self, f: F) -> Option<A>
    where
        F: Fn(&'a FieldConfig<N>, &'a BigInt<N>) -> A,
    {
        match self {
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");
                Some(f(config, value))
            },
            _ => None,
        }
    }

    pub fn with_init_value_or<'a, F, A>(&'a self, f: F, default: A) -> A
    where
        F: Fn(&'a FieldConfig<N>, &'a BigInt<N>) -> A,
    {
        match self {
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");
                f(config, value)
            },
            _ => default,
        }
    }

    pub fn with_either<'a, R, I, A>(&'a self, raw_fn: R, init_fn: I) -> A
    where
        I: Fn(&'a FieldConfig<N>, &'a BigInt<N>) -> A,
        R: Fn(&'a BigInt<N>) -> A,
    {
        match self {
            Raw { value } => raw_fn(value),
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                init_fn(config, value)
            },
        }
    }

    pub fn with_either_mut<'a, R, I, A>(&'a mut self, raw_fn: R, init_fn: I) -> A
    where
        I: Fn(&'a FieldConfig<N>, &'a mut BigInt<N>) -> A,
        R: Fn(&'a mut BigInt<N>) -> A,
    {
        match self {
            Raw { value } => raw_fn(value),
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                init_fn(config, value)
            },
        }
    }

    pub fn with_either_owned<R, I, A>(self, raw_fn: R, init_fn: I) -> A
    where
        I: Fn(&FieldConfig<N>, BigInt<N>) -> A,
        R: Fn(BigInt<N>) -> A,
    {
        match self {
            Raw { value } => raw_fn(value),
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                init_fn(config, value)
            },
        }
    }

    #[inline(always)]
    pub fn value(&self) -> &BigInt<N> {
        match self {
            Raw { value } => value,
            Initialized { value, .. } => value,
        }
    }

    #[inline(always)]
    pub fn value_mut(&mut self) -> &mut BigInt<N> {
        match self {
            Raw { value } => value,
            Initialized { value, .. } => value,
        }
    }

    pub fn with_aligned_config_mut<F, G, A>(
        &mut self,
        rhs: &Self,
        with_config: F,
        without_config: G,
    ) -> A
    where
        F: Fn(&mut BigInt<N>, &BigInt<N>, &FieldConfig<N>) -> A,
        G: Fn(&mut BigInt<N>, &BigInt<N>) -> A,
    {
        match (self, rhs) {
            (Raw { value: value_self }, Raw { value: rhs }) => without_config(value_self, rhs),
            (
                Initialized {
                    value: value_self,
                    config,
                },
                Initialized {
                    value: value_rhs, ..
                },
            ) => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                with_config(value_self, value_rhs, config)
            },
            (
                Initialized {
                    value: value_self,
                    config,
                },
                rhs @ Raw { .. },
            ) => unsafe {
                let rhs = (*rhs).set_config_owned(*config);
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                with_config(value_self, rhs.value(), config)
            },
            (
                lhs @ Raw { .. },
                Initialized {
                    value: value_rhs,
                    config,
                },
            ) => unsafe {
                lhs.set_config(*config);
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                with_config(lhs.value_mut(), value_rhs, config)
            },
        }
    }
}

impl<const N: usize> UniformRand for RandomField<N> {
    fn rand<R: ark_std::rand::Rng + ?Sized>(rng: &mut R) -> Self {
        let value = BigInt::rand(rng);

        Self::Raw { value }
    }
}

impl<const N: usize> RandomField<N> {
    pub fn set_config(&mut self, config: *const FieldConfig<N>) {
        self.with_raw_value_mut_or(
            |value| {
                // Ideally we should do something like:
                //
                // ```
                // let modulus: BigInt<N> = unsafe { (*config).modulus };
                // *value = *value % modulus;
                // ```
                //
                // but we don't have `mod` out of the box.
                // So let's hope we don't exceed the modulus.

                // TODO: prettify this
                *value = *Self::from_bigint(config, *value)
                    .expect("Should not end up with a None here.")
                    .value();
            },
            (),
        );

        let value = std::mem::take(self.value_mut());

        *self = Initialized { config, value }
    }

    /// Config setter that can be used after a `RandomField::rand(...)` call.
    pub fn set_config_owned(mut self, config: *const FieldConfig<N>) -> Self {
        self.set_config(config);
        self
    }
}

// TODO: Finalise this
//impl<const N: usize> AdditiveGroup for RandomField<N> {
//    type Scalar = ;
//
//    const ZERO: Self = Self { config: std::ptr::null(), BigInt::zero() };
//}

impl<const N: usize> Zeroize for RandomField<N> {
    fn zeroize(&mut self) {
        unsafe { *self = std::mem::zeroed() }
    }
}

impl<const N: usize> RandomField<N> {
    #[inline(always)]
    pub fn config_ref(&self) -> Option<&FieldConfig<N>> {
        self.with_init_value(|config, _| config)
    }
}

impl<const N: usize> RandomField<N> {
    fn new_unchecked(config: *const FieldConfig<N>, value: BigInt<N>) -> Self {
        Initialized { config, value }
    }

    /// Convert from `BigInteger` to `RandomField`
    ///
    /// If `BigInteger` is greater then field modulus return `None`
    pub fn from_bigint(config: *const FieldConfig<N>, value: BigInt<N>) -> Option<Self> {
        if config.is_null() {
            return Some(Raw { value });
        }

        unsafe {
            if value.is_zero() {
                Some(Self::zero())
            } else if value >= (*config).modulus {
                None
            } else {
                let mut r = value;
                (*config).mul_assign(&mut r, &(*config).r2);
                Some(Self::new_unchecked(config, r))
            }
        }
    }

    pub fn into_bigint(self) -> BigInt<N> {
        self.with_either_owned(|value| value, Self::demontgomery)
    }

    fn demontgomery(config: &FieldConfig<N>, value: BigInt<N>) -> BigInt<N> {
        let mut r = value.0;
        // Montgomery Reduction
        for i in 0..N {
            let k = r[i].wrapping_mul(config.inv);
            let mut carry = 0;

            field_config::mac_with_carry(r[i], k, config.modulus.0[0], &mut carry);
            for j in 1..N {
                r[(j + i) % N] = field_config::mac_with_carry(
                    r[(j + i) % N],
                    k,
                    config.modulus.0[j],
                    &mut carry,
                );
            }
            r[i % N] = carry;
        }

        BigInt::new(r)
    }
}

impl<const N: usize> SubAssign<RandomField<N>> for RandomField<N> {
    fn sub_assign(&mut self, rhs: RandomField<N>) {
        self.sub_assign(&rhs);
    }
}

impl<'a, const N: usize> SubAssign<&'a RandomField<N>> for RandomField<N> {
    fn sub_assign(&mut self, rhs: &'a RandomField<N>) {
        self.with_aligned_config_mut(
            rhs,
            |lhs, rhs, config| {
                config.sub_assign(lhs, rhs);
            },
            |lhs, rhs| {
                lhs.sub_with_borrow(rhs);
            },
        );
    }
}

impl<const N: usize> Sub<RandomField<N>> for RandomField<N> {
    type Output = RandomField<N>;

    fn sub(self, rhs: RandomField<N>) -> RandomField<N> {
        &self - &rhs
    }
}

impl<'a, const N: usize> Sub<&'a RandomField<N>> for &RandomField<N> {
    type Output = RandomField<N>;

    fn sub(self, rhs: &'a RandomField<N>) -> RandomField<N> {
        let mut res = *self;
        res.sub_assign(rhs);

        res
    }
}

impl<'a, const N: usize> AddAssign<&'a RandomField<N>> for RandomField<N> {
    fn add_assign(&mut self, rhs: &'a RandomField<N>) {
        self.with_aligned_config_mut(
            rhs,
            |lhs, rhs, config| {
                config.add_assign(lhs, rhs);
            },
            |lhs, rhs| {
                lhs.add_with_carry(rhs);
            },
        );
    }
}

impl<const N: usize> Add<RandomField<N>> for RandomField<N> {
    type Output = RandomField<N>;

    fn add(self, rhs: RandomField<N>) -> RandomField<N> {
        &self + &rhs
    }
}

impl<'a, const N: usize> Add<&'a RandomField<N>> for &RandomField<N> {
    type Output = RandomField<N>;

    fn add(self, rhs: &'a RandomField<N>) -> RandomField<N> {
        let mut res = *self;

        res.add_assign(rhs);

        res
    }
}

impl<'a, const N: usize> Sum<&'a RandomField<N>> for RandomField<N> {
    fn sum<I: Iterator<Item = &'a RandomField<N>>>(iter: I) -> Self {
        iter.fold(Self::zero(), |mut acc, x| {
            acc.add_assign(x);
            acc
        })
    }
}

impl<const N: usize> Mul<RandomField<N>> for RandomField<N> {
    type Output = RandomField<N>;

    fn mul(self, rhs: RandomField<N>) -> RandomField<N> {
        &self * &rhs
    }
}

impl<'a, const N: usize> Mul<&'a RandomField<N>> for &RandomField<N> {
    type Output = RandomField<N>;

    fn mul(self, rhs: &'a RandomField<N>) -> RandomField<N> {
        let mut res = *self;
        res.mul_assign(rhs);

        res
    }
}

impl<'a, const N: usize> Mul<&'a RandomField<N>> for RandomField<N> {
    type Output = RandomField<N>;

    fn mul(self, rhs: &'a RandomField<N>) -> RandomField<N> {
        let mut res = self;
        res.mul_assign(rhs);

        res
    }
}

impl<const N: usize> Div<RandomField<N>> for RandomField<N> {
    type Output = RandomField<N>;

    fn div(self, rhs: RandomField<N>) -> RandomField<N> {
        &self / &rhs
    }
}

impl<'a, const N: usize> Div<&'a RandomField<N>> for &RandomField<N> {
    type Output = RandomField<N>;
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: &'a RandomField<N>) -> RandomField<N> {
        let mut res = *self;
        res /= rhs;

        res
    }
}

impl<'a, const N: usize> MulAssign<&'a Self> for RandomField<N> {
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.with_aligned_config_mut(
            rhs,
            |lhs, rhs, config| {
                config.mul_assign(lhs, rhs);
            },
            |lhs, rhs| {
                lhs.mul(rhs);
            },
        );
    }
}

impl<const N: usize> MulAssign<Self> for RandomField<N> {
    fn mul_assign(&mut self, rhs: Self) {
        self.mul_assign(&rhs);
    }
}

impl<const N: usize> std::fmt::Debug for RandomField<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Raw { value } => write!(f, "{}, no config", value),
            self_ => write!(
                f,
                "{} in Z_{}",
                self_.value(),
                self.config_ref().unwrap().modulus
            ),
        }
    }
}

impl<const N: usize> std::fmt::Display for RandomField<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // TODO: we should go back from Montgomery here.
        write!(f, "{}", self.value())
    }
}

impl<const N: usize> Zero for RandomField<N> {
    fn zero() -> Self {
        Raw {
            value: BigInt::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.value().is_zero()
    }

    fn set_zero(&mut self) {
        *self.value_mut() = BigInt::zero()
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

impl<const N: usize> Neg for RandomField<N> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        if self.is_zero() {
            return self;
        }

        self.with_either_mut(
            |_| panic!("Cannot negate without a field config"),
            |config, value| {
                let tmp = *value;
                *value = config.modulus;
                value.sub_with_borrow(&tmp);
            },
        );

        self
    }
}

impl<const N: usize> DivAssign<Self> for RandomField<N> {
    fn div_assign(&mut self, rhs: Self) {
        self.div_assign(&rhs);
    }
}

impl<'a, const N: usize> DivAssign<&'a Self> for RandomField<N> {
    fn div_assign(&mut self, rhs: &'a Self) {
        if rhs.is_zero() {
            panic!("Attempt to divide by zero");
        }

        self.with_aligned_config_mut(
            rhs,
            |lhs, rhs, config| {
                config.mul_assign(lhs, &config.inverse(rhs).unwrap());
            },
            |_, _| panic!("Cannot divide without a field config"),
        );
    }
}

impl<'a, const N: usize> DivAssign<&'a mut Self> for RandomField<N> {
    fn div_assign(&mut self, rhs: &'a mut Self) {
        *self /= *rhs;
    }
}
impl<'a, const N: usize> Div<&'a Self> for RandomField<N> {
    type Output = Self;

    fn div(self, rhs: &'a Self) -> Self::Output {
        self / *rhs
    }
}
impl<'a, const N: usize> Div<&'a mut Self> for RandomField<N> {
    type Output = Self;

    fn div(self, rhs: &'a mut Self) -> Self::Output {
        self / *rhs
    }
}

impl<'a, const N: usize> core::iter::Product<&'a Self> for RandomField<N> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), core::ops::Mul::mul)
    }
}

unsafe impl<const N: usize> Send for RandomField<N> {}
unsafe impl<const N: usize> Sync for RandomField<N> {}

impl<const N: usize> From<u128> for RandomField<N> {
    fn from(other: u128) -> Self {
        let mut value = BigInt::default();
        if N == 1 {
            panic!("Integer is 128 bits but field is 64 bits")
        } else {
            value.0[0] = ((other << 64) >> 64) as u64;
            value.0[1] = (other >> 64) as u64;
        }
        Raw { value }
    }
}

impl<const N: usize> From<u64> for RandomField<N> {
    fn from(value: u64) -> Self {
        let value = BigInt::from(value);
        Raw { value }
    }
}

impl<const N: usize> From<u32> for RandomField<N> {
    fn from(value: u32) -> Self {
        let value = BigInt::from(value);
        Raw { value }
    }
}

impl<const N: usize> From<u16> for RandomField<N> {
    fn from(value: u16) -> Self {
        let value = BigInt::from(value);
        Raw { value }
    }
}
impl<const N: usize> From<u8> for RandomField<N> {
    fn from(value: u8) -> Self {
        let value = BigInt::from(value);
        Raw { value }
    }
}

impl<const N: usize> From<bool> for RandomField<N> {
    fn from(value: bool) -> Self {
        let value = BigInt::from(value as u8);
        Raw { value }
    }
}

impl<const N: usize> From<u128> for RandomField<N> {
    fn from(_value: u128) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::{One, Zero};

    use crate::{
        biginteger::{BigInt, BigInteger128, BigInteger256, BigInteger64},
        field_config::FieldConfig,
    };

    use super::RandomField;

    #[test]
    fn test_bigint_conversion() {
        let field_config = FieldConfig::new(
            BigInteger256::from_str("695962179703626800597079116051991347").unwrap(),
        );

        let bigint = BigInteger256::from_str("695962179703").unwrap();

        let field_elem = RandomField::from_bigint(&field_config, bigint).unwrap();
        assert_eq!(bigint, field_elem.into_bigint());
        let bigint = BigInteger256::from_str("695962179703626800597079116051991346").unwrap();

        let field_elem = RandomField::from_bigint(&field_config, bigint).unwrap();
        assert_eq!(bigint, field_elem.into_bigint())
    }

    #[test]
    fn test_addition() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("22").unwrap();
        let rhs = BigInteger64::from_str("2").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let sum = lhs + rhs;
        assert_eq!(sum.into_bigint(), BigInteger64::one());

        // Test 2
        let lhs = BigInteger64::from_str("20").unwrap();
        let rhs = BigInteger64::from_str("20").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let sum = lhs + rhs;
        assert_eq!(sum.into_bigint(), BigInteger64::from_str("17").unwrap())
    }

    #[test]
    fn test_add_one() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("22").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::one();

        let sum = lhs + rhs;
        assert_eq!(sum.into_bigint(), BigInteger64::zero());

        let sum = rhs + lhs;
        assert_eq!(sum.into_bigint(), BigInteger64::zero());
    }

    #[test]
    fn test_add_two_ones() {
        let lhs: RandomField<1> = RandomField::one();

        let rhs = RandomField::one();

        assert_eq!(
            lhs + rhs,
            RandomField::Raw {
                value: BigInt::from(2u32)
            }
        );
    }

    #[test]
    fn test_multiplication() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("22").unwrap();
        let rhs = BigInteger64::from_str("2").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let product = lhs * rhs;
        assert_eq!(product.into_bigint(), BigInteger64::from_str("21").unwrap());

        // Test 2
        let lhs = BigInteger64::from_str("20").unwrap();
        let rhs = BigInteger64::from_str("20").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs: RandomField<1> = RandomField::from_bigint(&field_config, rhs).unwrap();

        let product = lhs * rhs;
        assert_eq!(product.into_bigint(), BigInteger64::from_str("9").unwrap())
    }

    #[test]
    fn test_left_mul_by_zero() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("22").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::zero();

        let product = lhs * rhs;
        assert!(product.is_zero());
    }

    #[test]
    fn test_right_mul_by_zero() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let rhs = BigInteger64::from_str("22").unwrap();

        let lhs = RandomField::zero();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let product = lhs * rhs;
        assert!(product.is_zero());
    }

    #[test]
    fn test_division() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("22").unwrap();
        let rhs = BigInteger64::from_str("2").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let quotient = lhs / rhs;
        assert_eq!(
            quotient.into_bigint(),
            BigInteger64::from_str("11").unwrap()
        );

        // Test 2
        let lhs = BigInteger64::from_str("20").unwrap();
        let rhs = BigInteger64::from_str("20").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let quotient = lhs / rhs;
        assert_eq!(quotient.into_bigint(), BigInteger64::from_str("1").unwrap());

        // Test 3
        let lhs = BigInteger64::from_str("17").unwrap();
        let rhs = BigInteger64::from_str("4").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let quotient = lhs / rhs;
        assert_eq!(
            quotient.into_bigint(),
            BigInteger64::from_str("10").unwrap()
        )
    }

    #[test]
    #[should_panic]
    fn test_division_by_zero() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("17").unwrap();
        let rhs = BigInteger64::zero();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let _sum = lhs / rhs;
    }

    #[test]
    fn test_big_division() {
        let config = FieldConfig::new(
            BigInteger256::from_str("695962179703626800597079116051991347").unwrap(),
        );

        let a = RandomField::from_bigint(&config, BigInteger256::from_str("3").unwrap()).unwrap();
        let mut b = RandomField::from_bigint(&config, BigInteger256::one()).unwrap();
        b /= a;
        assert_eq!(
            b.into_bigint(),
            BigInteger256::from_str("231987393234542266865693038683997116").unwrap()
        );

        let a =
            RandomField::from_bigint(&config, BigInteger256::from_str("19382769832175").unwrap())
                .unwrap();

        let b =
            RandomField::from_bigint(&config, BigInteger256::from_str("97133987132135").unwrap())
                .unwrap();
        assert_eq!(
            BigInteger256::from_str("243043087159742188419721163456177516").unwrap(),
            (b / a).into_bigint()
        );
    }
    #[test]
    fn test_negation() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let op_val = BigInteger64::from_str("22").unwrap();

        let operand = RandomField::from_bigint(&field_config, op_val).unwrap();

        let negated = -operand;
        assert_eq!(negated.into_bigint(), BigInteger64::from_str("1").unwrap());

        // Test 2
        let op_val = BigInteger64::from_str("17").unwrap();

        let operand = RandomField::from_bigint(&field_config, op_val).unwrap();

        let negated = -operand;
        assert_eq!(negated.into_bigint(), BigInteger64::from_str("6").unwrap());

        // test with zero
        let op_val = BigInteger64::from_str("0").unwrap();

        let operand = RandomField::from_bigint(&field_config, op_val).unwrap();

        let negated = -operand;
        assert_eq!(negated.into_bigint(), BigInteger64::from_str("0").unwrap());
    }

    #[test]
    fn test_subtraction() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("2").unwrap();
        let rhs = BigInteger64::from_str("22").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let sum = lhs - rhs;
        assert_eq!(sum.into_bigint(), BigInteger64::from_str("3").unwrap());

        // Test 2
        let lhs = BigInteger64::from_str("20").unwrap();
        let rhs = BigInteger64::from_str("20").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let sum = lhs - rhs;
        assert_eq!(sum.into_bigint(), BigInteger64::zero())
    }

    #[test]
    fn test_init_sub_raw() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("2").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::one();
        let res = lhs - rhs;
        let mut expected = lhs;
        expected.set_one();
        assert_eq!(res, expected)
    }

    #[test]
    fn test_from_u128() {
        let int = 243043087159742188419721163456177516u128;
        let raw_elem = RandomField::<2>::from(int);
        assert_eq!(
            raw_elem,
            RandomField::Raw {
                value: BigInteger128::from_str("243043087159742188419721163456177516").unwrap()
            }
        )
    }

    #[test]
    fn test_from_u32() {
        let int = 23u32;
        let raw_elem = RandomField::<1>::from(int);
        assert_eq!(
            raw_elem,
            RandomField::Raw {
                value: BigInteger64::from(23u32)
            }
        )
    }

    #[should_panic]
    #[test]
    fn test_failing_from_u128() {
        let int = 243043087159742188419721163456177516u128;
        let _ = RandomField::<1>::from(int);
    }
}
