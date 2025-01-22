#![allow(clippy::not_unsafe_ptr_arg_deref)]
use std::{
    iter::Sum,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign},
};

use ark_ff::{One, UniformRand, Zero};
use crypto_bigint::NonZero;
use zeroize::Zeroize;

use crate::{
    biginteger::BigInt,
    field_config::{self, FieldConfig},
};

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct RandomField<const N: usize> {
    pub config: *const FieldConfig<N>,
    pub value: BigInt<N>,
}

impl<const N: usize> UniformRand for RandomField<N> {
    fn rand<R: ark_std::rand::Rng + ?Sized>(rng: &mut R) -> Self {
        let value = BigInt::rand(rng);

        // Super unsafe. Once a number has been generated,
        // the config should be provided.
        Self {
            config: std::ptr::null(),
            value,
        }
    }
}

impl<const N: usize> RandomField<N> {
    /// Config setter that can be used after a `RandomField::rand(...)` call.
    pub fn set_config(mut self, config: *const FieldConfig<N>) -> Self {
        let modulus: BigInt<N> = unsafe { (*config).modulus };
        self.value = self.value % NonZero::new(modulus).unwrap();

        Self::from_bigint(config, self.value).expect("Should not end up with a None here.")
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
        if self.config.is_null() {
            return None;
        }

        unsafe { self.config.as_ref() }
    }
}

impl<const N: usize> RandomField<N> {
    fn new_unchecked(config: *const FieldConfig<N>, value: BigInt<N>) -> Self {
        RandomField { config, value }
    }
    /// Convert from `BigInteger` to `RandomField`
    ///
    /// If `BigInteger` is greater then field modulus return `None`
    pub fn from_bigint(config: *const FieldConfig<N>, value: BigInt<N>) -> Option<Self> {
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

    pub fn into_bigint(&self) -> BigInt<N> {
        if self.is_zero() {
            return BigInt::zero();
        }

        if self.value == BigInt::one() && self.config.is_null() {
            return BigInt::one();
        }

        let config = self
            .config_ref()
            .expect("This field element has no associated field");
        let mut r = self.value.0;
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

    fn increment_by_one(&mut self) {
        let mut value = std::mem::take(&mut self.value);
        let config = self.config_ref().expect("Cannot add one, field is None");
        config.add_assign(&mut value, &config.r);

        self.value = value;
    }

    fn has_no_config(&self) -> bool {
        self.config.is_null()
    }
}

impl<const N: usize> SubAssign<RandomField<N>> for RandomField<N> {
    fn sub_assign(&mut self, rhs: RandomField<N>) {
        self.sub_assign(&rhs);
    }
}

impl<'a, const N: usize> SubAssign<&'a RandomField<N>> for RandomField<N> {
    fn sub_assign(&mut self, rhs: &'a RandomField<N>) {
        if rhs.is_zero() {
            return;
        }

        if self.is_zero() {
            *self = -*rhs;
        }

        let mut value = std::mem::take(&mut self.value);
        let config = check_equal_configs(self, rhs);

        config.sub_assign(&mut value, &rhs.value);
        self.value = value;
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
        if rhs.is_zero() {
            return;
        }

        if self.is_zero() {
            *self = *rhs;
            return;
        }

        if self.is_one() && self.has_no_config() {
            *self = *rhs;
            self.increment_by_one();
            return;
        }

        if rhs.is_one() && rhs.has_no_config() {
            self.increment_by_one();
            return;
        }

        let mut value = std::mem::take(&mut self.value);
        let config = check_equal_configs(self, rhs);

        config.add_assign(&mut value, &rhs.value);
        self.value = value;
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
        if rhs.is_zero() {
            panic!("Attempt to divide by zero");
        }

        if rhs.is_one() && rhs.has_no_config() {
            return *self;
        }

        let config = check_equal_configs(self, rhs);
        let mut res = *self;
        config.mul_assign(&mut res.value, &config.inverse(&rhs.value).unwrap());
        res
    }
}

impl<'a, const N: usize> MulAssign<&'a Self> for RandomField<N> {
    fn mul_assign(&mut self, rhs: &'a Self) {
        if self.is_zero() || rhs.is_zero() {
            *self = RandomField::zero();
            return;
        }

        if self.is_one() {
            self.value = rhs.value;

            // If we do have a config we don't care.
            if self.has_no_config() && !rhs.has_no_config() {
                self.config = rhs.config;
            }

            return;
        }

        if rhs.is_one() {
            // If we do have a config we don't care.
            if self.has_no_config() && !rhs.has_no_config() {
                self.config = rhs.config;
            }

            return;
        }

        check_equal_configs(self, rhs);

        rhs.config_ref()
            .unwrap()
            .mul_assign(&mut self.value, &rhs.value);
    }
}

impl<const N: usize> MulAssign<Self> for RandomField<N> {
    fn mul_assign(&mut self, rhs: Self) {
        self.mul_assign(&rhs);
    }
}

impl<const N: usize> std::fmt::Debug for RandomField<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.config_ref() {
            Some(config) => write!(f, "{} in the field Z_{}", self.value, config.modulus),
            None => write!(f, "{}", self.value),
        }
    }
}

impl<const N: usize> std::fmt::Display for RandomField<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // TODO: we should go back from Montgomery here.
        write!(f, "{}", self.value)
    }
}

impl<const N: usize> Zero for RandomField<N> {
    fn zero() -> Self {
        Self::new_unchecked(
            std::ptr::null::<FieldConfig<N>>().cast_mut(),
            BigInt::zero(),
        )
    }

    fn is_zero(&self) -> bool {
        self.value == BigInt::zero()
    }

    fn set_zero(&mut self) {
        self.value = BigInt::zero()
    }
}

impl<const N: usize> One for RandomField<N> {
    fn one() -> Self {
        Self::new_unchecked(std::ptr::null::<FieldConfig<N>>().cast_mut(), BigInt::one())
    }

    fn set_one(&mut self) {
        self.value = BigInt::one()
    }

    fn is_one(&self) -> bool {
        match self.config_ref() {
            Some(conf) => self.value == conf.r,
            None => self.value == BigInt::one(),
        }
    }
}

impl<const N: usize> Neg for RandomField<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.is_zero() {
            return self;
        }
        let config = self
            .config_ref()
            .expect("This field element has no associated field");

        let mut val = config.modulus;
        val.sub_with_borrow(&self.value);
        Self::new_unchecked(self.config, val)
    }
}

unsafe impl<const N: usize> Send for RandomField<N> {}
unsafe impl<const N: usize> Sync for RandomField<N> {}

/// Checks if field configs are equal
/// Panics otherwise
pub fn check_equal_configs<'a, const N: usize>(
    l_element: &'a RandomField<N>,
    r_element: &'a RandomField<N>,
) -> &'a FieldConfig<N> {
    let lconfig = l_element
        .config_ref()
        .expect("This field element has no associated field");
    let rconfig = r_element
        .config_ref()
        .expect("This field element has no associated field");

    if lconfig != rconfig {
        panic!("Cannot operate on field elements of different fields");
    }

    lconfig
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::{One, Zero};

    use crate::{
        biginteger::{BigInteger256, BigInteger64},
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

    #[should_panic]
    #[test]
    fn test_add_two_ones() {
        let lhs: RandomField<1> = RandomField::one();

        let rhs = RandomField::one();

        let _ = lhs + rhs;
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
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

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
        assert_eq!(product, rhs);
    }

    #[test]
    fn test_right_mul_by_zero() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let rhs = BigInteger64::from_str("22").unwrap();

        let lhs = RandomField::zero();
        let rhs = RandomField::from_bigint(&field_config, rhs).unwrap();

        let product = lhs * rhs;
        assert_eq!(product, lhs);
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
        b = b / a;
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
    #[should_panic]
    fn test_failing_subtraction() {
        let field_config = FieldConfig::new(BigInteger64::from_str("23").unwrap());

        let lhs = BigInteger64::from_str("2").unwrap();

        let lhs = RandomField::from_bigint(&field_config, lhs).unwrap();
        let rhs = RandomField::one();
        let _ = lhs - rhs;
    }
}
