use crate::field::RandomField;
use ark_ff::{One, Zero};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

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

impl<const N: usize> Mul<RandomField<N>> for RandomField<N> {
    type Output = RandomField<N>;

    fn mul(self, rhs: RandomField<N>) -> RandomField<N> {
        &self * &rhs
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

impl<'a, const N: usize> Mul<&'a RandomField<N>> for &RandomField<N> {
    type Output = RandomField<N>;

    fn mul(self, rhs: &'a RandomField<N>) -> RandomField<N> {
        let mut res = *self;
        res.mul_assign(rhs);

        res
    }
}

impl<const N: usize> MulAssign<Self> for RandomField<N> {
    fn mul_assign(&mut self, rhs: Self) {
        self.mul_assign(&rhs);
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

impl<'a, const N: usize> Sum<&'a RandomField<N>> for RandomField<N> {
    fn sum<I: Iterator<Item = &'a RandomField<N>>>(iter: I) -> Self {
        iter.fold(Self::zero(), |mut acc, x| {
            acc.add_assign(x);
            acc
        })
    }
}

impl<const N: usize> Sum<RandomField<N>> for RandomField<N> {
    fn sum<I: Iterator<Item = RandomField<N>>>(iter: I) -> Self {
        iter.fold(Self::zero(), |mut acc, x| {
            acc.add_assign(&x);
            acc
        })
    }
}

impl<'a, const N: usize> core::iter::Product<&'a Self> for RandomField<N> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), core::ops::Mul::mul)
    }
}

#[cfg(test)]
mod test {
    use crate::biginteger::{BigInt, BigInteger256, BigInteger64};
    use crate::field::RandomField;
    use crate::field_config::FieldConfig;
    use ark_ff::{One, Zero};
    use std::str::FromStr;

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
}
