use std::ops::{Add, Div, Mul, Neg, Sub};

use ark_ff::{BigInt, BigInteger, One, Zero};

use crate::field_config::{self, FieldConfig};

#[derive(Copy, Clone)]
pub struct RandomField<'config, const N: usize> {
    pub config: Option<&'config FieldConfig<N>>,
    pub value: BigInt<N>,
}

impl<'config, const N: usize> RandomField<'config, N> {
    fn new_unchecked(config: Option<&'config FieldConfig<N>>, value: BigInt<N>) -> Self {
        RandomField { config, value }
    }
    /// Convert from `BigInteger` to `RandomField`
    ///
    /// If `BigInteger` is greater then field modulus return `None`
    pub fn from_bigint(config: &'config FieldConfig<N>, value: BigInt<N>) -> Option<Self> {
        if value.is_zero() {
            Some(Self::new_unchecked(Some(config), value))
        } else if value >= config.modulus {
            None
        } else {
            let mut r = value;
            config.mul_assign(&mut r, &config.r2);
            Some(Self::new_unchecked(Some(config), r))
        }
    }

    pub fn into_bigint(&self) -> BigInt<N> {
        if self.is_zero() {
            return BigInt::zero();
        }

        if self.value == BigInt::one() && self.config.is_none() {
            return BigInt::one();
        }

        let config = self
            .config
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
}

impl<'config, const N: usize> Sub<RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn sub(self, rhs: RandomField<'config, N>) -> RandomField<'config, N> {
        &self - &rhs
    }
}

impl<'a, 'config, const N: usize> Sub<&'a RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn sub(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        &self - rhs
    }
}

impl<'a, 'config, const N: usize> Sub<&'a RandomField<'config, N>> for &RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn sub(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        if rhs.is_zero() {
            return *self;
        }
        if self.is_zero() {
            return -*rhs;
        }

        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.
        let lconfig = self
            .config
            .expect("This field element has no associated field");
        let rconfig = rhs
            .config
            .expect("This field element has no associated field");

        if lconfig != rconfig {
            panic!("cannot add field elements of different fields");
        }

        let mut res = RandomField::zero();
        res = res + *self;
        lconfig.sub_assign(&mut res.value, &rhs.value);
        res
    }
}

impl<'config, const N: usize> Add<RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: RandomField<'config, N>) -> RandomField<'config, N> {
        &self + &rhs
    }
}

impl<'a, 'config, const N: usize> Add<&'a RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        &self + rhs
    }
}

impl<'a, 'config, const N: usize> Add<&'a RandomField<'config, N>> for &RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        if rhs.is_zero() {
            return *self;
        }
        if self.is_zero() {
            return *rhs;
        }

        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.
        let lconfig = self
            .config
            .expect("This field element has no associated field");
        let rconfig = rhs
            .config
            .expect("This field element has no associated field");

        if lconfig != rconfig {
            panic!("cannot add field elements of different fields");
        }

        let mut res = RandomField::zero();
        res = res + *self;
        lconfig.add_assign(&mut res.value, &rhs.value);
        res
    }
}

impl<'config, const N: usize> Mul<RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn mul(self, rhs: RandomField<'config, N>) -> RandomField<'config, N> {
        &self * &rhs
    }
}

impl<'a, 'config, const N: usize> Mul<&'a RandomField<'config, N>> for &RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn mul(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        if self.is_one() {
            return *rhs;
        }
        if rhs.is_one() {
            return *self;
        }

        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.
        let lconfig = self
            .config
            .expect("This field element has no associated field");
        let rconfig = rhs
            .config
            .expect("This field element has no associated field");

        if lconfig != rconfig {
            panic!("cannot multiply field elements of different fields");
        }

        let mut res = RandomField::one();
        res = res * *self;
        lconfig.mul_assign(&mut res.value, &rhs.value);
        res
    }
}
impl<'config, const N: usize> Div<RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn div(self, rhs: RandomField<'config, N>) -> RandomField<'config, N> {
        &self / &rhs
    }
}

impl<'a, 'config, const N: usize> Div<&'a RandomField<'config, N>> for &RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn div(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        if rhs.is_zero() {
            panic!("Attempt to divide by zero");
        }

        if rhs.is_one() {
            return *self;
        }

        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.
        let lconfig = self
            .config
            .expect("This field element has no associated field");
        let rconfig = rhs
            .config
            .expect("This field element has no associated field");

        if lconfig != rconfig {
            panic!("cannot divide field elements of different fields");
        }

        let mut res = RandomField::one();
        res = res * *self;

        lconfig.mul_assign(&mut res.value, &lconfig.inverse(&rhs.value).unwrap());
        res
    }
}

impl<const N: usize> std::fmt::Debug for RandomField<'_, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.config {
            Some(config) => write!(f, "{} in the field Z_{}", self.value, config.modulus),
            None => write!(f, "{}", self.value),
        }
    }
}

impl<const N: usize> Zero for RandomField<'_, N> {
    fn zero() -> Self {
        Self::new_unchecked(None, BigInt::zero())
    }

    fn is_zero(&self) -> bool {
        self.value == BigInt::zero()
    }

    fn set_zero(&mut self) {
        self.value = BigInt::zero()
    }
}

impl<const N: usize> One for RandomField<'_, N> {
    fn one() -> Self {
        Self::new_unchecked(None, BigInt::one())
    }

    fn set_one(&mut self) {
        self.value = BigInt::one()
    }

    fn is_one(&self) -> bool {
        match self.config {
            Some(conf) => self.value == conf.r,
            None => self.value == BigInt::one(),
        }
    }
}

impl<const N: usize> PartialEq for RandomField<'_, N> {
    fn eq(&self, other: &Self) -> bool {
        if self.config.is_none() && other.config.is_none() {
            return self.is_one() && other.is_one() || self.is_zero() && other.is_zero();
        }

        let config_ptr_lhs: *const FieldConfig<N> = self.config.unwrap();
        let config_ptr_rhs: *const FieldConfig<N> = other.config.unwrap();

        self.value == other.value && config_ptr_lhs == config_ptr_rhs
    }
}

impl<const N: usize> Neg for RandomField<'_, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.is_zero() {
            return self;
        }
        let config = self
            .config
            .expect("This field element has no associated field");

        let mut val = config.modulus;
        val.sub_with_borrow(&self.value);
        Self::new_unchecked(self.config, val)
    }
}

impl<const N: usize> Eq for RandomField<'_, N> {}

unsafe impl<const N: usize> Send for RandomField<'_, N> {}

unsafe impl<const N: usize> Sync for RandomField<'_, N> {}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::{BigInteger256, BigInteger64, One};

    use crate::field_config::FieldConfig;

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
