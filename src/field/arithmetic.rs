use ark_ff::{One, Zero};
use ark_std::{
    iter::Sum,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{field::RandomField, traits::Config};

macro_rules! impl_ops {
    (
        impl($($gen:tt)*) for $type:ty,
        $trait:ident, $op:ident,
        $trait_assign:ident, $op_assign:ident
    ) => {
        impl<$($gen)*> $trait<Self> for $type {
            type Output = Self;

            fn $op(mut self, rhs: Self) -> Self::Output {
                self.$op_assign(&rhs);
                self
            }
        }

        impl<$($gen)*> $trait<&Self> for $type {
            type Output = Self;

            fn $op(mut self, rhs: &Self) -> Self::Output {
                self.$op_assign(rhs);
                self
            }
        }

        impl<$($gen)*> $trait<Self> for &$type {
            type Output = $type;

            fn $op(self, rhs: Self) -> Self::Output {
                let mut res = *self;
                res.$op_assign(rhs);
                res
            }
        }

        impl<$($gen)*> $trait<$type> for &$type {
            type Output = $type;

            fn $op(self, rhs: $type) -> Self::Output {
                let mut res = *self;
                res.$op_assign(&rhs);
                res
            }
        }

        impl<$($gen)*> $trait_assign<Self> for $type {
            fn $op_assign(&mut self, rhs: Self) {
                self.$op_assign(&rhs);
            }
        }
    };
}

impl_ops!(impl('cfg, const N: usize) for RandomField<'cfg, N>, Add, add, AddAssign, add_assign);
impl_ops!(impl('cfg, const N: usize) for RandomField<'cfg, N>, Sub, sub, SubAssign, sub_assign);
impl_ops!(impl('cfg, const N: usize) for RandomField<'cfg, N>, Mul, mul, MulAssign, mul_assign);
impl_ops!(impl('cfg, const N: usize) for RandomField<'cfg, N>, Div, div, DivAssign, div_assign);

impl<const N: usize> AddAssign<&Self> for RandomField<'_, N> {
    fn add_assign(&mut self, rhs: &Self) {
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

impl<const N: usize> SubAssign<&Self> for RandomField<'_, N> {
    fn sub_assign(&mut self, rhs: &Self) {
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

impl<const N: usize> MulAssign<&Self> for RandomField<'_, N> {
    fn mul_assign(&mut self, rhs: &Self) {
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

impl<const N: usize> DivAssign<&Self> for RandomField<'_, N> {
    fn div_assign(&mut self, rhs: &Self) {
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

impl<const N: usize> DivAssign<&mut Self> for RandomField<'_, N> {
    fn div_assign(&mut self, rhs: &mut Self) {
        *self /= *rhs;
    }
}

impl<const N: usize> Neg for RandomField<'_, N> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        if self.is_zero() {
            return self;
        }

        self.with_either_mut(
            |_| panic!("Cannot negate without a field config"),
            |config, value| {
                let tmp = *value;
                *value = *config.modulus();
                value.sub_with_borrow(&tmp);
            },
        );

        self
    }
}

impl<'a, const N: usize> Sum<&'a Self> for RandomField<'_, N> {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |mut acc, x| {
            acc.add_assign(x);
            acc
        })
    }
}

impl<const N: usize> Sum for RandomField<'_, N> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |mut acc, x| {
            acc.add_assign(&x);
            acc
        })
    }
}

impl<'a, const N: usize> core::iter::Product<&'a Self> for RandomField<'_, N> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), core::ops::Mul::mul)
    }
}

#[cfg(test)]
mod test {
    use ark_ff::{One, Zero};
    use ark_std::str::FromStr;

    use crate::{
        biginteger::BigInt,
        create_bigint, create_random_field,
        field::RandomField,
        field_config::{ConfigRef, FieldConfig},
        traits::Config,
    };

    #[test]
    fn test_add_wrapping_around_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 22);
        let rhs: RandomField<1> = create_random_field!(config, 2);

        let sum = lhs + rhs;
        assert_eq!(sum.into_bigint(), BigInt::one());
    }

    #[test]
    fn test_add_without_wrapping() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 20);
        let rhs: RandomField<1> = create_random_field!(config, 20);

        let sum = lhs + rhs;
        assert_eq!(sum.into_bigint(), create_bigint!(17));
    }

    #[test]
    fn test_add_one() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 22);
        let rhs: RandomField<1> = RandomField::one();

        let sum = lhs + rhs;
        assert_eq!(sum.into_bigint(), BigInt::zero());
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
    fn test_sub_wrapping_around_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 2);
        let rhs: RandomField<1> = create_random_field!(config, 22);

        let difference = lhs - rhs;
        assert_eq!(difference.into_bigint(), create_bigint!(3));
    }

    #[test]
    fn test_sub_identical_values_results_in_zero() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 20);
        let rhs: RandomField<1> = create_random_field!(config, 20);

        let difference = lhs - rhs;
        assert_eq!(difference.into_bigint(), BigInt::zero());
    }

    #[test]
    fn test_init_sub_raw() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 2);
        let rhs = RandomField::one();
        let res = lhs - rhs;
        let mut expected = lhs;
        expected.set_one();
        assert_eq!(res, expected)
    }

    #[test]
    fn test_sub_assign_works() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let mut lhs: RandomField<1> = create_random_field!(config, 10);
        let rhs: RandomField<1> = create_random_field!(config, 7);

        lhs -= rhs;

        assert_eq!(lhs.into_bigint(), create_bigint!(3));
    }

    #[test]
    fn test_sub_assign_wraps_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let mut lhs: RandomField<1> = create_random_field!(config, 3);
        let rhs: RandomField<1> = create_random_field!(config, 7);

        lhs -= rhs;

        assert_eq!(lhs.into_bigint(), create_bigint!(19)); // 3 - 7 mod 23 = 19
    }

    #[test]
    fn test_mul_wraps_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 22);
        let rhs: RandomField<1> = create_random_field!(config, 2);

        let product = lhs * rhs;
        assert_eq!(product.into_bigint(), create_bigint!(21));
    }

    #[test]
    fn test_mul_without_wrapping() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 20);
        let rhs: RandomField<1> = create_random_field!(config, 20);

        let product = lhs * rhs;
        assert_eq!(product.into_bigint(), create_bigint!(9));
    }

    #[test]
    fn test_left_mul_by_zero() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 22);
        let rhs = RandomField::zero();

        let product = lhs * rhs;
        assert!(product.is_zero());
    }

    #[test]
    fn test_right_mul_by_zero() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let lhs: RandomField<1> = RandomField::zero();
        let rhs: RandomField<1> = create_random_field!(config, 22);

        let product = lhs * rhs;
        assert!(product.is_zero());
    }

    #[test]
    fn test_mul_assign_works() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let mut lhs: RandomField<1> = create_random_field!(config, 5);
        let rhs: RandomField<1> = create_random_field!(config, 4);

        lhs *= rhs;

        assert_eq!(lhs.into_bigint(), create_bigint!(20));
    }

    #[test]
    fn test_mul_assign_wraps_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let mut lhs: RandomField<1> = create_random_field!(config, 6);
        let rhs: RandomField<1> = create_random_field!(config, 4);

        lhs *= rhs;

        assert_eq!(lhs.into_bigint(), create_bigint!(1)); // 6 * 4 mod 23 = 1
    }

    #[test]
    fn test_div_wraps_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 22);
        let rhs: RandomField<1> = create_random_field!(config, 2);

        let quotient = lhs / rhs;
        assert_eq!(quotient.into_bigint(), create_bigint!(11));
    }

    #[test]
    fn test_div_identical_values_results_in_one() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 20);
        let rhs: RandomField<1> = create_random_field!(config, 20);

        let quotient = lhs / rhs;
        assert_eq!(quotient.into_bigint(), create_bigint!(1));
    }

    #[test]
    fn test_div_without_wrapping() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 17);
        let rhs: RandomField<1> = create_random_field!(config, 4);

        let quotient = lhs / rhs;
        assert_eq!(quotient.into_bigint(), create_bigint!(10));
    }

    #[test]
    #[should_panic]
    fn test_div_by_zero_should_panic() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 17);
        let rhs: RandomField<1> = create_random_field!(config, 0);

        let _sum = lhs / rhs;
    }

    #[test]
    fn test_div_bigint256() {
        let config =
            FieldConfig::new(BigInt::from_str("695962179703626800597079116051991347").unwrap());
        let config = ConfigRef::<4>::from(&config);

        let a: RandomField<4> = create_random_field!(config, 3);
        let mut b = RandomField::one();
        b /= a;
        assert_eq!(
            b.into_bigint(),
            create_bigint!(4, 231987393234542266865693038683997116)
        );

        let a: RandomField<4> = create_random_field!(config, 19382769832175);

        let b: RandomField<4> = create_random_field!(config, 97133987132135);

        assert_eq!(
            create_bigint!(4, 243043087159742188419721163456177516),
            (b / a).into_bigint()
        );
    }

    #[test]
    fn test_div_by_reference_works() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 15);
        let rhs = create_random_field!(config, 3);

        #[allow(clippy::op_ref)] // This implementation could be removed?
        let quotient = lhs / &rhs;

        assert_eq!(quotient.into_bigint(), create_bigint!(5));
    }

    #[test]
    fn test_div_by_mutable_reference_works() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let lhs: RandomField<1> = create_random_field!(config, 9);
        let rhs = create_random_field!(config, 3);

        #[allow(clippy::op_ref)] // This implementation could be removed?
        let quotient = lhs / &rhs;

        assert_eq!(quotient.into_bigint(), create_bigint!(3));
    }

    #[test]
    fn test_div_assign_works() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let mut lhs: RandomField<1> = create_random_field!(config, 15);
        let rhs: RandomField<1> = create_random_field!(config, 3);

        lhs /= rhs;

        assert_eq!(lhs.into_bigint(), create_bigint!(5));
    }

    #[test]
    #[should_panic(expected = "Attempt to divide by zero")]
    fn test_div_assign_by_zero_should_panic() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let mut lhs: RandomField<1> = create_random_field!(config, 15);
        let rhs = RandomField::zero();

        lhs /= rhs;
    }

    #[test]
    fn test_div_assign_by_mutable_reference() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let mut lhs: RandomField<1> = create_random_field!(config, 18);
        let mut rhs = create_random_field!(config, 3);

        lhs /= &mut rhs;

        assert_eq!(lhs.into_bigint(), create_bigint!(6)); // 18 / 3 mod 23 = 6
    }

    #[test]
    fn test_neg_large_value() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let operand: RandomField<1> = create_random_field!(config, 22);
        let negated = -operand;

        assert_eq!(negated.into_bigint(), create_bigint!(1));
    }

    #[test]
    fn test_neg_mid_value() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let operand: RandomField<1> = create_random_field!(config, 17);
        let negated = -operand;

        assert_eq!(negated.into_bigint(), create_bigint!(6));
    }

    #[test]
    fn test_neg_zero() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::<1>::from(&config);

        let operand: RandomField<1> = create_random_field!(config, 0);
        let negated = -operand;

        assert_eq!(negated.into_bigint(), BigInt::zero());
    }

    #[test]
    fn test_sum_of_multiple_values() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [
            create_random_field!(config, 2),
            create_random_field!(config, 4),
            create_random_field!(config, 6),
        ];

        let sum: RandomField<1> = values.iter().sum();

        assert_eq!(sum.into_bigint(), create_bigint!(12));
    }

    #[test]
    fn test_sum_with_zero() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [
            RandomField::zero(),
            create_random_field!(config, 5),
            create_random_field!(config, 7),
        ];

        let sum: RandomField<1> = values.iter().sum();

        assert_eq!(sum.into_bigint(), create_bigint!(12));
    }

    #[test]
    fn test_sum_wraps_modulus() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [
            create_random_field!(config, 10),
            create_random_field!(config, 15),
            create_random_field!(config, 21),
        ];

        let sum: RandomField<1> = values.iter().sum();

        assert_eq!(sum.into_bigint(), create_bigint!(0));
    }

    #[test]
    fn test_sum_empty_iterator() {
        let sum: RandomField<1> = ark_std::iter::empty::<&RandomField<1>>().sum();
        assert!(sum.is_zero()); // Empty sum should return zero
    }

    #[test]
    fn test_sum_single_element() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [create_random_field!(config, 9)];

        let sum: RandomField<1> = values.iter().sum();

        assert_eq!(sum.into_bigint(), create_bigint!(9));
    }

    #[test]
    fn test_sum_with_modulus_wrapping() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [
            create_random_field!(config, 12),
            create_random_field!(config, 15),
        ];

        let sum: RandomField<1> = values.iter().sum();

        assert_eq!(sum.into_bigint(), create_bigint!(4));
    }

    #[test]
    fn test_product_of_multiple_values() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [
            create_random_field!(config, 2),
            create_random_field!(config, 4),
            create_random_field!(config, 6),
        ];

        let product: RandomField<1> = values.iter().product();

        assert_eq!(product.into_bigint(), create_bigint!(2));
    }

    #[test]
    fn test_product_with_one() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [
            RandomField::one(),
            create_random_field!(config, 5),
            create_random_field!(config, 7),
        ];

        let product: RandomField<1> = values.iter().product();

        assert_eq!(product.into_bigint(), create_bigint!(12));
    }

    #[test]
    fn test_product_with_zero() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);
        let values = [
            create_random_field!(config, 3),
            RandomField::zero(),
            create_random_field!(config, 9),
        ];

        let product: RandomField<1> = values.iter().product();

        assert!(product.is_zero());
    }

    #[test]
    fn test_product_negative_modular_complements() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [
            create_random_field!(config, 10),
            create_random_field!(config, 15),
            create_random_field!(config, 21),
        ];

        let product: RandomField<1> = values.iter().product();

        assert_eq!(product.into_bigint(), create_bigint!(22));
    }

    #[test]
    fn test_product_empty_iterator() {
        let product: RandomField<1> = ark_std::iter::empty::<&RandomField<1>>().product();
        assert!(product.is_one()); // Empty product should return one
    }

    #[test]
    fn test_product_single_element() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [create_random_field!(config, 9)];

        let product: RandomField<1> = values.iter().product();

        assert_eq!(product.into_bigint(), create_bigint!(9));
    }

    #[test]
    fn test_product_with_modulus_wrapping() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);

        let values = [
            create_random_field!(config, 12),
            create_random_field!(config, 15),
        ];

        let product: RandomField<1> = values.iter().product();

        assert_eq!(product.into_bigint(), create_bigint!(19));
    }
}
