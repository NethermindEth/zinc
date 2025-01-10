use std::ops::Add;

use ark_ff::{BigInt, BigInteger};

use crate::field_config::{self, FieldConfig};

pub struct RandomField<'config, const N: usize> {
    pub config: &'config FieldConfig<N>,
    pub value: BigInt<N>,
}

impl<'config, const N: usize> RandomField<'config, N> {
    pub fn new_unchecked(config: &'config FieldConfig<N>, value: BigInt<N>) -> Self {
        RandomField { config, value }
    }

    pub fn from_bigint(config: &'config FieldConfig<N>, value: BigInt<N>) -> Option<Self> {
        if value.is_zero() {
            Some(Self::new_unchecked(config, value))
        } else if value >= config.modulus {
            None
        } else {
            let mut r = value;
            config.mul_assign(&mut r, &config.r2);
            Some(Self::new_unchecked(config, value))
        }
    }

    pub fn into_bigint(&self) -> BigInt<N> {
        let mut r = self.value.0;
        // Montgomery Reduction
        for i in 0..N {
            let k = r[i].wrapping_mul(self.config.inv);
            let mut carry = 0;

            field_config::mac_with_carry(r[i], k, self.config.modulus.0[0], &mut carry);
            for j in 1..N {
                r[(j + i) % N] = field_config::mac_with_carry(
                    r[(j + i) % N],
                    k,
                    self.config.modulus.0[j],
                    &mut carry,
                );
            }
            r[i % N] = carry;
        }

        BigInt::new(r)
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
        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.

        let config_ptr_lhs: *const FieldConfig<N> = self.config;
        let config_ptr_rhs: *const FieldConfig<N> = rhs.config;

        if config_ptr_lhs != config_ptr_rhs {
            panic!("cannot add field elements of different fields");
        }

        todo!()
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_add() {
        // TODO: fill this in
    }
}
