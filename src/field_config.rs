//! Finite field arithmetic using Montgomery reduction.
//!
//! This module defines the [`FieldConfig`] struct, which stores precomputed
//! constants for performing modular arithmetic efficiently in a finite field
//! modulo a given integer.
//!
//! Montgomery reduction is a method for accelerating modular multiplication by
//! avoiding costly division operations. This module precomputes constants
//! such as `R`, `R^2`, and `INV` to enable fast reductions for a fixed modulus.
use crate::biginteger::BigInt;

macro_rules! mac {
    ($a:expr, $b:expr, $c:expr, &mut $carry:expr$(,)?) => {{
        let tmp = ($a as u128) + widening_mul($b, $c);
        $carry = (tmp >> 64) as u64;
        tmp as u64
    }};
}
#[macro_export]
#[doc(hidden)]
macro_rules! mac_with_carry {
    ($a:expr, $b:expr, $c:expr, &mut $carry:expr$(,)?) => {{
        let tmp = ($a as u128) + widening_mul($b, $c) + ($carry as u128);
        $carry = (tmp >> 64) as u64;
        tmp as u64
    }};
}

macro_rules! adc {
    ($a:expr, $b:expr, &mut $carry:expr$(,)?) => {{
        let tmp = ($a as u128) + ($b as u128) + ($carry as u128);
        $carry = (tmp >> 64) as u64;
        tmp as u64
    }};
}
pub(super) fn mac_with_carry(a: u64, b: u64, c: u64, carry: &mut u64) -> u64 {
    let tmp = (a as u128) + widening_mul(b, c) + (*carry as u128);
    *carry = (tmp >> 64) as u64;
    tmp as u64
}
/// Defines a prime fie
#[derive(Clone, Copy, Default)]
/// Configuration parameters for a finite prime field.
///
/// This struct contains precomputed constants used for modular arithmetic
/// operations in Montgomery form. It is parameterized by `N`, the number
/// of 64-bit limbs representing the modulus.
pub struct FieldConfig<const N: usize> {
    /// The modulus of the field.
    pub modulus: BigInt<N>,

    /// Let `M` be the power of 2^64 nearest to `Self::MODULUS_BITS`. Then
    /// `R = M % Self::MODULUS`.
    pub r: BigInt<N>,

    /// R2 = R^2 % Self::MODULUS
    pub r2: BigInt<N>,

    /// INV = -MODULUS^{-1} mod 2^64
    pub inv: u64,

    /// Does the modulus have a spare unused bit
    ///
    /// This condition applies if
    /// (a) `Self::MODULUS[N-1] >> 63 == 0`
    #[doc(hidden)]
    modulus_has_spare_bit: bool,
}

impl<const N: usize> FieldConfig<N> {
    /// Creates a new [`FieldConfig`] instance for the given modulus.
    ///
    /// Precomputes the constants required for Montgomery arithmetic:
    /// - `r`: the Montgomery constant `R = 2^M mod modulus`
    /// - `r2`: the Montgomery constant `R^2 mod modulus`
    /// - `inv`: the modular inverse `INV = -modulus^{-1} mod 2^64`
    ///
    /// Also computes whether the modulus has a spare high bit.
    ///
    /// # Arguments
    /// - `modulus`: The modulus defining the field. Must be a prime number.
    ///
    /// # Example
    /// ```
    /// let modulus = BigInt::<2>::from(23u32);
    /// let config = FieldConfig::new(modulus);
    /// ```
    pub fn new(modulus: BigInt<N>) -> Self {
        let modulus_has_spare_bit = modulus.0[N - 1] >> 63 == 0;
        Self {
            modulus,
            r: modulus.montgomery_r(),
            r2: modulus.montgomery_r2(),
            inv: inv(modulus),

            modulus_has_spare_bit,
        }
    }

    pub(crate) fn add_assign(&self, a: &mut BigInt<N>, b: &BigInt<N>) {
        // This cannot exceed the backing capacity.
        let c = a.add_with_carry(b);
        // However, it may need to be reduced
        if self.modulus_has_spare_bit {
            if *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        } else if c || *a >= self.modulus {
            a.sub_with_borrow(&self.modulus);
        }
    }

    pub(crate) fn sub_assign(&self, a: &mut BigInt<N>, b: &BigInt<N>) {
        // If `other` is larger than `self`, add the modulus to self first.
        if b > a {
            a.add_with_carry(&self.modulus);
        }
        a.sub_with_borrow(b);
    }

    pub(crate) fn mul_assign(&self, a: &mut BigInt<N>, b: &BigInt<N>) {
        let (mut lo, mut hi) = ([0u64; N], [0u64; N]);
        crate::const_for!((i in 0..N) {
            let mut carry = 0;
            crate::const_for!((j in 0..N) {
                let k = i + j;
                if k >= N {
                    hi[k - N] = mac_with_carry!(hi[k - N], (a).0[i], (b).0[j], &mut carry);
                } else {
                    lo[k] = mac_with_carry!(lo[k], (a).0[i], (b).0[j], &mut carry);
                }
            });
            hi[i] = carry;
        });

        // Montgomery reduction
        let mut carry2 = 0;
        crate::const_for!((i in 0..N) {
            let tmp = lo[i].wrapping_mul(self.inv);
            let mut carry;
            mac!(lo[i], tmp, self.modulus.0[0], &mut carry);
            crate::const_for!((j in 1..N) {
                let k = i + j;
                if k >= N {
                    hi[k - N] = mac_with_carry!(hi[k - N], tmp, self.modulus.0[j], &mut carry);
                }  else {
                    lo[k] = mac_with_carry!(lo[k], tmp, self.modulus.0[j], &mut carry);
                }
            });
            hi[i] = adc!(hi[i], carry, &mut carry2);
        });

        crate::const_for!((i in 0..N) {
            (a).0[i] = hi[i];
        });

        let carry = carry2 != 0;

        if self.modulus_has_spare_bit {
            if *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        } else if carry || *a >= self.modulus {
            a.sub_with_borrow(&self.modulus);
        }
    }

    pub(crate) fn inverse(&self, a: &BigInt<N>) -> Option<BigInt<N>> {
        if a.is_zero() {
            return None;
        }

        // Guajardo Kumar Paar Pelzl
        // Efficient Software-Implementation of Finite Fields with Applications to
        // Cryptography
        // Algorithm 16 (BEA for Inversion in Fp)

        let one = BigInt::one();

        let mut u = *a;
        let mut v = self.modulus;
        let mut b = self.r2; // Avoids unnecessary reduction step.
        let mut c = BigInt::zero();

        while u != one && v != one {
            while u.is_even() {
                u.div2();

                if b.is_even() {
                    b.div2();
                } else {
                    let carry = b.add_with_carry(&self.modulus);
                    b.div2();
                    if !self.modulus_has_spare_bit && carry {
                        (b).0[N - 1] |= 1 << 63;
                    }
                }
            }

            while v.is_even() {
                v.div2();

                if c.is_even() {
                    c.div2();
                } else {
                    let carry = c.add_with_carry(&self.modulus);
                    c.div2();

                    if !self.modulus_has_spare_bit && carry {
                        (c).0[N - 1] |= 1 << 63;
                    }
                }
            }

            if v < u {
                u.sub_with_borrow(&v);

                if c > b {
                    b.add_with_carry(&self.modulus);
                }
                b.sub_with_borrow(&c);
            } else {
                v.sub_with_borrow(&u);

                if b > c {
                    c.add_with_carry(&self.modulus);
                }

                c.sub_with_borrow(&b);
            }
        }

        if u == one {
            Some(b)
        } else {
            Some(c)
        }
    }
}

impl<const N: usize> std::fmt::Debug for FieldConfig<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, " Z_{}", self.modulus,)
    }
}

/// Compute -M^{-1} mod 2^64.
pub const fn inv<const N: usize>(modulus: BigInt<N>) -> u64 {
    // We compute this as follows.
    // First, MODULUS mod 2^64 is just the lower 64 bits of MODULUS.
    // Hence MODULUS mod 2^64 = MODULUS.0[0] mod 2^64.
    //
    // Next, computing the inverse mod 2^64 involves exponentiating by
    // the multiplicative group order, which is euler_totient(2^64) - 1.
    // Now, euler_totient(2^64) = 1 << 63, and so
    // euler_totient(2^64) - 1 = (1 << 63) - 1 = 1111111... (63 digits).
    // We compute this powering via standard square and multiply.
    let mut inv = 1u64;
    crate::const_for!((_i in 0..63) {
        // Square
        inv = inv.wrapping_mul(inv);
        // Multiply
        inv = inv.wrapping_mul(modulus.0[0]);
    });
    inv.wrapping_neg()
}

fn widening_mul(a: u64, b: u64) -> u128 {
    a as u128 * b as u128
}

impl<const N: usize> PartialEq for FieldConfig<N> {
    fn eq(&self, other: &Self) -> bool {
        self.modulus == other.modulus
    }
}

impl<const N: usize> Eq for FieldConfig<N> {}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use crate::biginteger::{BigInteger128, BigInteger256};

    use super::FieldConfig;

    //BIGINTS ARE LITTLE ENDIAN!!
    #[test]
    fn test_addition() {
        let field = FieldConfig::new(BigInteger128::new([
            9307119299070690521,
            9320126393725433252,
        ]));
        let mut a = BigInteger128::new([2, 0]);
        let b = BigInteger128::new([2, 0]);
        field.add_assign(&mut a, &b);
        assert_eq!(a, BigInteger128::new([4, 0]));
    }

    #[test]
    fn test_subtraction() {
        let field = FieldConfig::new(BigInteger128::new([
            9307119299070690521,
            9320126393725433252,
        ]));
        let mut a = BigInteger128::new([2, 0]);
        let b = BigInteger128::new([2, 0]);
        field.sub_assign(&mut a, &b);
        assert_eq!(a, BigInteger128::zero());
    }

    #[test]
    fn test_multiplication() {
        let field = FieldConfig::new(
            BigInteger256::from_str("695962179703626800597079116051991347").unwrap(),
        );
        let mut a = BigInteger256::from_str("423024736033").unwrap();
        let b = BigInteger256::from_str("246308734").unwrap();
        field.mul_assign(&mut a, &b);

        assert_eq!(
            BigInteger256::from_str("504579159360957705315139767875358506").unwrap(),
            a
        );
    }
}
