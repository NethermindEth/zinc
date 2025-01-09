#![allow(dead_code)]
use ark_ff::{BigInt, BigInteger};

macro_rules! mac {
    ($a:expr, $b:expr, $c:expr, &mut $carry:expr$(,)?) => {{
        let tmp = ($a as u128) + widening_mul($b, $c);
        $carry = (tmp >> 64) as u64;
        tmp as u64
    }};
}

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

pub struct MontConfig<const N: usize> {
    /// The modulus of the field.
    modulus: BigInt<N>,

    /// Let `M` be the power of 2^64 nearest to `Self::MODULUS_BITS`. Then
    /// `R = M % Self::MODULUS`.
    r: BigInt<N>,

    /// R2 = R^2 % Self::MODULUS
    r2: BigInt<N>,

    /// INV = -MODULUS^{-1} mod 2^64
    inv: u64,

    /// A multiplicative generator of the field.
    /// `Self::GENERATOR` is an element having multiplicative order
    /// `Self::MODULUS - 1`.
    generator: BigInt<N>,

    /// Does the modulus have a spare unused bit
    ///
    /// This condition applies if
    /// (a) `Self::MODULUS[N-1] >> 63 == 0`
    #[doc(hidden)]
    modulus_has_spare_bit: bool,
}

impl<const N: usize> MontConfig<N> {
    fn new(modulus: BigInt<N>, generator: BigInt<N>) -> Self {
        let modulus_has_spare_bit = modulus.0[N - 1] >> 63 == 0;
        Self {
            modulus,
            r: modulus.montgomery_r(),
            r2: modulus.montgomery_r2(),
            inv: inv(modulus),
            generator,
            modulus_has_spare_bit,
        }
    }

    fn add_assign(&self, a: &mut BigInt<N>, b: &BigInt<N>) {
        // This cannot exceed the backing capacity.
        let c = a.add_with_carry(&b);
        // However, it may need to be reduced
        if self.modulus_has_spare_bit {
            if *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        } else {
            if c || *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        }
    }

    fn sub_assign(&self, a: &mut BigInt<N>, b: &BigInt<N>) {
        // If `other` is larger than `self`, add the modulus to self first.
        if b > a {
            a.add_with_carry(&self.modulus);
        }
        a.sub_with_borrow(&b);
    }

    fn double_in_place(&self, a: &mut BigInt<N>) {
        // This cannot exceed the backing capacity.
        let carry = a.mul2();
        // However, it may need to be reduced.
        if self.modulus_has_spare_bit {
            if *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        } else {
            if carry || *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        }
    }

    /// Sets `a = -a`.
    #[inline(always)]
    fn neg_in_place(&self, a: &mut BigInt<N>) {
        if !a.is_zero() {
            let mut tmp = self.modulus;
            tmp.sub_with_borrow(&a);
            *a = tmp;
        }
    }

    fn mul_assign(&self, a: &mut BigInt<N>, b: &BigInt<N>) {
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
        } else {
            if carry || *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        }
    }

    pub fn inverse(&self, a: &BigInt<N>) -> Option<BigInt<N>> {
        if a.is_zero() {
            return None;
        }

        // Guajardo Kumar Paar Pelzl
        // Efficient Software-Implementation of Finite Fields with Applications to
        // Cryptography
        // Algorithm 16 (BEA for Inversion in Fp)

        let one = BigInt::from(1u64);

        let mut u = a.clone();
        let mut v = self.modulus;
        let mut b = one.clone(); // Avoids unnecessary reduction step.
        let mut c = BigInt::<N>::zero();

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
                b.sub_with_borrow(&c);
            } else {
                v.sub_with_borrow(&u);
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
#[macro_export]
macro_rules! const_for {
    (($i:ident in $start:tt..$end:tt)  $code:expr ) => {{
        let mut $i = $start;
        while $i < $end {
            $code
            $i += 1;
        }
    }};
}

fn widening_mul(a: u64, b: u64) -> u128 {
    a as u128 * b as u128
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInteger64;

    use super::MontConfig;

    #[test]
    fn test_addition() {
        let field = MontConfig::new(BigInteger64::from(83 as u64), BigInteger64::from(2 as u64));
        let mut a = BigInteger64::from(6 as u64);
        let b = BigInteger64::from(81 as u64);
        field.add_assign(&mut a, &b);
        assert_eq!(a, BigInteger64::from(4 as u32));
    }

    #[test]
    fn test_subtraction() {
        let field = MontConfig::new(BigInteger64::from(83 as u64), BigInteger64::from(2 as u64));
        let mut a = BigInteger64::from(5 as u64);
        let b = BigInteger64::from(8 as u64);
        field.sub_assign(&mut a, &b);
        assert_eq!(a, BigInteger64::from(80 as u32));
    }

    #[test]
    fn test_multiplication() {
        let field = MontConfig::new(BigInteger64::from(83 as u64), BigInteger64::from(2 as u64));

        let mut a = BigInteger64::from(4 as u64);
        let b = BigInteger64::from(8 as u64);
        field.mul_assign(&mut a, &b);
        assert_eq!(a, BigInteger64::from(47 as u32));

        let mut a = BigInteger64::from(2 as u64);
        let b = BigInteger64::from(7 as u64);
        field.mul_assign(&mut a, &b);
        assert_eq!(a, BigInteger64::from(5 as u32));
    }

    #[test]
    fn test_division() {
        let field = MontConfig::new(BigInteger64::from(83 as u64), BigInteger64::from(2 as u64));

        let a = BigInteger64::from(2 as u64);
        let b = field.inverse(&a).unwrap();
        assert_eq!(b, BigInteger64::from(42 as u32));

        let a = BigInteger64::from(3 as u64);
        let b = field.inverse(&a).unwrap();
        assert_eq!(b, BigInteger64::from(28 as u32));

        let a = BigInteger64::from(4 as u64);
        let b = field.inverse(&a).unwrap();
        assert_eq!(b, BigInteger64::from(21 as u32));
    }
}
