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
    MODULUS_HAS_SPARE_BIT: bool,

    /// Can we use the no-carry optimization for multiplication
    /// outlined [here](https://hackmd.io/@gnark/modular_multiplication)?
    ///
    /// This optimization applies if
    /// (a) `Self::MODULUS[N-1] < u64::MAX >> 1`, and
    /// (b) the bits of the modulus are not all 1.
    CAN_USE_NO_CARRY_MUL_OPT: bool,

    /// Can we use the no-carry optimization for squaring
    /// outlined [here](https://hackmd.io/@gnark/modular_multiplication)?
    ///
    /// This optimization applies if
    /// (a) `Self::MODULUS[N-1] < u64::MAX >> 2`, and
    /// (b) the bits of the modulus are not all 1.
    CAN_USE_NO_CARRY_SQUARE_OPT: bool,

    /// 2^s root of unity computed by GENERATOR^t
    TWO_ADIC_ROOT_OF_UNITY: BigInt<N>,

    /// An integer `b` such that there exists a multiplicative subgroup
    /// of size `b^k` for some integer `k`.
    SMALL_SUBGROUP_BASE: Option<u32>,

    /// The integer `k` such that there exists a multiplicative subgroup
    /// of size `Self::SMALL_SUBGROUP_BASE^k`.
    SMALL_SUBGROUP_BASE_ADICITY: Option<u32>,

    /// GENERATOR^((MODULUS-1) / (2^s *
    /// SMALL_SUBGROUP_BASE^SMALL_SUBGROUP_BASE_ADICITY)).
    /// Used for mixed-radix FFT.
    LARGE_SUBGROUP_ROOT_OF_UNITY: Option<BigInt<N>>,
}

impl<const N: usize> MontConfig<N> {
    fn add_assign(&self, a: &mut BigInt<N>, b: &BigInt<N>) {
        // This cannot exceed the backing capacity.
        let c = a.add_with_carry(&b);
        // However, it may need to be reduced
        if self.MODULUS_HAS_SPARE_BIT {
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
        if self.MODULUS_HAS_SPARE_BIT {
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

    fn mul_assign(self, a: &mut BigInt<N>, b: &BigInt<N>) {
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

        if self.MODULUS_HAS_SPARE_BIT {
            if *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        } else {
            if carry || *a >= self.modulus {
                a.sub_with_borrow(&self.modulus);
            }
        }
    }

    fn inverse(self, a: &BigInt<N>) -> Option<BigInt<N>> {
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
        let mut b = self.r2; // Avoids unnecessary reduction step.
        let mut c = BigInt::<N>::zero();

        while u != one && v != one {
            while u.is_even() {
                u.div2();

                if b.is_even() {
                    b.div2();
                } else {
                    let carry = b.add_with_carry(&self.modulus);
                    b.div2();
                    if !self.MODULUS_HAS_SPARE_BIT && carry {
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
                    if !self.MODULUS_HAS_SPARE_BIT && carry {
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
