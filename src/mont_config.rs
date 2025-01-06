use ark_ff::BigInt;

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

    /// Does the modulus have a spare unused bit
    ///
    /// This condition applies if
    /// (a) `Self::MODULUS[N-1] >> 63 == 0`
    MODULUS_HAS_SPARE_BIT: bool,

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
