use ark_ff::BigInt;

pub struct MontConfig<const N: usize> {
    modulus: BigInt<N>,

    r: BigInt<N>,
    r2: BigInt<N>,
    inv: u64,

    generator: BigInt<N>,
    CAN_USE_NO_CARRY_MUL_OPT: bool,
    CAN_USE_NO_CARRY_SQUARE_OPT: bool,
    MODULUS_HAS_SPARE_BIT: bool,
    TWO_ADIC_ROOT_OF_UNITY: BigInt<N>,
    SMALL_SUBGROUP_BASE: Option<u32>,
    SMALL_SUBGROUP_BASE_ADICITY: Option<u32>,
}
