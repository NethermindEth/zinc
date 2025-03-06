#![allow(clippy::not_unsafe_ptr_arg_deref, non_snake_case)]

use ark_ff::UniformRand;

use i256::I256;

use crate::{
    biginteger::BigInt,
    field_config::{self, FieldConfig},
};

pub mod arithmetic;
pub mod comparison;
pub mod constant;
pub mod conversion;

#[derive(Copy, Clone)]
pub enum RandomField<const N: usize> {
    Raw {
        value: BigInt<N>,
    },
    Initialized {
        config: *const FieldConfig<N>,
        value: BigInt<N>,
    },
}

use RandomField::*;

impl<const N: usize> RandomField<N> {
    pub fn is_raw(&self) -> bool {
        matches!(self, Raw { .. })
    }

    pub fn is_initialized(&self) -> bool {
        matches!(self, Initialized { .. })
    }

    pub fn with_raw_value_or<F, A>(&self, f: F, default: A) -> A
    where
        F: Fn(&BigInt<N>) -> A,
    {
        match self {
            Raw { value } => f(value),
            _ => default,
        }
    }

    pub fn with_raw_value_mut_or<F, A>(&mut self, f: F, default: A) -> A
    where
        F: Fn(&mut BigInt<N>) -> A,
    {
        match self {
            Raw { value } => f(value),
            _ => default,
        }
    }

    pub fn with_init_value<'a, F, A>(&'a self, f: F) -> Option<A>
    where
        F: Fn(&'a FieldConfig<N>, &'a BigInt<N>) -> A,
    {
        match self {
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");
                Some(f(config, value))
            },
            _ => None,
        }
    }

    pub fn with_init_value_or<'a, F, A>(&'a self, f: F, default: A) -> A
    where
        F: Fn(&'a FieldConfig<N>, &'a BigInt<N>) -> A,
    {
        match self {
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");
                f(config, value)
            },
            _ => default,
        }
    }

    pub fn with_either<'a, R, I, A>(&'a self, raw_fn: R, init_fn: I) -> A
    where
        I: Fn(&'a FieldConfig<N>, &'a BigInt<N>) -> A,
        R: Fn(&'a BigInt<N>) -> A,
    {
        match self {
            Raw { value } => raw_fn(value),
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                init_fn(config, value)
            },
        }
    }

    pub fn with_either_mut<'a, R, I, A>(&'a mut self, raw_fn: R, init_fn: I) -> A
    where
        I: Fn(&'a FieldConfig<N>, &'a mut BigInt<N>) -> A,
        R: Fn(&'a mut BigInt<N>) -> A,
    {
        match self {
            Raw { value } => raw_fn(value),
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                init_fn(config, value)
            },
        }
    }

    pub fn with_either_owned<R, I, A>(self, raw_fn: R, init_fn: I) -> A
    where
        I: Fn(&FieldConfig<N>, BigInt<N>) -> A,
        R: Fn(BigInt<N>) -> A,
    {
        match self {
            Raw { value } => raw_fn(value),
            Initialized { config, value } => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                init_fn(config, value)
            },
        }
    }

    #[inline(always)]
    pub fn value(&self) -> &BigInt<N> {
        match self {
            Raw { value } => value,
            Initialized { value, .. } => value,
        }
    }

    #[inline(always)]
    pub fn value_mut(&mut self) -> &mut BigInt<N> {
        match self {
            Raw { value } => value,
            Initialized { value, .. } => value,
        }
    }

    pub fn with_aligned_config_mut<F, G, A>(
        &mut self,
        rhs: &Self,
        with_config: F,
        without_config: G,
    ) -> A
    where
        F: Fn(&mut BigInt<N>, &BigInt<N>, &FieldConfig<N>) -> A,
        G: Fn(&mut BigInt<N>, &BigInt<N>) -> A,
    {
        match (self, rhs) {
            (Raw { value: value_self }, Raw { value: rhs }) => without_config(value_self, rhs),
            (
                Initialized {
                    value: value_self,
                    config,
                },
                Initialized {
                    value: value_rhs, ..
                },
            ) => unsafe {
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                with_config(value_self, value_rhs, config)
            },
            (
                Initialized {
                    value: value_self,
                    config,
                },
                rhs @ Raw { .. },
            ) => unsafe {
                let rhs = (*rhs).set_config_owned(*config);
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                with_config(value_self, rhs.value(), config)
            },
            (
                lhs @ Raw { .. },
                Initialized {
                    value: value_rhs,
                    config,
                },
            ) => unsafe {
                lhs.set_config(*config);
                let config = config
                    .as_ref()
                    .expect("Cannot have a null config for Initialized");

                with_config(lhs.value_mut(), value_rhs, config)
            },
        }
    }
}

impl<const N: usize> UniformRand for RandomField<N> {
    fn rand<R: ark_std::rand::Rng + ?Sized>(rng: &mut R) -> Self {
        let value = BigInt::rand(rng);

        Self::Raw { value }
    }
}

pub fn rand_with_config<const N: usize, R: ark_std::rand::Rng + ?Sized>(
    rng: &mut R,
    config: *const FieldConfig<N>,
) -> RandomField<N> {
    loop {
        let mut value = BigInt::rand(rng);
        let modulus = unsafe { (*config).modulus };
        let shave_bits = 64 * N - modulus.num_bits() as usize;
        // Mask away the unused bits at the beginning.
        assert!(shave_bits <= 64);
        let mask = if shave_bits == 64 {
            0
        } else {
            u64::MAX >> shave_bits
        };

        if let Some(val) = value.0.last_mut() {
            *val &= mask
        }

        if value < modulus {
            return RandomField::from_bigint(config, value).unwrap();
        }
    }
}

impl<const N: usize> RandomField<N> {
    pub fn set_config(&mut self, config: *const FieldConfig<N>) {
        self.with_raw_value_mut_or(
            |value| {
                // Ideally we should do something like:
                //
                // ```
                // let modulus: BigInt<N> = unsafe { (*config).modulus };
                // *value = *value % modulus;
                // ```
                //
                // but we don't have `mod` out of the box.
                // So let's hope we don't exceed the modulus.

                // TODO: prettify this

                *value = *Self::from_bigint(config, *value)
                    .expect("Should not end up with a None here.")
                    .value();
            },
            (),
        );

        let value = std::mem::take(self.value_mut());

        *self = Initialized { config, value }
    }

    /// Config setter that can be used after a `RandomField::rand(...)` call.
    pub fn set_config_owned(mut self, config: *const FieldConfig<N>) -> Self {
        self.set_config(config);
        self
    }
}

// TODO: Finalise this
//impl<const N: usize> AdditiveGroup for RandomField<N> {
//    type Scalar = ;
//
//    const ZERO: Self = Self { config: std::ptr::null(), BigInt::zero() };
//}

impl<const N: usize> RandomField<N> {
    #[inline(always)]
    pub fn config_ref(&self) -> Option<&FieldConfig<N>> {
        self.with_init_value(|config, _| config)
    }

    #[inline(always)]
    pub fn config_ptr(&self) -> *const FieldConfig<N> {
        match self {
            Raw { .. } => std::ptr::null(),
            Initialized { config, .. } => *config,
        }
    }
}

impl<const N: usize> RandomField<N> {
    pub fn new_unchecked(config: *const FieldConfig<N>, value: BigInt<N>) -> Self {
        Initialized { config, value }
    }

    /// Convert from `BigInteger` to `RandomField`
    ///
    /// If `BigInteger` is greater then field modulus return `None`
    pub fn from_bigint(config: *const FieldConfig<N>, value: BigInt<N>) -> Option<Self> {
        if config.is_null() {
            return Some(Raw { value });
        }

        unsafe {
            if value >= (*config).modulus {
                None
            } else {
                let mut r = value;

                (*config).mul_assign(&mut r, &(*config).r2);

                Some(Self::new_unchecked(config, r))
            }
        }
    }

    pub fn from_i64(value: i64, config: *const FieldConfig<N>) -> Option<RandomField<N>> {
        if config.is_null() {
            panic!("Cannot convert signed integer to prime field element without a modulus")
        }
        unsafe {
            if BigInt::from(value.unsigned_abs()) >= (*config).modulus {
                None
            } else {
                let mut r = (value.unsigned_abs()).into();

                (*config).mul_assign(&mut r, &(*config).r2);

                let mut elem = Self::new_unchecked(config, r);
                if value.is_negative() {
                    elem = -elem;
                }
                Some(elem)
            }
        }
    }
    // TODO This meethod is rather bad, needs to be improved
    pub(crate) fn from_I256(value: I256, config: *const FieldConfig<N>) -> RandomField<N> {
        if config.is_null() {
            panic!("Cannot convert signed integer to prime field element without a modulus")
        }
        unsafe {
            let mut modulus: [u64; N] = (*config).modulus.0;
            modulus.reverse();
            let val: [u64; 4] = value.abs().to_le_u64();
            let mut r = match N {
                n if n < 4 => {
                    let mut wider_modulus: [u64; 4] = [0; 4];
                    wider_modulus[..(4 - N)].copy_from_slice(&modulus);
                    let mut value = crypto_bigint::Uint::<4>::from_words(val);

                    let modu = crypto_bigint::Uint::<4>::from_words(wider_modulus);

                    value %= crypto_bigint::NonZero::from_uint(modu);

                    let mut result = [0u64; N];
                    result.copy_from_slice(&value.to_words()[..4 - N]);

                    BigInt(result)
                }
                4 => {
                    let mut value_N: [u64; N] = [0; N];
                    value_N.copy_from_slice(&val);

                    let mut value = crypto_bigint::Uint::<N>::from_words(value_N);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    value %= crypto_bigint::NonZero::from_uint(modu);
                    BigInt(value.to_words())
                }
                _ => {
                    let mut wider_value: [u64; N] = [0; N];
                    wider_value[(N - 4)..].copy_from_slice(&val);
                    let mut wider = crypto_bigint::Uint::<N>::from_words(wider_value);
                    let modu = crypto_bigint::Uint::<N>::from_words(modulus);
                    wider %= crypto_bigint::NonZero::from_uint(modu);
                    BigInt(wider.to_words())
                }
            };

            (*config).mul_assign(&mut r, &(*config).r2);

            let mut elem = Self::new_unchecked(config, r);
            if value.is_negative() {
                elem = -elem;
            }

            elem
        }
    }

    pub fn into_bigint(self) -> BigInt<N> {
        self.with_either_owned(|value| value, Self::demontgomery)
    }

    fn demontgomery(config: &FieldConfig<N>, value: BigInt<N>) -> BigInt<N> {
        let mut r = value.0;
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

impl<const N: usize> std::fmt::Debug for RandomField<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Raw { value } => write!(f, "{}, no config", value),
            self_ => write!(
                f,
                "{} in Z_{}",
                self_.into_bigint(),
                self.config_ref().unwrap().modulus
            ),
        }
    }
}

impl<const N: usize> std::fmt::Display for RandomField<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // TODO: we should go back from Montgomery here.
        match self {
            Raw { value } => {
                write!(f, "{}", value)
            }
            self_ @ Initialized { .. } => {
                write!(f, "{}", self_.into_bigint())
            }
        }
    }
}
impl<const N: usize> Default for RandomField<N> {
    fn default() -> Self {
        Raw {
            value: BigInt::zero(),
        }
    }
}

pub fn zero_with_config<const N: usize>(config: *const FieldConfig<N>) -> RandomField<N> {
    Initialized {
        config,
        value: BigInt::zero(),
    }
}

unsafe impl<const N: usize> Send for RandomField<N> {}
unsafe impl<const N: usize> Sync for RandomField<N> {}

#[cfg(test)]
mod tests {
    use crate::{biginteger::BigInt, field::RandomField, field_config::FieldConfig};
    use ark_std::str::FromStr;

    /// Helper macro to create a field config with a given modulus
    #[macro_export]
    macro_rules! create_field_config {
        ($N:expr, $modulus:expr) => {{
            FieldConfig::<$N>::new(BigInt::<$N>::from_str(stringify!($modulus)).unwrap())
        }};

        ($modulus:expr) => {{
            FieldConfig::<1>::new(BigInt::<1>::from_str(stringify!($modulus)).unwrap())
        }};
    }

    /// Helper macro to create a BigInt with a given modulus
    #[macro_export]
    macro_rules! create_bigint {
        ($N:expr, $value:expr) => {{
            BigInt::<$N>::from_str(stringify!($value)).unwrap()
        }};

        ($value:expr) => {{
            BigInt::from_str(stringify!($value)).unwrap()
        }};
    }

    /// Helper macro to create a RandomField with config and value.
    #[macro_export]
    macro_rules! create_random_field {
        ($config:expr, $value:expr) => {{
            RandomField::from_bigint($config, create_bigint!($value)).unwrap()
        }};
    }

    #[test]
    fn test_with_raw_value_or_for_raw_variant() {
        let raw_field = RandomField::<1>::Raw {
            value: create_bigint!(42),
        };

        assert_eq!(
            raw_field.with_raw_value_or(|v| *v, create_bigint!(99)),
            create_bigint!(42)
        );
    }

    #[test]
    fn test_with_raw_value_or_for_initialized_variant() {
        let config = create_field_config!(23);
        let init_field = RandomField::<1>::Initialized {
            config: &config,
            value: create_bigint!(10),
        };

        assert_eq!(
            init_field.with_raw_value_or(|v| *v, create_bigint!(99)),
            create_bigint!(99)
        );
    }
    #[test]
    fn test_with_init_value_or_initialized() {
        let config = create_field_config!(23);
        let init_field = RandomField::<1>::Initialized {
            config: &config,
            value: create_bigint!(10),
        };

        assert_eq!(
            init_field.with_init_value_or(|_, v| *v, create_bigint!(99)),
            create_bigint!(10)
        );
    }

    #[test]
    fn test_with_init_value_or_raw() {
        let raw_field = RandomField::<1>::Raw {
            value: create_bigint!(42),
        };

        assert_eq!(
            raw_field.with_init_value_or(|_, v| *v, create_bigint!(99)),
            create_bigint!(99)
        );
    }
}
