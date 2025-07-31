#![allow(non_snake_case)]

use ark_ff::UniformRand;
use crypto_bigint::{Int, Random, Uint};

use crate::{
    biginteger::BigInt,
    field_config::FieldConfig,
    traits::{Config, ConfigReference, FieldMap},
};

pub mod arithmetic;
pub mod comparison;
pub mod constant;
pub mod conversion;

#[derive(Copy, Clone)]
pub enum RandomField<'cfg, const N: usize> {
    Raw {
        value: BigInt<N>,
    },
    Initialized {
        config: ConfigRef<'cfg, N>,
        value: BigInt<N>,
    },
}

use RandomField::*;

use crate::{
    biginteger::Words,
    field_config::{ConfigRef, DebugFieldConfig},
    traits::{Field, Integer},
    transcript::KeccakTranscript,
};

impl<'cfg, const N: usize> RandomField<'cfg, N> {
    pub type Cfg = ConfigRef<'cfg, N>;

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
            Initialized { config, value } => Some(f(
                config.reference().expect("Field config cannot be none"),
                value,
            )),
            _ => None,
        }
    }

    pub fn with_init_value_or<'a, F, A>(&'a self, f: F, default: A) -> A
    where
        F: Fn(&'a FieldConfig<N>, &'a BigInt<N>) -> A,
    {
        match self {
            Initialized { config, value } => f(
                config.reference().expect("Field config cannot be none"),
                value,
            ),
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
            Initialized { config, value } => init_fn(
                config.reference().expect("Field config cannot be none"),
                value,
            ),
        }
    }

    pub fn with_either_mut<'a, R, I, A>(&'a mut self, raw_fn: R, init_fn: I) -> A
    where
        I: Fn(&'a FieldConfig<N>, &'a mut BigInt<N>) -> A,
        R: Fn(&'a mut BigInt<N>) -> A,
    {
        match self {
            Raw { value } => raw_fn(value),
            Initialized { config, value } => init_fn(
                config.reference().expect("Field config cannot be none"),
                value,
            ),
        }
    }

    pub fn with_either_owned<R, I, A>(self, raw_fn: R, init_fn: I) -> A
    where
        I: Fn(&FieldConfig<N>, BigInt<N>) -> A,
        R: Fn(BigInt<N>) -> A,
    {
        match self {
            Raw { value } => raw_fn(value),
            Initialized { config, value } => init_fn(
                config.reference().expect("Field config cannot be none"),
                value,
            ),
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
            ) => with_config(
                value_self,
                value_rhs,
                config.reference().expect("Field config cannot be none"),
            ),
            (
                Initialized {
                    value: value_self,
                    config,
                },
                rhs @ Raw { .. },
            ) => {
                let rhs = (*rhs).set_config_owned(*config);
                with_config(
                    value_self,
                    rhs.value(),
                    config.reference().expect("Field config cannot be none"),
                )
            }
            (
                lhs @ Raw { .. },
                Initialized {
                    value: value_rhs,
                    config,
                },
            ) => {
                lhs.set_config(*config);

                with_config(
                    lhs.value_mut(),
                    value_rhs,
                    config.reference().expect("Field config cannot be none"),
                )
            }
        }
    }

    pub fn config_copied(&self) -> Option<FieldConfig<N>> {
        match self {
            Raw { .. } => None,
            Initialized { config, .. } => config.reference().copied(),
        }
    }

    pub fn zero_with_config(config: Self::Cfg) -> Self {
        Initialized {
            config,
            value: BigInt::zero(),
        }
    }

    /// Config setter that can be used after a `RandomField::rand(...)` call.
    pub fn set_config_owned(mut self, config: Self::Cfg) -> Self {
        self.set_config(config);
        self
    }

    #[inline(always)]
    pub fn config_ptr(&self) -> Self::Cfg {
        match self {
            Raw { .. } => ConfigRef::NONE,
            Initialized { config, .. } => *config,
        }
    }

    /// Convert from `BigInteger` to `RandomField`
    ///
    /// If `BigInteger` is greater then field modulus return `None`
    pub fn from_bigint(config: ConfigRef<N>, value: BigInt<N>) -> Option<RandomField<N>> {
        let config_ref = match config.reference() {
            Some(config) => config,
            None => return Some(Raw { value }),
        };

        if value >= *config_ref.modulus() {
            None
        } else {
            let mut r = value;
            config_ref.mul_assign(&mut r, config_ref.r2());

            Some(RandomField::new_unchecked(config, r))
        }
    }

    pub fn from_i64(value: i64, config: ConfigRef<N>) -> Option<RandomField<N>> {
        let config_ref = match config.reference() {
            Some(config) => config,
            None => {
                panic!("Cannot convert signed integer to prime field element without a modulus")
            }
        };

        if BigInt::from(value.unsigned_abs()) >= *config_ref.modulus() {
            None
        } else {
            let mut r = value.unsigned_abs().into();
            config_ref.mul_assign(&mut r, config_ref.r2());

            let mut elem = RandomField::new_unchecked(config, r);
            if value.is_negative() {
                elem = -elem;
            }
            Some(elem)
        }
    }

    #[inline]
    pub fn into_bigint(self) -> BigInt<N> {
        self.with_either_owned(|value| value, Self::demontgomery)
    }

    #[inline]
    fn demontgomery(config: &FieldConfig<N>, value: BigInt<N>) -> BigInt<N> {
        value.demontgomery(config.modulus(), config.inv())
    }
}

impl<'cfg, const N: usize> Field for RandomField<'cfg, N> {
    type I = BigInt<N>;
    type C = FieldConfig<N>;
    type Cr = ConfigRef<'cfg, N>;
    type W = Words<N>;
    type CryptoInt = Int<N>;
    type CryptoUint = Uint<N>;
    type DebugField = DebugRandomField;

    fn new_unchecked(config: ConfigRef<'cfg, N>, value: BigInt<N>) -> Self {
        Initialized { config, value }
    }

    fn without_config(value: Self::I) -> Self {
        Raw { value }
    }

    fn rand_with_config<R: ark_std::rand::Rng + ?Sized>(rng: &mut R, config: Self::Cr) -> Self {
        loop {
            let mut value = BigInt::rand(rng);
            let modulus = config
                .reference()
                .expect("Field config cannot be none")
                .modulus();
            let shave_bits = 64 * N - modulus.num_bits() as usize;
            // Mask away the unused bits at the beginning.
            assert!(shave_bits <= 64);
            let mask = if shave_bits == 64 {
                0
            } else {
                u64::MAX >> shave_bits
            };

            let val = value.last_mut();
            *val &= mask;

            if value < *modulus {
                return value.map_to_field(config);
            }
        }
    }

    fn set_config(&mut self, config: Self::Cr) {
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

        let value = ark_std::mem::take(self.value_mut());

        *self = Initialized { config, value }
    }

    #[inline(always)]
    fn value(&self) -> &BigInt<N> {
        match self {
            Raw { value } => value,
            Initialized { value, .. } => value,
        }
    }

    #[inline(always)]
    fn value_mut(&mut self) -> &mut BigInt<N> {
        match self {
            Raw { value } => value,
            Initialized { value, .. } => value,
        }
    }

    fn absorb_into_transcript(&self, transcript: &mut KeccakTranscript) {
        match self {
            Raw { value } => {
                transcript.absorb(&[0x1]);
                transcript.absorb(&value.to_bytes_be());
                transcript.absorb(&[0x3])
            }
            Initialized { config, value } => {
                let config = config.reference().expect("Field config cannot be none");

                transcript.absorb(&[0x3]);
                transcript.absorb(&config.modulus().to_bytes_be());
                transcript.absorb(&[0x5]);

                transcript.absorb(&[0x1]);
                transcript.absorb(&value.to_bytes_be());
                transcript.absorb(&[0x3])
            }
        }
    }
}

impl<const N: usize> UniformRand for RandomField<'_, N> {
    fn rand<R: ark_std::rand::Rng + ?Sized>(rng: &mut R) -> Self {
        let value = BigInt::rand(rng);

        Self::Raw { value }
    }
}

impl<const N: usize> Random for RandomField<'_, N> {
    fn random(rng: &mut (impl ark_std::rand::RngCore + ?Sized)) -> Self {
        let value = BigInt::rand(rng);

        Self::Raw { value }
    }
}

impl<const N: usize> ark_std::fmt::Debug for RandomField<'_, N> {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        match self {
            Raw { value } => write!(f, "{}, no config", value),
            self_ => write!(
                f,
                "{} in Z_{}",
                self_.into_bigint(),
                self.config_ptr().reference().unwrap().modulus()
            ),
        }
    }
}

impl<const N: usize> ark_std::fmt::Display for RandomField<'_, N> {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
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

impl<const N: usize> Default for RandomField<'_, N> {
    fn default() -> Self {
        Raw {
            value: BigInt::zero(),
        }
    }
}

unsafe impl<const N: usize> Send for RandomField<'_, N> {}
unsafe impl<const N: usize> Sync for RandomField<'_, N> {}

#[derive(Debug)]
pub enum DebugRandomField {
    Raw {
        value: num_bigint::BigInt,
    },

    Initialized {
        config: DebugFieldConfig,
        value: num_bigint::BigInt,
    },
}

impl<const N: usize> From<RandomField<'_, N>> for DebugRandomField {
    fn from(value: RandomField<'_, N>) -> Self {
        match value {
            Raw { value } => Self::Raw {
                value: value.into(),
            },
            Initialized { config, value } => Self::Initialized {
                config: (*config.reference().unwrap()).into(),
                value: value.into(),
            },
        }
    }
}

impl ark_std::fmt::Display for DebugRandomField {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        match self {
            Self::Raw { value } => {
                write!(f, "{}", value)
            }
            self_ @ Self::Initialized { .. } => {
                write!(f, "{}", self_)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_std::str::FromStr;

    use crate::{
        biginteger::BigInt,
        field::RandomField,
        field_config::{ConfigRef, FieldConfig},
        traits::Config,
    };

    /// Helper macro to create a field config with a given modulus
    #[macro_export]
    macro_rules! create_field_config {
        ($N:expr, $modulus:expr) => {{
            let bigint = BigInt::<$N>::from_str(stringify!($modulus))
                .expect("Failed to parse modulus into BigInt");
            let cfg = FieldConfig::new(bigint);
            (cfg, ConfigRef::from(&cfg))
        }};

        ($modulus:expr) => {{
            let bigint = BigInt::<1>::from_str(&$modulus.to_string())
                .expect("Failed to parse modulus into BigInt");

            let cfg = FieldConfig::new(bigint);
            (cfg, ConfigRef::from(&cfg))
        }};
    }

    /// Helper macro to create a BigInt with a given modulus
    #[macro_export]
    macro_rules! create_bigint {
        ($N:expr, $value:expr) => {{
            BigInt::<$N>::from_str(stringify!($value)).unwrap()
        }};

        ($value:expr) => {{
            BigInt::<1>::from_str(stringify!($value)).unwrap()
        }};
    }

    /// Helper macro to create a RandomField with config and value.
    #[macro_export]
    macro_rules! create_random_field {
        ($config:expr, $value:expr) => {{
            use $crate::traits::FieldMap;
            create_bigint!($value).map_to_field($config)
        }};
    }

    #[test]
    fn test_with_raw_value_or_for_raw_variant() {
        let raw_field: RandomField<'_, 1> = RandomField::Raw {
            value: create_bigint!(42),
        };

        assert_eq!(
            raw_field.with_raw_value_or(|v| *v, create_bigint!(99)),
            create_bigint!(42)
        );
    }

    #[test]
    fn test_with_raw_value_or_for_initialized_variant() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);
        let init_field: RandomField<'_, 1> = RandomField::Initialized {
            config,
            value: create_bigint!(10),
        };

        assert_eq!(
            init_field.with_raw_value_or(|v| *v, create_bigint!(99)),
            create_bigint!(99)
        );
    }
    #[test]
    fn test_with_init_value_or_initialized() {
        let config = FieldConfig::new(BigInt::from_str("23").unwrap());
        let config = ConfigRef::from(&config);
        let init_field: RandomField<'_, 1> = RandomField::Initialized {
            config,
            value: create_bigint!(10),
        };

        assert_eq!(
            init_field.with_init_value_or(|_, v| *v, create_bigint!(99)),
            create_bigint!(10)
        );
    }

    #[test]
    fn test_with_init_value_or_raw() {
        let raw_field: RandomField<'_, 1> = RandomField::Raw {
            value: create_bigint!(42),
        };

        assert_eq!(
            raw_field.with_init_value_or(|_, v| *v, create_bigint!(99)),
            create_bigint!(99)
        );
    }
}
