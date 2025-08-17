#![allow(non_snake_case)]

use ark_ff::UniformRand;
use ark_std::{format, string::String};
use crypto_bigint::Random;

use crate::traits::{Config, ConfigReference, FieldMap, FromBytes, FromRef, MapsToField};

mod arithmetic;
mod biginteger;
mod comparison;
mod config;
mod constant;
mod int;
mod uint;

pub use biginteger::{
    BigInt, BigInteger64, BigInteger128, BigInteger256, BigInteger320, BigInteger384,
    BigInteger448, BigInteger768, BigInteger832, Words, signed_mod_reduction,
};
pub use config::{ConfigRef, FieldConfig};
pub use int::Int;
pub use uint::Uint;
#[derive(Copy, Clone)]
pub enum RandomField<C: ConfigReference> {
    Raw { value: C::B },
    Initialized { config: C, value: C::B },
}

use RandomField::*;

use crate::{traits::BigInteger, transcript::KeccakTranscript};

impl<C: ConfigReference> RandomField<C> {
    pub fn is_raw(&self) -> bool {
        matches!(self, Raw { .. })
    }

    pub fn is_initialized(&self) -> bool {
        matches!(self, Initialized { .. })
    }

    pub fn with_raw_value_or<F, A>(&self, f: F, default: A) -> A
    where
        F: Fn(&C::B) -> A,
    {
        match self {
            Raw { value } => f(value),
            _ => default,
        }
    }

    pub fn with_raw_value_mut_or<F, A>(&mut self, f: F, default: A) -> A
    where
        F: Fn(&mut C::B) -> A,
    {
        match self {
            Raw { value } => f(value),
            _ => default,
        }
    }

    pub fn with_init_value<'a, F, A>(&'a self, f: F) -> Option<A>
    where
        F: Fn(&'a C::C, &'a C::B) -> A,
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
        F: Fn(&'a C::C, &'a C::B) -> A,
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
        I: Fn(&'a C::C, &'a C::B) -> A,
        R: Fn(&'a C::B) -> A,
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
        I: Fn(&'a C::C, &'a mut C::B) -> A,
        R: Fn(&'a mut C::B) -> A,
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
        I: Fn(&C::C, C::B) -> A,
        R: Fn(C::B) -> A,
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
        F: Fn(&mut C::B, &C::B, &C::C) -> A,
        G: Fn(&mut C::B, &C::B) -> A,
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
                let rhs = rhs.clone().set_config_owned(*config);
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

    pub fn config_copied(&self) -> Option<C::C> {
        match self {
            Raw { .. } => None,
            Initialized { config, .. } => config.reference().cloned(),
        }
    }

    pub fn zero_with_config(config: C) -> Self {
        Initialized {
            config,
            value: C::B::zero(),
        }
    }

    /// Config setter that can be used after a `RandomField::rand(...)` call.
    pub fn set_config_owned(mut self, config: C) -> Self {
        self.set_config(config);
        self
    }

    #[inline(always)]
    pub fn config_ptr(&self) -> C {
        match self {
            Raw { .. } => C::NONE,
            Initialized { config, .. } => *config,
        }
    }

    /// Convert from `BigInteger` to `RandomField`
    ///
    /// If `BigInteger` is greater then field modulus return `None`
    pub fn from_bigint(config: C, value: C::B) -> Option<RandomField<C>> {
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

    pub fn from_i64(value: i64, config: C) -> Option<RandomField<C>> {
        let config_ref = match config.reference() {
            Some(config) => config,
            None => {
                panic!("Cannot convert signed integer to prime field element without a modulus")
            }
        };

        if &C::B::from(value.unsigned_abs()) >= config_ref.modulus() {
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
    pub fn into_bigint(self) -> C::B {
        self.with_either_owned(|value| value, Self::demontgomery)
    }

    #[inline]
    fn demontgomery(config: &C::C, value: C::B) -> C::B {
        value.demontgomery(config.modulus(), config.inv())
    }
}

impl<C: ConfigReference> RandomField<C> {
    pub fn new_unchecked(config: C, value: C::B) -> Self {
        Initialized { config, value }
    }

    pub fn without_config(value: C::B) -> Self {
        Raw { value }
    }

    pub fn rand_with_config<R: ark_std::rand::Rng + ?Sized>(rng: &mut R, config: C) -> Self {
        loop {
            let mut value = C::B::rand(rng);
            let modulus = config
                .reference()
                .expect("Field config cannot be none")
                .modulus();
            let shave_bits = 64 * C::N - modulus.num_bits() as usize;
            // Mask away the unused bits at the beginning.
            assert!(shave_bits <= 64);
            let mask = if shave_bits == 64 {
                0
            } else {
                u64::MAX >> shave_bits
            };

            let val = value.last_mut();
            *val &= mask;

            if &value < modulus {
                return value.map_to_field(config);
            }
        }
    }

    pub fn set_config(&mut self, config: C) {
        self.with_raw_value_mut_or(
            |value| {
                // Ideally we should do something like:
                //
                // ```
                // let modulus: C::B = unsafe { (*config).modulus };
                // *value = *value % modulus;
                // ```
                //
                // but we don't have `mod` out of the box.
                // So let's hope we don't exceed the modulus.

                // TODO: prettify this

                *value = Self::from_bigint(config, value.clone())
                    .expect("Should not end up with a None here.")
                    .value()
                    .clone();
            },
            (),
        );

        let value = ark_std::mem::take(self.value_mut());

        *self = Initialized { config, value }
    }

    #[inline(always)]
    pub fn value(&self) -> &C::B {
        match self {
            Raw { value } => value,
            Initialized { value, .. } => value,
        }
    }

    #[inline(always)]
    pub fn value_mut(&mut self) -> &mut C::B {
        match self {
            Raw { value } => value,
            Initialized { value, .. } => value,
        }
    }

    pub fn absorb_into_transcript(&self, transcript: &mut KeccakTranscript) {
        match self {
            Raw { value } => {
                transcript.absorb(&[0x1]);
                transcript.absorb(&value.clone().to_bytes_be());
                transcript.absorb(&[0x3])
            }
            Initialized { config, value } => {
                let config = config.reference().expect("Field config cannot be none");

                transcript.absorb(&[0x3]);
                transcript.absorb(&config.modulus().clone().to_bytes_be());
                transcript.absorb(&[0x5]);

                transcript.absorb(&[0x1]);
                transcript.absorb(&value.clone().to_bytes_be());
                transcript.absorb(&[0x3])
            }
        }
    }
}

impl<C: ConfigReference> UniformRand for RandomField<C> {
    fn rand<R: ark_std::rand::Rng + ?Sized>(rng: &mut R) -> Self {
        let value = C::B::rand(rng);

        Self::Raw { value }
    }
}

impl<C: ConfigReference> Random for RandomField<C> {
    fn random(rng: &mut (impl ark_std::rand::RngCore + ?Sized)) -> Self {
        let value = C::B::rand(rng);

        Self::Raw { value }
    }
}

impl<C: ConfigReference> ark_std::fmt::Debug for RandomField<C> {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        match self {
            Raw { value } => write!(f, "{value:?}, no config"),
            self_ => write!(
                f,
                "{:?} in Z_{:?}",
                self_.clone().into_bigint(),
                self.config_ptr().reference().unwrap().modulus()
            ),
        }
    }
}

impl<C: ConfigReference> ark_std::fmt::Display for RandomField<C> {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        // TODO: we should go back from Montgomery here.
        match self {
            Raw { value } => {
                write!(f, "{value:?}")
            }
            self_ @ Initialized { .. } => {
                write!(f, "{:?}", self_.clone().into_bigint())
            }
        }
    }
}

impl<C: ConfigReference> Default for RandomField<C> {
    fn default() -> Self {
        Raw {
            value: C::B::zero(),
        }
    }
}

unsafe impl<C: ConfigReference> Send for RandomField<C> {}
unsafe impl<C: ConfigReference> Sync for RandomField<C> {}

impl<C: ConfigReference> From<u128> for RandomField<C> {
    fn from(value: u128) -> Self {
        let value = C::B::from(value);

        Raw { value }
    }
}

macro_rules! impl_from_uint {
    ($type:ty) => {
        impl<C: ConfigReference> From<$type> for RandomField<C> {
            fn from(value: $type) -> Self {
                let value = C::B::from(value);
                Raw { value }
            }
        }
    };
}

impl_from_uint!(u64);
impl_from_uint!(u32);
impl_from_uint!(u16);
impl_from_uint!(u8);

impl<C: ConfigReference> From<bool> for RandomField<C> {
    fn from(value: bool) -> Self {
        let value = C::B::from(value as u8);
        Raw { value }
    }
}

impl<C: ConfigReference> FromBytes for RandomField<C> {
    fn from_bytes_le(bytes: &[u8]) -> Option<Self> {
        Some(Raw {
            value: C::B::from_bytes_le(bytes)?,
        })
    }

    fn from_bytes_be(bytes: &[u8]) -> Option<Self> {
        Some(Raw {
            value: C::B::from_bytes_be(bytes)?,
        })
    }
}

impl<C: ConfigReference> From<RandomField<C>> for String {
    fn from(value: RandomField<C>) -> Self {
        format!("{value:?}")
    }
}

impl<C: ConfigReference> RandomField<C> {
    pub fn from_bytes_le_with_config(config: C, bytes: &[u8]) -> Option<Self> {
        let value = C::B::from_bytes_le(bytes);

        Self::from_bigint(config, value?)
    }

    pub fn from_bytes_be_with_config(config: C, bytes: &[u8]) -> Option<Self> {
        let value = C::B::from_bytes_be(bytes);

        Self::from_bigint(config, value?)
    }
}

// Implementation of FieldMap for C::B
impl<C: ConfigReference, const M: usize> FieldMap<C> for BigInt<M>
where
    C::I: FromRef<BigInt<M>>,
    Int<M>: FromRef<C::B>,
    C::B: From<Int<M>>,
{
    type Output = RandomField<C>;

    fn map_to_field(&self, config_ref: C) -> Self::Output {
        let config = match config_ref.reference() {
            Some(config) => config,
            None => panic!("Cannot convert BigInt to prime field element without a modulus"),
        };

        let mut value = if M > C::N {
            let modulus: Int<M> = config.modulus().into();
            let mut value: Int<M> = self.into();
            value %= modulus;

            C::B::from(value)
        } else {
            let modulus: C::I = config.modulus().into();
            let mut value: C::I = self.into();
            value %= modulus;

            value.into()
        };

        config.mul_assign(&mut value, config.r2());

        RandomField::new_unchecked(config_ref, value)
    }
}

// Implementation of FieldMap for reference to C::B
impl<C: ConfigReference, const M: usize> FieldMap<C> for &BigInt<M>
where
    BigInt<M>: MapsToField<C>,
{
    type Output = RandomField<C>;
    fn map_to_field(&self, config_ref: C) -> Self::Output {
        (*self).map_to_field(config_ref)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        big_int,
        field::{RandomField, config::ConfigRef},
        field_config,
    };

    #[test]
    fn test_with_raw_value_or_for_raw_variant() {
        let raw_field = RandomField::<ConfigRef<1>>::Raw {
            value: big_int!(42),
        };

        assert_eq!(
            raw_field.with_raw_value_or(|v| *v, big_int!(99)),
            big_int!(42)
        );
    }

    #[test]
    fn test_with_raw_value_or_for_initialized_variant() {
        let config = field_config!(23);
        let config = ConfigRef::<1>::from(&config);
        let init_field = RandomField::Initialized {
            config,
            value: big_int!(10),
        };

        assert_eq!(
            init_field.with_raw_value_or(|v| *v, big_int!(99)),
            big_int!(99)
        );
    }
    #[test]
    fn test_with_init_value_or_initialized() {
        let config = field_config!(23);
        let config = ConfigRef::<1>::from(&config);
        let init_field = RandomField::Initialized {
            config,
            value: big_int!(10),
        };

        assert_eq!(
            init_field.with_init_value_or(|_, v| *v, big_int!(99)),
            big_int!(10)
        );
    }

    #[test]
    fn test_with_init_value_or_raw() {
        let raw_field = RandomField::<ConfigRef<1>>::Raw {
            value: big_int!(42),
        };

        assert_eq!(
            raw_field.with_init_value_or(|_, v| *v, big_int!(99)),
            big_int!(99)
        );
    }
}
