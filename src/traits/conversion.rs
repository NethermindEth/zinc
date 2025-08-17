use ark_std::vec::Vec;

use crate::{field::RandomField, traits::ConfigReference};

/// A trait for converting from little-endian and big-endian byte slices into a concrete type.
pub trait FromBytes: Sized {
    /// Constructs an instance from a little-endian byte slice.
    fn from_bytes_le(bytes: &[u8]) -> Option<Self>;

    /// Constructs an instance from a big-endian byte slice.
    fn from_bytes_be(bytes: &[u8]) -> Option<Self>;
}

pub trait ToBytes {
    fn to_bytes(&self) -> Vec<u8>;
}

pub trait FieldMap<C: ConfigReference> {
    type Output;
    fn map_to_field(&self, config_ref: C) -> Self::Output;
}

pub trait MapsToField<C: ConfigReference>
where
    Self: FieldMap<C, Output = RandomField<C>>,
{
}

impl<C: ConfigReference, T> MapsToField<C> for T where T: FieldMap<C, Output = RandomField<C>> {}

pub trait FromRef<T>: for<'a> From<&'a T> {}

impl<T, U> FromRef<T> for U where U: for<'a> From<&'a T> {}
