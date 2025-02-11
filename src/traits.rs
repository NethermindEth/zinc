/// A trait for converting from little-endian and big-endian byte slices into a concrete type.
pub trait FromBytes: Sized {
    /// Constructs an instance from a little-endian byte slice.
    fn from_bytes_le(bytes: &[u8]) -> Option<Self>;

    /// Constructs an instance from a byte-endian byte slice.
    fn from_bytes_be(bytes: &[u8]) -> Option<Self>;
}
