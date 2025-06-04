use num_traits::{Signed as NumTraitsSigned, Unsigned as NumTraitsUnsigned};

/// Marker trait for signed integer types that links to the corresponding unsigned type.
pub trait Signed: NumTraitsSigned {
    /// The associated unsigned type (e.g., `i32::Unsigned` is `u32`).
    type Unsigned: Unsigned;
}

pub trait Abs {
    type Unsigned: Unsigned + Copy;
    type Signed: Signed + Copy;
    fn unsigned_abs(self) -> Self::Unsigned;

    fn is_negative(&self) -> bool;
}

/// Marker trait for unsigned integer types that links to the corresponding signed type.
pub trait Unsigned: NumTraitsUnsigned {
    fn as_array<const N: usize>(&self) -> [u64; N];
    fn limbs() -> usize;
}

/// Implements [`Signed`] and [`Unsigned`] traits for pairs of signed and unsigned integer types.
macro_rules! impl_signed_and_unsigned {
    ($(($s:ty, $u:ty, $b:expr))*) => {
        $(
            impl Signed for $s {
                type Unsigned = $u;
            }

            impl Abs for $s {
                type Unsigned = $u;
                type Signed = $s;

                #[inline]
                fn unsigned_abs(self) -> Self::Unsigned {
                    (self as $s).unsigned_abs()
                }

                #[inline]
                fn is_negative(&self) -> bool {
                    (*self as $s).is_negative()
                }
            }

            impl Abs for $u {
                type Unsigned = $u;
                type Signed = $s;

                #[inline]
                fn unsigned_abs(self) -> Self::Unsigned {
                    self
                }

                #[inline]
                fn is_negative(&self) -> bool {
                    false
                }
            }

            impl Unsigned for $u {
                fn as_array<const N: usize>(&self) -> [u64; N] {
                    let mut res = [0u64; N];

                    res[0] = *self as u64;

                    if size_of::<Self>() > 8 && N > 1 {
                        res[1] = (self >> 64) as u64;
                    }

                    res
                }

                fn limbs() -> usize {
                    (size_of::<Self>() + 7) / 8
                }
            }
        )*
    };
}

impl_signed_and_unsigned!((isize, usize, std::mem::size_of::<usize>() * 8usize)(
    i8, u8, 8
)(i16, u16, 16)(i32, u32, 32)(i64, u64, 64)(i128, u128, 128));
