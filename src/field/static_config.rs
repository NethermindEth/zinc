pub trait StaticReference {
    type C;
    const CONFIG: Self::C;
}

#[macro_export]
macro_rules! impl_static_ref {
    ($name:ident, $n:expr, $modulus:expr) => {
        #[derive(Debug, Copy, Clone, Eq, PartialEq)]
        pub struct $name;

        impl $crate::field::StaticReference for $name {
            type C = $crate::field::FieldConfig<{ $n }>;

            const CONFIG: Self::C = $crate::field::FieldConfig::const_new($modulus);
        }

        impl $crate::traits::ConfigReference for $name {
            type C = $crate::field::FieldConfig<{ $n }>;
            type B = $crate::field::BigInt<{ $n }>;
            type I = $crate::field::Int<{ $n }>;
            type U = $crate::field::Uint<{ $n }>;
            type W = $crate::field::Words<{ $n }>;
            const N: usize = { $n };

            fn reference(&self) -> Option<&'static $crate::field::FieldConfig<{ $n }>> {
                use $crate::field::StaticReference;
                Some(&Self::CONFIG)
            }

            unsafe fn new(_: *mut $crate::field::FieldConfig<{ $n }>) -> Self {
                // Passionately do nothing.
                Self
            }

            fn pointer(&self) -> Option<*mut $crate::field::FieldConfig<{ $n }>> {
                use $crate::field::StaticReference;
                Some(&Self::CONFIG as *const _ as *mut _)
            }

            const NONE: Self = Self;
        }

        impl Default for $name {
            fn default() -> Self {
                Self
            }
        }
    };
}
