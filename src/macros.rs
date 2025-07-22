#[macro_export]
macro_rules! new_bigint {
    ($s:expr) => {
        zinc::biginteger::BigInt::from_str(stringify!($s)).unwrap()
    };

    ($n:expr, $s:expr) => {
        zinc::biginteger::BigInt::<$n>::from_str(stringify!($s)).unwrap()
    };
}

#[macro_export]
macro_rules! new_config {
    ($s:expr) => {
        $crate::config::FieldConfig::new($crate::new_bigint!($s).unwrap())
    };

    ($n:expr, $s:expr) => {
        $crate::config::FieldConfig::<$n>::new($crate::new_bigint!($s).unwrap())
    };
}
