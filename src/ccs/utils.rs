//! Provides operations used for working with constraint systems
#![allow(non_snake_case)]

use std::ops::{Add, Mul};

use crate::sparse_matrix::SparseMatrix;

use super::error::CSError as Error;

//  Computes the hadamard product of two ring
#[allow(dead_code)]
pub(crate) fn hadamard_vec<R: Clone + Mul<R, Output = R>>(lhs: &[R], rhs: &[R]) -> Vec<R> {
    lhs.iter()
        .zip(rhs)
        .map(|(lhs, rhs)| lhs.clone() * rhs.clone())
        .collect()
}

// Multiplies Vector of rings by another ring
#[allow(dead_code)]
pub(crate) fn vec_value_mul<R: Clone + Mul<R, Output = R>>(lhs: &[R], rhs: &R) -> Vec<R> {
    lhs.iter()
        .map(|lhs_i| lhs_i.clone() * rhs.clone())
        .collect()
}

// Adds two ring vectors
pub(crate) fn vec_add<R: Clone + Add<R, Output = R>>(a: &[R], b: &[R]) -> Result<Vec<R>, Error> {
    if a.len() != b.len() {
        return Err(Error::LengthsNotEqual(
            "a".to_string(),
            "b".to_string(),
            a.len(),
            b.len(),
        ));
    }
    Ok(a.iter()
        .zip(b.iter())
        .map(|(x, y)| x.clone() + y.clone())
        .collect())
}

pub(crate) fn vec_scalar_mul<R: Clone + Mul<R, Output = R>>(vec: &[R], c: &R) -> Vec<R> {
    vec.iter().map(|a| a.clone() * c.clone()).collect()
}

pub(crate) fn hadamard<R: Clone + Mul<R, Output = R>>(a: &[R], b: &[R]) -> Result<Vec<R>, Error> {
    if a.len() != b.len() {
        return Err(Error::LengthsNotEqual(
            "a".to_string(),
            "b".to_string(),
            a.len(),
            b.len(),
        ));
    }
    Ok(a.iter()
        .zip(b)
        .map(|(a, b)| a.clone() * b.clone())
        .collect())
}

pub(crate) fn mat_vec_mul<R>(M: &SparseMatrix<R>, z: &[R]) -> Result<Vec<R>, Error>
where
    R: Clone + Send + Sync + Mul<R, Output = R> + Add<Output = R> + Default,
    for<'a> &'a R: Mul<&'a R, Output = R>,
{
    if M.n_cols != z.len() {
        return Err(Error::LengthsNotEqual(
            "M".to_string(),
            "z".to_string(),
            M.n_cols,
            z.len(),
        ));
    }

    let mut result = Vec::with_capacity(M.coeffs.len());

    for row in &M.coeffs {
        let mut acc = R::default(); // Assuming Default gives the additive identity (e.g., 0)
        for (value, col_i) in row {
            acc = acc + (value * &z[*col_i]);
        }
        result.push(acc);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;
    use crate::biginteger::BigInt;
    use crate::field::conversion::FieldMap;
    use crate::field_config::FieldConfig;
    use crate::sparse_matrix::dense_matrix_to_sparse;

    const N: usize = 3;
    fn get_config() -> FieldConfig<N> {
        FieldConfig::new(BigInt::<N>::from_str("695962179703626800597079116051991347").unwrap())
    }

    #[test]
    fn test_hadamard_vec() {
        let config = get_config();

        let a = [
            2u64.map_to_field(&config),
            3u64.map_to_field(&config),
            4u64.map_to_field(&config),
        ];
        let b = [
            5u64.map_to_field(&config),
            6u64.map_to_field(&config),
            7u64.map_to_field(&config),
        ];
        let result = hadamard_vec(&a, &b);
        let expected = vec![
            10u64.map_to_field(&config),
            18u64.map_to_field(&config),
            28u64.map_to_field(&config),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_vec_value_mul() {
        let config = get_config();

        let a = [
            2u64.map_to_field(&config),
            3u64.map_to_field(&config),
            4u64.map_to_field(&config),
        ];
        let scalar = 2u64.map_to_field(&config);
        let result = vec_value_mul(&a, &scalar);
        let expected = vec![
            4u64.map_to_field(&config),
            6u64.map_to_field(&config),
            8u64.map_to_field(&config),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_vec_add() {
        let config = get_config();

        let a = [
            1u64.map_to_field(&config),
            2u64.map_to_field(&config),
            3u64.map_to_field(&config),
        ];
        let b = [
            4u64.map_to_field(&config),
            5u64.map_to_field(&config),
            6u64.map_to_field(&config),
        ];
        let result = vec_add(&a, &b);
        let expected = vec![
            5u64.map_to_field(&config),
            7u64.map_to_field(&config),
            9u64.map_to_field(&config),
        ];
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_vec_add_error_case() {
        let config = get_config();

        let a = [1u64.map_to_field(&config), 2u64.map_to_field(&config)];
        let b = [
            3u64.map_to_field(&config),
            4u64.map_to_field(&config),
            5u64.map_to_field(&config),
        ];
        let result = vec_add(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_vec_scalar_mul() {
        let config = get_config();

        let vec = [
            1u64.map_to_field(&config),
            2u64.map_to_field(&config),
            3u64.map_to_field(&config),
        ];
        let c = 3u64.map_to_field(&config);
        let result = vec_scalar_mul(&vec, &c);
        let expected = vec![
            3u64.map_to_field(&config),
            6u64.map_to_field(&config),
            9u64.map_to_field(&config),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hadamard() {
        let config = get_config();

        let a = [
            2u64.map_to_field(&config),
            3u64.map_to_field(&config),
            4u64.map_to_field(&config),
        ];
        let b = [
            5u64.map_to_field(&config),
            6u64.map_to_field(&config),
            7u64.map_to_field(&config),
        ];
        let result = hadamard(&a, &b);
        let expected = vec![
            10u64.map_to_field(&config),
            18u64.map_to_field(&config),
            28u64.map_to_field(&config),
        ];
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_hadamard_error_case() {
        let config = get_config();

        let a = [2u64.map_to_field(&config), 3u64.map_to_field(&config)];
        let b = [
            5u64.map_to_field(&config),
            6u64.map_to_field(&config),
            7u64.map_to_field(&config),
        ];
        let result = hadamard(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_mat_vec_mul() {
        let config = get_config();

        let dense_matrix = vec![
            vec![
                1u64.map_to_field(&config),
                0u64.map_to_field(&config),
                0u64.map_to_field(&config),
            ],
            vec![
                0u64.map_to_field(&config),
                2u64.map_to_field(&config),
                1u64.map_to_field(&config),
            ],
            vec![
                0u64.map_to_field(&config),
                0u64.map_to_field(&config),
                3u64.map_to_field(&config),
            ],
        ];
        let M = dense_matrix_to_sparse(dense_matrix);

        let z = [
            1u64.map_to_field(&config),
            1u64.map_to_field(&config),
            1u64.map_to_field(&config),
        ];
        let result = mat_vec_mul(&M, &z);
        let expected = vec![
            1u64.map_to_field(&config),
            3u64.map_to_field(&config),
            3u64.map_to_field(&config),
        ];
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_mat_vec_mul_error_case() {
        let config = get_config();

        let dense_matrix = vec![
            vec![
                1u64.map_to_field(&config),
                0u64.map_to_field(&config),
                0u64.map_to_field(&config),
            ],
            vec![
                0u64.map_to_field(&config),
                2u64.map_to_field(&config),
                1u64.map_to_field(&config),
            ],
            vec![
                0u64.map_to_field(&config),
                0u64.map_to_field(&config),
                3u64.map_to_field(&config),
            ],
        ];
        let M = dense_matrix_to_sparse(dense_matrix);

        let z = [1u32.map_to_field(&config), 1u32.map_to_field(&config)];
        let result = mat_vec_mul(&M, &z);
        assert!(result.is_err());
    }
}
