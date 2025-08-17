//! Provides operations used for working with constraint systems
#![allow(non_snake_case)]

use ark_std::{
    ops::{Add, Mul},
    string::String,
    vec::Vec,
};

use super::error::CSError as Error;
use crate::sparse_matrix::SparseMatrix;

// Adds two ring vectors
pub(crate) fn vec_add<R: Clone + Add<R, Output = R>>(a: &[R], b: &[R]) -> Result<Vec<R>, Error> {
    if a.len() != b.len() {
        return Err(Error::LengthsNotEqual(
            String::from("a"),
            String::from("b"),
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
            "a".into(),
            "b".into(),
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
    for<'a> R: Mul<&'a R, Output = R>,
{
    if M.n_cols != z.len() {
        return Err(Error::LengthsNotEqual(
            "M".into(),
            "z".into(),
            M.n_cols,
            z.len(),
        ));
    }

    let mut result = Vec::with_capacity(M.coeffs.len());

    for row in &M.coeffs {
        let mut acc = R::default(); // Assuming Default gives the additive identity (e.g., 0)
        for (value, col_i) in row {
            acc = acc + (z[*col_i].clone() * value);
        }
        result.push(acc);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use ark_std::vec;

    use super::*;
    use crate::{
        big_int,
        field::{ConfigRef, FieldConfig},
        field_config,
        sparse_matrix::dense_matrix_to_sparse,
        traits::FieldMap,
    };

    const N: usize = 3;
    fn get_config() -> FieldConfig<N> {
        field_config!(695962179703626800597079116051991347, N)
    }

    #[test]
    fn test_vec_add() {
        let config = get_config();
        let config_ptr = ConfigRef::from(&config);

        let a = vec![1u64, 2u64, 3u64].map_to_field(config_ptr);
        let b = vec![4u64, 5u64, 6u64].map_to_field(config_ptr);
        let result = vec_add(&a, &b);
        let expected = vec![5u64, 7u64, 9u64].map_to_field(config_ptr);
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_vec_add_error_case() {
        let config = get_config();
        let config = ConfigRef::from(&config);

        let a = vec![1u64, 2].map_to_field(config);
        let b = vec![3u64, 4u64, 5u64].map_to_field(config);
        let result = vec_add(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_vec_scalar_mul() {
        let config = get_config();
        let config = ConfigRef::from(&config);

        let vec = vec![1u64, 2u64, 3u64].map_to_field(config);
        let c = 3u64.map_to_field(config);
        let result = vec_scalar_mul(&vec, &c);
        let expected = vec![3u64, 6u64, 9u64].map_to_field(config);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hadamard() {
        let config = get_config();
        let config = ConfigRef::from(&config);

        let a = vec![2u64, 3u64, 4u64].map_to_field(config);
        let b = vec![5u64, 6u64, 7u64].map_to_field(config);
        let result = hadamard(&a, &b);
        let expected = vec![10u64, 18u64, 28u64].map_to_field(config);
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_hadamard_error_case() {
        let config = get_config();
        let config = ConfigRef::from(&config);

        let a = vec![2u64, 3u64].map_to_field(config);
        let b = vec![5u64, 6u64, 7u64].map_to_field(config);
        let result = hadamard(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_mat_vec_mul() {
        let config = get_config();
        let config = ConfigRef::from(&config);

        let dense_matrix = vec![
            vec![1u64, 0u64, 0u64].map_to_field(config),
            vec![0u64, 2u64, 1u64].map_to_field(config),
            vec![0u64, 0u64, 3u64].map_to_field(config),
        ];
        let M = dense_matrix_to_sparse(dense_matrix);

        let z = vec![1u64, 1u64, 1u64].map_to_field(config);
        let result = mat_vec_mul(&M, &z);
        let expected = vec![1u64, 3u64, 3u64].map_to_field(config);
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_mat_vec_mul_error_case() {
        let config = get_config();
        let config = ConfigRef::from(&config);

        let dense_matrix = vec![
            vec![1u64, 0u64, 0u64].map_to_field(config),
            vec![0u64, 2u64, 1u64].map_to_field(config),
            vec![0u64, 0u64, 3u64].map_to_field(config),
        ];
        let M = dense_matrix_to_sparse(dense_matrix);

        let z = vec![1u32, 1u32].map_to_field(config);
        let result = mat_vec_mul(&M, &z);
        assert!(result.is_err());
    }
}
