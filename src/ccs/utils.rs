//! Provides operations used for working with constraint systems
#![allow(non_snake_case)]
#[cfg(feature = "parallel")]
use rayon::iter::*;
use std::{
    iter::Sum,
    ops::{Add, Mul},
};

use ark_std::cfg_iter;

use crate::{field_config::FieldConfig, sparse_matrix::SparseMatrix};

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
    R: Clone + Send + Sync + Mul<R, Output = R> + for<'a> Sum<R>,
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

    Ok(cfg_iter!(M.coeffs)
        .map(|row| row.iter().map(|(value, col_i)| value * &z[*col_i]).sum())
        .collect())
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::Zero;

    use super::*;
    use crate::biginteger::BigInt;
    use crate::field::RandomField;
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
            RandomField::from_bigint(&config, 2u64.into()).unwrap(),
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
            RandomField::from_bigint(&config, 4u64.into()).unwrap(),
        ];
        let b = [
            RandomField::from_bigint(&config, 5u64.into()).unwrap(),
            RandomField::from_bigint(&config, 6u64.into()).unwrap(),
            RandomField::from_bigint(&config, 7u64.into()).unwrap(),
        ];
        let result = hadamard_vec(&a, &b);
        let expected = vec![
            RandomField::from_bigint(&config, 10u64.into()).unwrap(),
            RandomField::from_bigint(&config, 18u64.into()).unwrap(),
            RandomField::from_bigint(&config, 28u64.into()).unwrap(),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_vec_value_mul() {
        let config = get_config();

        let a = [
            RandomField::from_bigint(&config, 2u64.into()).unwrap(),
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
            RandomField::from_bigint(&config, 4u64.into()).unwrap(),
        ];
        let scalar = RandomField::from_bigint(&config, 2u64.into()).unwrap();
        let result = vec_value_mul(&a, &scalar);
        let expected = vec![
            RandomField::from_bigint(&config, 4u64.into()).unwrap(),
            RandomField::from_bigint(&config, 6u64.into()).unwrap(),
            RandomField::from_bigint(&config, 8u64.into()).unwrap(),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_vec_add() {
        let config = get_config();

        let a = [
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            RandomField::from_bigint(&config, 2u64.into()).unwrap(),
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
        ];
        let b = [
            RandomField::from_bigint(&config, 4u64.into()).unwrap(),
            RandomField::from_bigint(&config, 5u64.into()).unwrap(),
            RandomField::from_bigint(&config, 6u64.into()).unwrap(),
        ];
        let result = vec_add(&a, &b);
        let expected = vec![
            RandomField::from_bigint(&config, 5u64.into()).unwrap(),
            RandomField::from_bigint(&config, 7u64.into()).unwrap(),
            RandomField::from_bigint(&config, 9u64.into()).unwrap(),
        ];
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_vec_add_error_case() {
        let config = get_config();

        let a = [
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            RandomField::from_bigint(&config, 2u64.into()).unwrap(),
        ];
        let b = [
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
            RandomField::from_bigint(&config, 4u64.into()).unwrap(),
            RandomField::from_bigint(&config, 5u64.into()).unwrap(),
        ];
        let result = vec_add(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_vec_scalar_mul() {
        let config = get_config();

        let vec = [
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            RandomField::from_bigint(&config, 2u64.into()).unwrap(),
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
        ];
        let c = RandomField::from_bigint(&config, 3u64.into()).unwrap();
        let result = vec_scalar_mul(&vec, &c);
        let expected = vec![
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
            RandomField::from_bigint(&config, 6u64.into()).unwrap(),
            RandomField::from_bigint(&config, 9u64.into()).unwrap(),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hadamard() {
        let config = get_config();

        let a = [
            RandomField::from_bigint(&config, 2u64.into()).unwrap(),
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
            RandomField::from_bigint(&config, 4u64.into()).unwrap(),
        ];
        let b = [
            RandomField::from_bigint(&config, 5u64.into()).unwrap(),
            RandomField::from_bigint(&config, 6u64.into()).unwrap(),
            RandomField::from_bigint(&config, 7u64.into()).unwrap(),
        ];
        let result = hadamard(&a, &b);
        let expected = vec![
            RandomField::from_bigint(&config, 10u64.into()).unwrap(),
            RandomField::from_bigint(&config, 18u64.into()).unwrap(),
            RandomField::from_bigint(&config, 28u64.into()).unwrap(),
        ];
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_hadamard_error_case() {
        let config = get_config();

        let a = [
            RandomField::from_bigint(&config, 2u64.into()).unwrap(),
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
        ];
        let b = [
            RandomField::from_bigint(&config, 5u64.into()).unwrap(),
            RandomField::from_bigint(&config, 6u64.into()).unwrap(),
            RandomField::from_bigint(&config, 7u64.into()).unwrap(),
        ];
        let result = hadamard(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_mat_vec_mul() {
        let config = get_config();

        let dense_matrix = vec![
            vec![
                RandomField::from_bigint(&config, 1u64.into()).unwrap(),
                RandomField::zero(),
                RandomField::zero(),
            ],
            vec![
                RandomField::zero(),
                RandomField::from_bigint(&config, 2u64.into()).unwrap(),
                RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            ],
            vec![
                RandomField::zero(),
                RandomField::zero(),
                RandomField::from_bigint(&config, 3u64.into()).unwrap(),
            ],
        ];
        let M = dense_matrix_to_sparse(dense_matrix);

        let z = [
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
        ];
        let result = mat_vec_mul(&M, &z);
        let expected = vec![
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
            RandomField::from_bigint(&config, 3u64.into()).unwrap(),
        ];
        assert_eq!(result.unwrap(), expected);
    }

    #[test]
    fn test_mat_vec_mul_error_case() {
        let config = get_config();

        let dense_matrix = vec![
            vec![
                RandomField::from_bigint(&config, 1u64.into()).unwrap(),
                RandomField::zero(),
                RandomField::zero(),
            ],
            vec![
                RandomField::zero(),
                RandomField::from_bigint(&config, 2u64.into()).unwrap(),
                RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            ],
            vec![
                RandomField::zero(),
                RandomField::zero(),
                RandomField::from_bigint(&config, 3u64.into()).unwrap(),
            ],
        ];
        let M = dense_matrix_to_sparse(dense_matrix);

        let z = [
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
            RandomField::from_bigint(&config, 1u64.into()).unwrap(),
        ];
        let result = mat_vec_mul(&M, &z);
        assert!(result.is_err());
    }
}