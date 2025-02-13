#![allow(dead_code)]
use ark_ff::One;

use crate::{field::RandomField, sparse_matrix::SparseMatrix};

pub(crate) fn create_dummy_identity_sparse_matrix<const N: usize>(
    rows: usize,
    columns: usize,
) -> SparseMatrix<RandomField<N>> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((RandomField::one(), i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
pub(crate) fn create_dummy_squaring_sparse_matrix<const N: usize>(
    rows: usize,
    columns: usize,
    witness: &[RandomField<N>],
) -> SparseMatrix<RandomField<N>> {
    assert_eq!(
        rows,
        witness.len(),
        "Length of witness vector must be equal to ccs width"
    );
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((witness[i], i));
    }
    matrix
}
