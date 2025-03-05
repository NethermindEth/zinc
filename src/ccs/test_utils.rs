#![allow(dead_code, non_snake_case)]
use std::{ops::Neg, sync::atomic::AtomicPtr, vec};

use ark_std::{log2, rand::Rng};

use crate::{
    biginteger::BigInt,
    field::{rand_with_config, RandomField},
    field_config::FieldConfig,
    sparse_matrix::SparseMatrix,
};

use super::{
    ccs_f::{Statement_F, Witness_F, CCS_F},
    ccs_z::{Statement_Z, Witness_Z, CCS_Z},
};

pub(crate) fn create_dummy_identity_sparse_matrix(
    rows: usize,
    columns: usize,
) -> SparseMatrix<i128> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((1i128, i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
pub(crate) fn create_dummy_squaring_sparse_matrix(
    rows: usize,
    columns: usize,
    witness: &[i64],
) -> SparseMatrix<i128> {
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
        row.push((witness[i] as i128, i));
    }
    matrix
}

fn get_dummy_ccs_from_z(z: &[i64]) -> (CCS_Z, Statement_Z, Witness_Z) {
    let ccs = CCS_Z {
        m: z.len(),
        n: z.len(),
        l: 1,
        t: 3,
        q: 2,
        d: 2,
        s: log2(z.len()) as usize,
        s_prime: log2(z.len()) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![1, -1],
    };

    let A = create_dummy_identity_sparse_matrix(z.len(), z.len());
    let B = A.clone();
    let C = create_dummy_squaring_sparse_matrix(z.len(), z.len(), z);

    let statement = Statement_Z {
        constraints: vec![A, B, C],
        public_input: Vec::new(),
    };

    let wit = Witness_Z {
        w_ccs: z[ccs.l..].to_vec(),
    };

    (ccs, statement, wit)
}

pub fn get_dummy_ccs_from_z_length(
    n: usize,
    rng: &mut impl Rng,
) -> (Vec<i64>, CCS_Z, Statement_Z, Witness_Z) {
    let z: Vec<_> = (0..n).map(|_| rng.gen_range(i64::MIN..=i64::MAX)).collect();
    let (ccs, statement, wit) = get_dummy_ccs_from_z(&z);

    (z, ccs, statement, wit)
}
