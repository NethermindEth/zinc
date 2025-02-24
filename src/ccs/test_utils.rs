#![allow(dead_code, non_snake_case)]
use std::{ops::Neg, sync::atomic::AtomicPtr, vec};

use ark_std::{log2, rand::Rng};

use crate::{
    biginteger::BigInt,
    field::{rand_with_config, RandomField},
    field_config::FieldConfig,
    sparse_matrix::SparseMatrix,
};

use super::ccs_f::{Statement, Witness, CCS_F};

pub(crate) fn create_dummy_identity_sparse_matrix<const N: usize>(
    rows: usize,
    columns: usize,
    config: *const FieldConfig<N>,
) -> SparseMatrix<RandomField<N>> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((RandomField::from_bigint(config, 1u32.into()).unwrap(), i));
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

fn get_dummy_ccs_from_z<const N: usize>(
    z: &[RandomField<N>],
    config: *const FieldConfig<N>,
) -> (CCS_F<N>, Statement<N>, Witness<N>) {
    let ccs = CCS_F {
        m: z.len(),
        n: z.len(),
        l: 1,
        t: 3,
        q: 2,
        d: 2,
        s: log2(z.len()) as usize,
        s_prime: log2(z.len()) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![
            RandomField::from_bigint(config, BigInt::one()).unwrap(),
            RandomField::from_bigint(config, BigInt::one())
                .unwrap()
                .neg(),
        ],
        config: AtomicPtr::new(config as *mut FieldConfig<N>),
    };

    let A = create_dummy_identity_sparse_matrix(z.len(), z.len(), config);
    let B = A.clone();
    let C = create_dummy_squaring_sparse_matrix(z.len(), z.len(), z);

    let statement = Statement {
        constraints: vec![A, B, C],
        public_input: Vec::new(),
    };

    let wit = Witness {
        w_ccs: z[ccs.l..].to_vec(),
    };

    (ccs, statement, wit)
}

pub fn get_dummy_ccs_from_z_length<const N: usize>(
    n: usize,
    rng: &mut impl Rng,
    config: *const FieldConfig<N>,
) -> (Vec<RandomField<N>>, CCS_F<N>, Statement<N>, Witness<N>) {
    let z: Vec<_> = (0..n).map(|_| rand_with_config(rng, config)).collect();
    let (ccs, statement, wit) = get_dummy_ccs_from_z(&z, config);

    (z, ccs, statement, wit)
}
