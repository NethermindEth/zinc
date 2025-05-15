#![allow(non_snake_case)]
use std::{sync::atomic::AtomicPtr, vec};

use ark_ff::One;
use ark_std::{log2, rand::Rng};
use crypto_bigint::{Int, Random};

use super::{
    ccs_f::{Statement_F, Witness_F, CCS_F},
    ccs_z::{Statement_Z, Witness_Z, CCS_Z},
};
use crate::field_config::ConfigPtr;
use crate::{
    field::{conversion::FieldMap, rand_with_config, RandomField},
    sparse_matrix::SparseMatrix,
};

pub(crate) fn create_dummy_identity_sparse_matrix_Z<const N: usize>(
    rows: usize,
    columns: usize,
) -> SparseMatrix<Int<N>> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((Int::one(), i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
pub(crate) fn create_dummy_squaring_sparse_matrix_Z<const N: usize>(
    rows: usize,
    columns: usize,
    witness: &[Int<N>],
) -> SparseMatrix<Int<N>> {
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

pub(crate) fn create_dummy_identity_sparse_matrix_F<const N: usize>(
    rows: usize,
    columns: usize,
    config: ConfigPtr<N>,
) -> SparseMatrix<RandomField<N>> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((1u64.map_to_field(config), i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
pub(crate) fn create_dummy_squaring_sparse_matrix_F<const N: usize>(
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

fn get_dummy_ccs_Z_from_z<const N: usize>(
    z: &[Int<N>],
    pub_io_len: usize,
) -> (CCS_Z<N>, Statement_Z<N>, Witness_Z<N>) {
    let ccs = CCS_Z {
        m: z.len(),
        n: z.len(),
        l: pub_io_len,
        t: 3,
        q: 2,
        d: 2,
        s: log2(z.len()) as usize,
        s_prime: log2(z.len()) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![1, -1],
    };

    let A = create_dummy_identity_sparse_matrix_Z(z.len(), z.len());
    let B = A.clone();
    let C = create_dummy_squaring_sparse_matrix_Z(z.len(), z.len(), z);

    let statement = Statement_Z {
        constraints: vec![A, B, C],
        public_input: z[..pub_io_len].to_vec(),
    };

    let wit = Witness_Z {
        w_ccs: z[pub_io_len + 1..].to_vec(),
    };

    (ccs, statement, wit)
}

fn get_dummy_ccs_F_from_z<const N: usize>(
    z: &[RandomField<N>],
    pub_io_len: usize,
    config: ConfigPtr<N>,
) -> (CCS_F<N>, Statement_F<N>, Witness_F<N>) {
    let ccs = match config.pointer() {
        None => panic!("FieldConfig cannot be null"),
        Some(config_ptr) => CCS_F::<N> {
            m: z.len(),
            n: z.len(),
            l: pub_io_len,
            t: 3,
            q: 2,
            d: 2,
            s: log2(z.len()) as usize,
            s_prime: log2(z.len()) as usize,
            S: vec![vec![0, 1], vec![2]],
            c: vec![1u32.map_to_field(config), (-1i32).map_to_field(config)],
            config: AtomicPtr::new(config_ptr),
        },
    };

    let A = create_dummy_identity_sparse_matrix_F::<N>(z.len(), z.len(), config);
    let B = A.clone();
    let C = create_dummy_squaring_sparse_matrix_F::<N>(z.len(), z.len(), z);

    let statement = Statement_F::<N> {
        constraints: vec![A, B, C],
        public_input: z[..pub_io_len].to_vec(),
    };

    let wit = Witness_F::<N> {
        w_ccs: z[pub_io_len + 1..].to_vec(),
    };

    (ccs, statement, wit)
}

pub fn get_dummy_ccs_Z_from_z_length<const N: usize>(
    n: usize,
    rng: &mut impl Rng,
) -> (Vec<Int<N>>, CCS_Z<N>, Statement_Z<N>, Witness_Z<N>) {
    let mut z: Vec<_> = (0..n).map(|_| Int::<N>::random(rng)).collect();
    let pub_io_len = 1;
    z[pub_io_len] = Int::<N>::one();
    let (ccs, statement, wit) = get_dummy_ccs_Z_from_z(&z, pub_io_len);

    (z, ccs, statement, wit)
}

pub fn get_dummy_ccs_F_from_z_length<const N: usize>(
    n: usize,
    rng: &mut impl Rng,
    config: ConfigPtr<N>,
) -> (Vec<RandomField<N>>, CCS_F<N>, Statement_F<N>, Witness_F<N>) {
    let mut z: Vec<_> = (0..n).map(|_| rand_with_config(rng, config)).collect();
    let pub_io_len = 1;
    z[pub_io_len] = 1u64.map_to_field(config);
    let (ccs, statement, wit) = get_dummy_ccs_F_from_z(&z, pub_io_len, config);

    (z, ccs, statement, wit)
}
