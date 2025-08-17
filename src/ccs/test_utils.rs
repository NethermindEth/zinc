#![allow(non_snake_case)]

use ark_std::{log2, marker::PhantomData, rand::Rng, sync::atomic::AtomicPtr, vec, vec::Vec};

use super::{
    ccs_f::{CCS_F, Statement_F, Witness_F},
    ccs_z::{CCS_Z, Statement_Z, Witness_Z},
};
use crate::{
    field::RandomField,
    sparse_matrix::SparseMatrix,
    traits::{ConfigReference, FieldMap, Integer},
};

pub(crate) fn create_dummy_identity_sparse_matrix_Z<I: Integer>(
    rows: usize,
    columns: usize,
) -> SparseMatrix<I> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((I::one(), i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
pub(crate) fn create_dummy_squaring_sparse_matrix_Z<I: Integer>(
    rows: usize,
    columns: usize,
    witness: &[I],
) -> SparseMatrix<I> {
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
        row.push((witness[i].clone(), i));
    }
    matrix
}

pub(crate) fn create_dummy_identity_sparse_matrix_F<C: ConfigReference>(
    rows: usize,
    columns: usize,
    config: C,
) -> SparseMatrix<RandomField<C>> {
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
pub(crate) fn create_dummy_squaring_sparse_matrix_F<C: ConfigReference>(
    rows: usize,
    columns: usize,
    witness: &[RandomField<C>],
) -> SparseMatrix<RandomField<C>> {
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
        row.push((witness[i].clone(), i));
    }
    matrix
}

fn get_dummy_ccs_Z_from_z<I: Integer>(
    z: &[I],
    pub_io_len: usize,
) -> (CCS_Z<I>, Statement_Z<I>, Witness_Z<I>) {
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
        _phantom: PhantomData,
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

fn get_dummy_ccs_F_from_z<C: ConfigReference>(
    z: &[RandomField<C>],
    pub_io_len: usize,
    config: C,
) -> (CCS_F<C>, Statement_F<C>, Witness_F<C>) {
    let ccs = match config.pointer() {
        None => panic!("FieldConfig cannot be null"),
        Some(config_ptr) => CCS_F {
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

    let A = create_dummy_identity_sparse_matrix_F(z.len(), z.len(), config);
    let B = A.clone();
    let C = create_dummy_squaring_sparse_matrix_F(z.len(), z.len(), z);

    let statement = Statement_F::<C> {
        constraints: vec![A, B, C],
        public_input: z[..pub_io_len].to_vec(),
    };

    let wit = Witness_F {
        w_ccs: z[pub_io_len + 1..].to_vec(),
    };

    (ccs, statement, wit)
}

pub fn get_dummy_ccs_Z_from_z_length<I: Integer>(
    n: usize,
    rng: &mut impl Rng,
) -> (Vec<I>, CCS_Z<I>, Statement_Z<I>, Witness_Z<I>) {
    let mut z: Vec<_> = (0..n).map(|_| I::random(rng)).collect();
    let pub_io_len = 1;
    z[pub_io_len] = I::one();
    let (ccs, statement, wit) = get_dummy_ccs_Z_from_z(&z, pub_io_len);

    (z, ccs, statement, wit)
}

pub fn get_dummy_ccs_F_from_z_length<C: ConfigReference>(
    n: usize,
    rng: &mut impl Rng,
    config: C,
) -> (Vec<RandomField<C>>, CCS_F<C>, Statement_F<C>, Witness_F<C>) {
    let mut z: Vec<_> = (0..n)
        .map(|_| RandomField::rand_with_config(rng, config))
        .collect();
    let pub_io_len = 1;
    z[pub_io_len] = 1u64.map_to_field(config);

    let (ccs, statement, wit) = get_dummy_ccs_F_from_z(&z, pub_io_len, config);

    (z, ccs, statement, wit)
}
