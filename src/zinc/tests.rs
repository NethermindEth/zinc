use std::marker::PhantomData;

use ark_std::log2;
use crypto_bigint::{Int, U192, Zero, const_monty_params};
use rand::rng;
use rand_core::RngCore;

use crate::{
    ccs,
    ccs::{CcsF, CcsZ, Statement, Witness},
    field::{F192, WORD_FACTOR},
    sparse_matrix::{SparseMatrix, to_Z_matrix},
    traits::Field,
    transcript::KeccakTranscript,
    zinc::{
        prover::SpartanProver,
        structs::{ZincProver, ZincVerifier},
        verifier::SpartanVerifier,
    },
    zip::code::DefaultLinearCodeSpec,
};

const INT_LIMBS: usize = WORD_FACTOR;
const FIELD_LIMBS: usize = 3 * WORD_FACTOR;

const N: usize = INT_LIMBS;
const L: usize = INT_LIMBS * 2;
const K: usize = INT_LIMBS * 4;
const M: usize = INT_LIMBS * 8;

const_monty_params!(
    ModP,
    U192,
    "0000000000000000EB58CBFB80010B1661534393E6C4FA11"
);

type F = F192<ModP>;

#[test]
fn test_dummy_spartan_prover() {
    let n = 1 << 13;
    let mut rng = rng();

    let (_, ccs, statement, wit) = ccs::get_dummy_ccs_Z_from_z_length(n, &mut rng);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);

    let (z_ccs, z_mle, ccs_f, statement_f) = ZincProver::<
        N,
        L,
        K,
        M,
        F,
        FIELD_LIMBS,
        DefaultLinearCodeSpec,
    >::prepare_for_random_field_piop(
        &statement, &wit, &ccs
    )
    .expect("Failed to prepare for random field PIOP");

    let proof = SpartanProver::<INT_LIMBS, F, FIELD_LIMBS>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
    );

    assert!(proof.is_ok())
}

#[test]
fn test_spartan_verifier() {
    let input = 3;

    let (ccs, statement, wit, _) = get_test_ccs_stuff_Z(input);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);

    let (z_ccs, z_mle, ccs_f, statement_f) = ZincProver::<
        N,
        L,
        K,
        M,
        F,
        FIELD_LIMBS,
        DefaultLinearCodeSpec,
    >::prepare_for_random_field_piop(
        &statement, &wit, &ccs
    )
    .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<INT_LIMBS, F, FIELD_LIMBS>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);
    let mut verifier_transcript = KeccakTranscript::new();

    let res = SpartanVerifier::<F, FIELD_LIMBS>::verify(
        &verifier,
        &spartan_proof,
        &ccs_f,
        &mut verifier_transcript,
    );

    assert!(res.is_ok())
}

#[test]
fn test_dummy_spartan_verifier() {
    let n = 1 << 13;
    let mut rng = rng();

    let (_, ccs, statement, wit) = ccs::get_dummy_ccs_Z_from_z_length(n, &mut rng);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);

    let (z_ccs, z_mle, ccs_f, statement_f) = ZincProver::<
        N,
        L,
        K,
        M,
        F,
        FIELD_LIMBS,
        DefaultLinearCodeSpec,
    >::prepare_for_random_field_piop(
        &statement, &wit, &ccs
    )
    .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<INT_LIMBS, F, FIELD_LIMBS>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);
    let mut verifier_transcript = KeccakTranscript::new();
    let res = SpartanVerifier::<F, FIELD_LIMBS>::verify(
        &verifier,
        &spartan_proof,
        &ccs_f,
        &mut verifier_transcript,
    );

    assert!(res.is_ok(), "{res:?}");
}

#[test]
fn test_failing_spartan_verifier() {
    let input = 3;

    let (ccs, statement, mut wit, _) = get_test_ccs_stuff_Z(input);
    // Change the witness such that it is no longer valid
    assert_ne!(wit.w_ccs[3], Int::zero());
    wit.w_ccs[3] = Int::zero();
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);

    let (z_ccs, z_mle, ccs_f, statement_f) = ZincProver::<
        N,
        L,
        K,
        M,
        F,
        FIELD_LIMBS,
        DefaultLinearCodeSpec,
    >::prepare_for_random_field_piop(
        &statement, &wit, &ccs
    )
    .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<INT_LIMBS, F, FIELD_LIMBS>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);
    let mut verifier_transcript = KeccakTranscript::new();

    let res = SpartanVerifier::<F, FIELD_LIMBS>::verify(
        &verifier,
        &spartan_proof,
        &ccs_f,
        &mut verifier_transcript,
    );

    assert!(res.is_err())
}

#[allow(unused)]
pub(crate) fn create_dummy_identity_sparse_matrix_F<F: Field<LIMBS>, const LIMBS: usize>(
    rows: usize,
    columns: usize,
) -> SparseMatrix<F> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((F::one(), i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
#[allow(unused)]
pub(crate) fn create_dummy_squaring_sparse_matrix_F<F: Field<LIMBS>, const LIMBS: usize>(
    rows: usize,
    columns: usize,
    witness: &[F],
) -> SparseMatrix<F> {
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

#[allow(unused)]
fn get_dummy_ccs_F_from_z<F: Field<LIMBS>, const LIMBS: usize>(
    z: &[F],
    pub_io_len: usize,
) -> (CcsF<F>, Statement<F>, Witness<F>) {
    let ccs = CcsF {
        m: z.len(),
        n: z.len(),
        l: pub_io_len,
        t: 3,
        q: 2,
        d: 2,
        s: log2(z.len()) as usize,
        s_prime: log2(z.len()) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![F::one(), F::from(-1i64)],
    };

    let A = create_dummy_identity_sparse_matrix_F(z.len(), z.len());
    let B = A.clone();
    let C = create_dummy_squaring_sparse_matrix_F(z.len(), z.len(), z);

    let statement = Statement::<F> {
        constraints: vec![A, B, C],
        public_input: z[..pub_io_len].to_vec(),
    };

    let wit = Witness {
        w_ccs: z[pub_io_len + 1..].to_vec(),
    };

    (ccs, statement, wit)
}

#[allow(unused)]
pub fn get_dummy_ccs_F_from_z_length<F: Field<LIMBS>, const LIMBS: usize>(
    n: usize,
    rng: &mut impl RngCore,
) -> (Vec<F>, CcsF<F>, Statement<F>, Witness<F>) {
    let mut z: Vec<_> = (0..n).map(|_| F::random(rng)).collect();
    let pub_io_len = 1;
    z[pub_io_len] = F::one();

    let (ccs, statement, wit) = get_dummy_ccs_F_from_z(&z, pub_io_len);

    (z, ccs, statement, wit)
}

pub(crate) fn get_test_ccs_stuff_Z<const I: usize>(
    input: i64,
) -> (
    CcsZ<Int<I>>,
    Statement<Int<I>>,
    Witness<Int<I>>,
    Vec<Int<I>>,
) {
    let mut ccs = get_test_ccs_Z();
    let mut statement = get_test_ccs_Z_statement(input);
    let witness = get_test_wit_Z(input);
    let z = get_test_z_Z(input);
    let len = usize::max(ccs.m.next_power_of_two(), ccs.n.next_power_of_two());
    ccs.pad(&mut statement, len);
    (ccs, statement, witness, z)
}

pub(crate) fn get_test_ccs_Z<const I: usize>() -> CcsZ<Int<I>> {
    // R1CS for: x^3 + x + 5 = y (example from article
    // https://www.vitalik.ca/general/2016/12/10/qap.html )

    let m = 4;
    let n = 6;
    CcsZ {
        m,
        n,
        l: 1,
        t: 3,
        q: 2,
        d: 2,
        s: log2(m) as usize,
        s_prime: log2(n) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![1, -1],
        _phantom: PhantomData,
    }
}

pub(crate) fn get_test_ccs_Z_statement<const I: usize>(input: i64) -> Statement<Int<I>> {
    let A = to_Z_matrix(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 1, 0, 0],
        vec![1, 0, 0, 0, 1, 0],
        vec![0, 5, 0, 0, 0, 1],
    ]);
    let B = to_Z_matrix(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
    ]);
    let C = to_Z_matrix(vec![
        vec![0, 0, 0, 1, 0, 0],
        vec![0, 0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 0, 1],
        vec![0, 0, 1, 0, 0, 0],
    ]);
    let constraints = vec![A, B, C];
    let public_input = vec![Int::from_i64(input)];
    Statement {
        constraints,
        public_input,
    }
}

pub(crate) fn get_test_wit_Z<const I: usize>(input: i64) -> Witness<Int<I>> {
    Witness::new(vec![
        Int::from_i64(input * input * input + input + 5), // x^3 + x + 5
        Int::from_i64(input * input),                     // x^2
        Int::from_i64(input * input * input),             // x^2 * x
        Int::from_i64(input * input * input + input),     // x^3 + x
    ])
}

pub(crate) fn get_test_z_Z<const I: usize>(input: i64) -> Vec<Int<I>> {
    // z = (io, 1, w)
    vec![
        Int::from_i64(input), // io
        Int::from_i64(1),
        Int::from_i64(input * input * input + input + 5), // x^3 + x + 5
        Int::from_i64(input * input),                     // x^2
        Int::from_i64(input * input * input),             // x^2 * x
        Int::from_i64(input * input * input + input),     // x^3 + x
    ]
}
