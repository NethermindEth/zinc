#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use crypto_bigint::{U256, const_monty_params};
use rand::rng;

use crate::{
    field::{RandomField, WORD_FACTOR},
    poly::{dense::DenseMultilinearExtension, mle::MultilinearExtension},
    transcript::KeccakTranscript,
    zip::{
        code::DefaultLinearCodeSpec, code_raa::RaaCode, pcs::structs::MultilinearZip,
        pcs_transcript::PcsTranscript,
    },
};

const INT_LIMBS: usize = WORD_FACTOR;
const FIELD_LIMBS: usize = 4 * WORD_FACTOR;

const N: usize = INT_LIMBS;
const L: usize = INT_LIMBS * 2;
const K: usize = INT_LIMBS * 4;
const M: usize = INT_LIMBS * 8;

type LC = RaaCode<N, L, K, M>;
type TestZip = MultilinearZip<N, L, K, M, LC>;

const_monty_params!(
    ModP,
    U256,
    "EB0E9F20F7BFC231327A11792F585AC6C20C74ACCCAB538BE6B0C3AB2E3D176F"
);
type F = RandomField<ModP, FIELD_LIMBS>;

fn one_run<const P: usize>() {
    let mut rng = rng();
    // Match the benchmark’s transcript usage for linear code construction
    let mut keccak_transcript = KeccakTranscript::new();
    let poly_size = 1 << P;
    let linear_code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
    let params = TestZip::setup(poly_size, linear_code);

    let poly = DenseMultilinearExtension::rand(P, &mut rng);
    let (data, commitment) = TestZip::commit::<F, FIELD_LIMBS>(&params, &poly).expect("commit");

    // Same point choice as the bench
    let point = vec![1i64; P];
    let eval = *poly.evaluations.last().expect("nonempty evals");

    // Prover produces a proof once (exactly as in the bench)
    let mut prover_tx = PcsTranscript::<F, FIELD_LIMBS>::new();
    TestZip::open(
        &params,
        &poly,
        &data,
        &point.iter().map(F::from).collect::<Vec<_>>(),
        &mut prover_tx,
    )
    .expect("open");
    let proof = prover_tx.into_proof();

    // Verifier replays verification from the same proof (also like the bench)
    let mut verifier_tx = PcsTranscript::<F, FIELD_LIMBS>::from_proof(&proof);
    TestZip::verify(
        &params,
        &commitment,
        &point.iter().map(F::from).collect::<Vec<_>>(),
        eval.resize().into(),
        &mut verifier_tx,
    )
    .expect("verify");
}

#[test]
fn zip_verify_equivalent_to_bench_once_p12() {
    // Mirrors: Zip/Verify: RandomField<4>, poly_size = 2^12 (Int limbs = 1)
    one_run::<12>();
}

#[test]
fn zip_verify_equivalent_to_bench_reuse_proof_multiple_times_p12() {
    // Bench calls verify repeatedly with the SAME proof; reproduce that here.
    // If prover/verifier transcripts diverge on column selection, this will fail with
    // InvalidPcsOpen(\"Proximity failure\") like your bench did.
    for _ in 0..10 {
        one_run::<12>();
    }
}
