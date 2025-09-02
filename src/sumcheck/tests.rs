use crypto_bigint::{U128, const_monty_params};
use num_traits::{ConstZero, One, Zero};
use rand::rng;
use rand_core::RngCore;

use super::{
    IPForMLSumcheck, MLSumcheck, SumcheckProof,
    utils::{rand_poly, rand_poly_comb_fn},
};
use crate::{
    field::{F128, WORD_FACTOR},
    poly::dense::DenseMultilinearExtension,
    sumcheck::prover::ProverState,
    traits::Field,
    transcript::KeccakTranscript,
};

const N: usize = 2 * WORD_FACTOR;

const_monty_params!(ModP, U128, "00000000B933426489189CB5B47D567F");

type F = F128<ModP>;

fn generate_sumcheck_proof<F: Field<LIMBS>, const LIMBS: usize, Rn: RngCore + ?Sized>(
    num_vars: usize,
    mut rng: &mut Rn,
) -> (usize, F, SumcheckProof<F>) {
    let mut transcript = KeccakTranscript::default();

    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );
    (poly_degree, sum, proof)
}
#[test]
fn full_sumcheck_protocol_works_correctly() {
    let num_vars = 3;
    let mut rn = rng();
    for _ in 0..20 {
        let (poly_degree, sum, proof) =
            generate_sumcheck_proof::<F, { 2 * WORD_FACTOR }, _>(num_vars, &mut rn);

        let mut transcript = KeccakTranscript::default();
        let res =
            MLSumcheck::verify_as_subprotocol(&mut transcript, num_vars, poly_degree, sum, &proof);
        assert!(res.is_ok())
    }
}

#[test]
fn verifier_rejects_proof_with_incorrect_claimed_sum() {
    let mut rng = rng();
    let num_vars = 3;

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    let one = F::one();
    let incorrect_sum = sum + one;

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        incorrect_sum,
        &proof,
    );

    assert!(matches!(
        res,
        Err(super::SumCheckError::SumCheckFailed(_, _))
    ));
}

#[test]
fn verifier_rejects_proof_with_tampered_prover_message() {
    let mut rng = rng();
    let num_vars = 3;

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    let mut tampered_proof = proof.clone();
    let one: F = F::one();
    tampered_proof.0[0].evaluations[0] += one;

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &tampered_proof,
    );

    assert!(matches!(
        res,
        Err(super::SumCheckError::SumCheckFailed(_, _))
    ));
}

#[test]
fn verifier_rejects_proof_with_wrong_degree() {
    let mut rng = rng();
    let num_vars = 3;

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    let incorrect_degree = poly_degree - 1;

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        incorrect_degree,
        sum,
        &proof,
    );

    assert!(res.is_err());
}

#[test]
fn protocol_is_deterministic_with_same_transcript() {
    let mut rng = rng();
    let num_vars = 3;

    let ((poly_mles, poly_degree), products, _) = rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let mut transcript1 = KeccakTranscript::default();
    let (proof1, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript1,
        poly_mles.clone(),
        num_vars,
        poly_degree,
        comb_fn.clone(),
    );

    let mut transcript2 = KeccakTranscript::default();
    let (proof2, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript2,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    assert_eq!(proof1, proof2);
}

#[test]
fn different_polynomials_produce_different_proofs() {
    let mut rng = rng();
    let num_vars = 3;

    let ((poly_mles1, poly_degree1), products1, _) =
        rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn1 = {
        let products = products1.clone();
        move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) }
    };

    let mut transcript1 = KeccakTranscript::default();
    let (proof1, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript1,
        poly_mles1.clone(),
        num_vars,
        poly_degree1,
        comb_fn1,
    );

    let mut poly_mles2 = poly_mles1;
    let one: F = F::one();
    poly_mles2[0].evaluations[0] += one;

    let comb_fn2 = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products1) };

    let mut transcript2 = KeccakTranscript::default();
    let (proof2, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript2,
        poly_mles2,
        num_vars,
        poly_degree1,
        comb_fn2,
    );

    assert_ne!(proof1, proof2);
}

#[test]
fn sumcheck_with_zero_polynomial() {
    let num_vars = 3;

    let poly_degree = 2;
    let num_mles = 2;
    let zero_evals = vec![F::ZERO; 1 << num_vars];
    let poly_mles: Vec<DenseMultilinearExtension<F>> = (0..num_mles)
        .map(|_| DenseMultilinearExtension::from_evaluations_vec(num_vars, zero_evals.clone()))
        .collect();

    let sum: F = F::ZERO;

    let comb_fn = |vals: &[F]| -> F { vals.iter().product() };

    let mut transcript = KeccakTranscript::default();
    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    assert!(MLSumcheck::extract_sum(&proof).is_zero());

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
    );

    assert!(res.is_ok());
}

#[test]
fn sumcheck_with_constant_polynomial() {
    let num_vars = 3;

    let poly_degree = 2;
    let num_mles = 2;
    let one: F = F::one();
    let const_evals = vec![one; 1 << num_vars];
    let poly_mles: Vec<DenseMultilinearExtension<F>> = (0..num_mles)
        .map(|_| DenseMultilinearExtension::from_evaluations_vec(num_vars, const_evals.clone()))
        .collect();

    let num_evals = 1 << num_vars;
    let sum = F::from(num_evals);

    let comb_fn = |vals: &[F]| -> F { vals.iter().product() };

    let mut transcript = KeccakTranscript::default();
    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
    );

    assert!(res.is_ok());
}

#[test]
fn sumcheck_with_single_variable() {
    let mut rng = rng();
    let num_vars = 1;

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
    );

    assert!(res.is_ok());
}

#[test]
fn verifier_rejects_proof_if_transcript_is_tampered() {
    let mut rng = rng();
    let num_vars = 3;

    let mut prover_transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut prover_transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    let mut clean_transcript = KeccakTranscript::default();
    let clean_res = MLSumcheck::verify_as_subprotocol(
        &mut clean_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
    );
    assert!(clean_res.is_ok());

    let mut tampered_transcript = KeccakTranscript::default();
    tampered_transcript.absorb(b"tampering the transcript");
    let tampered_res = MLSumcheck::verify_as_subprotocol(
        &mut tampered_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
    );
    assert!(tampered_res.is_err());
}

#[test]
#[should_panic(expected = "Prover is not active")]
fn prover_panics_if_round_exceeds_num_vars() {
    let num_vars = 3;

    let mut prover_state = ProverState {
        randomness: vec![F::ZERO; num_vars],
        mles: Vec::new(),
        num_vars,
        max_degree: 2,
        round: num_vars, // Set to the last valid round
    };

    let comb_fn = |_vals: &[F]| F::ZERO;

    let verifier_msg = Some(super::verifier::VerifierMsg {
        randomness: F::ZERO,
    });

    IPForMLSumcheck::prove_round(&mut prover_state, &verifier_msg, comb_fn);
}

#[test]
fn verifier_errors_on_incomplete_proof() {
    let mut rng = rng();
    let num_vars = 3;

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    let mut incomplete_proof = proof.clone();
    incomplete_proof.0.pop();

    let mut verifier_transcript = KeccakTranscript::default();

    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &incomplete_proof,
    );

    assert!(
        matches!(res, Err(super::SumCheckError::InvalidProofLength { expected, got }) if expected == num_vars && got == num_vars - 1),
        "expected IncorrectRoundCount error"
    );
}

#[test]
fn prover_handles_empty_mle_list() {
    let num_vars = 3;

    let poly_mles: Vec<DenseMultilinearExtension<F>> = Vec::new();
    let poly_degree = 0;
    let sum = F::ZERO;

    let comb_fn = |_vals: &[F]| -> F { F::ZERO };

    let mut transcript = KeccakTranscript::default();
    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
    );

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
    );

    assert!(res.is_ok());
}

#[test]
#[should_panic(expected = "Attempt to prove a constant.")]
fn prover_panics_with_zero_variables() {
    let num_vars = 0;
    let degree = 2;

    IPForMLSumcheck::<F, N>::prover_init(Vec::new(), num_vars, degree);
}

#[test]
fn verifier_errors_on_mismatched_nvars() {
    let mut rng = rng();
    let nvars_prover = 3;
    let nvars_verifier = 4;

    let (poly_degree, sum, proof) = generate_sumcheck_proof::<F, N, _>(nvars_prover, &mut rng);

    let mut transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut transcript,
        nvars_verifier, // verifier expects more rounds than the proof contains
        poly_degree,
        sum,
        &proof,
    );

    assert!(
        matches!(res, Err(super::SumCheckError::InvalidProofLength { expected, got })
            if expected == nvars_verifier && got == nvars_prover),
        "expected IncorrectRoundCount: expected {nvars_verifier}, got {res:?}"
    );
}

#[test]
fn verifier_produces_correct_subclaim() {
    let mut rng = rng();
    let nvars = 3;

    let mut prover_transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) = rand_poly(nvars, (2, 5), 7, &mut rng).unwrap();

    let original_mles = poly_mles.clone();
    let products_for_verification = products.clone();

    let comb_fn = move |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut prover_transcript,
        poly_mles,
        nvars,
        poly_degree,
        comb_fn,
    );

    let mut verifier_transcript = KeccakTranscript::default();
    let subclaim = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        nvars,
        poly_degree,
        sum,
        &proof,
    )
    .unwrap();

    let mle_evals_at_point: Vec<F> = original_mles
        .iter()
        .map(|mle| mle.evaluate(&subclaim.point).unwrap())
        .collect();

    let manual_eval = rand_poly_comb_fn(&mle_evals_at_point, &products_for_verification);

    assert_eq!(manual_eval, subclaim.expected_evaluation);
}

#[test]
fn zero_variable_case_returns_correct_subclaim() {
    let num_vars = 0;
    let degree = 2;

    // No prover rounds for zero-variable case
    let proof = SumcheckProof::<F>(Vec::new());

    // Let's pick some arbitrary "claimed sum"
    let claimed_sum: F = F::from(42i32);

    let mut transcript = KeccakTranscript::default();
    let subclaim =
        MLSumcheck::verify_as_subprotocol(&mut transcript, num_vars, degree, claimed_sum, &proof)
            .expect("zero-variable verification should succeed");

    // Point should be empty, and expected evaluation should match claimed_sum
    assert!(
        subclaim.point.is_empty(),
        "point should be empty for nvars=0"
    );
    assert_eq!(
        subclaim.expected_evaluation, claimed_sum,
        "expected evaluation should equal claimed sum"
    );
}
