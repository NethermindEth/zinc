use ark_std::boxed::Box;
use ark_std::rand;
use ark_std::vec;
use ark_std::vec::Vec;
use num_traits::Zero;
use rand::Rng;

use super::{
    utils::{rand_poly, rand_poly_comb_fn},
    IPForMLSumcheck, MLSumcheck, SumcheckProof,
};
use crate::poly_f::mle::DenseMultilinearExtension;
use crate::sumcheck::prover::ProverState;
use crate::traits::FieldMap;
use crate::{
    big_int,
    field::{ConfigRef, RandomField},
    field_config,
    traits::{ConfigReference, Field},
    transcript::KeccakTranscript,
};

const N: usize = 2;
type F<'cfg> = RandomField<'cfg, N>;

fn get_config() -> ConfigRef<'static, 2> {
    let config: &'static _ = Box::leak(Box::new(field_config!(57316695564490278656402085503, N)));
    ConfigRef::from(config)
}

fn generate_sumcheck_proof<F: Field>(
    num_vars: usize,
    mut rng: &mut (impl Rng + Sized),
    config: F::R,
) -> (usize, F, SumcheckProof<F>) {
    let mut transcript = KeccakTranscript::default();

    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, config, &mut rng).unwrap();

    let comb_fn = |vals: &[F]| -> F { rand_poly_comb_fn(vals, &products, config) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config,
    );
    (poly_degree, sum, proof)
}
#[test]
fn full_sumcheck_protocol_works_correctly() {
    let mut rng = ark_std::test_rng();
    let num_vars = 3;
    let config_ref = get_config();

    config_ref.reference().expect("FieldConfig cannot be null");
    for _ in 0..20 {
        let (poly_degree, sum, proof) =
            generate_sumcheck_proof::<F>(num_vars, &mut rng, config_ref);

        let mut transcript = KeccakTranscript::default();
        let res = MLSumcheck::verify_as_subprotocol(
            &mut transcript,
            num_vars,
            poly_degree,
            sum,
            &proof,
            config_ref,
        );
        assert!(res.is_ok())
    }
}

#[test]
fn verifier_rejects_proof_with_incorrect_claimed_sum() {
    let mut rng = ark_std::test_rng();
    let num_vars = 3;

    let config_ref = get_config();

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let comb_fn =
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let one: F = 1i32.map_to_field(config_ref);
    let incorrect_sum = sum + one;

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        incorrect_sum,
        &proof,
        config_ref,
    );

    assert!(matches!(
        res,
        Err(super::SumCheckError::SumCheckFailed(_, _))
    ));
}

#[test]
fn verifier_rejects_proof_with_tampered_prover_message() {
    let mut rng = ark_std::test_rng();
    let num_vars = 3;

    let config_ref = get_config();

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let comb_fn =
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let mut tampered_proof = proof.clone();
    let one: F = 1i32.map_to_field(config_ref);
    tampered_proof.0[0].evaluations[0] += one;

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &tampered_proof,
        config_ref,
    );

    assert!(matches!(
        res,
        Err(super::SumCheckError::SumCheckFailed(_, _))
    ));
}

#[test]
fn verifier_rejects_proof_with_wrong_degree() {
    let mut rng = ark_std::test_rng();
    let num_vars = 3;

    let config_ref = get_config();

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let comb_fn =
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let incorrect_degree = poly_degree - 1;

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        incorrect_degree,
        sum,
        &proof,
        config_ref,
    );

    assert!(res.is_err());
}

#[test]
fn protocol_is_deterministic_with_same_transcript() {
    let mut rng = ark_std::test_rng();
    let num_vars = 3;

    let config_ref = get_config();

    let ((poly_mles, poly_degree), products, _) =
        rand_poly(num_vars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let comb_fn =
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) };

    let mut transcript1 = KeccakTranscript::default();
    let (proof1, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript1,
        poly_mles.clone(),
        num_vars,
        poly_degree,
        comb_fn.clone(),
        config_ref,
    );

    let mut transcript2 = KeccakTranscript::default();
    let (proof2, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript2,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    assert_eq!(proof1, proof2);
}

#[test]
fn different_polynomials_produce_different_proofs() {
    let mut rng = ark_std::test_rng();
    let num_vars = 3;

    let config_ref = get_config();

    let ((poly_mles1, poly_degree1), products1, _) =
        rand_poly(num_vars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let comb_fn1 = {
        let products = products1.clone();
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) }
    };

    let mut transcript1 = KeccakTranscript::default();
    let (proof1, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript1,
        poly_mles1.clone(),
        num_vars,
        poly_degree1,
        comb_fn1,
        config_ref,
    );

    let mut poly_mles2 = poly_mles1;
    let one: F = 1i32.map_to_field(config_ref);
    poly_mles2[0].evaluations[0] += one;

    let comb_fn2 = move |vals: &[F<'static>]| -> F<'static> {
        rand_poly_comb_fn(vals, &products1, config_ref)
    };

    let mut transcript2 = KeccakTranscript::default();
    let (proof2, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript2,
        poly_mles2,
        num_vars,
        poly_degree1,
        comb_fn2,
        config_ref,
    );

    assert_ne!(proof1, proof2);
}

#[test]
fn sumcheck_with_zero_polynomial() {
    let num_vars = 3;

    let config_ref = get_config();

    let poly_degree = 2;
    let num_mles = 2;
    let zero_evals = vec![0i32; 1 << num_vars].map_to_field(config_ref);
    let poly_mles: Vec<DenseMultilinearExtension<F>> = (0..num_mles)
        .map(|_| {
            DenseMultilinearExtension::from_evaluations_vec(
                num_vars,
                zero_evals.clone(),
                config_ref,
            )
        })
        .collect();

    let sum: F = 0i32.map_to_field(config_ref);

    let comb_fn = |vals: &[F<'static>]| -> F<'static> { vals.iter().product() };

    let mut transcript = KeccakTranscript::default();
    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    assert!(MLSumcheck::extract_sum(&proof).is_zero());

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
        config_ref,
    );

    assert!(res.is_ok());
}

#[test]
fn sumcheck_with_constant_polynomial() {
    let num_vars = 3;

    let config_ref = get_config();

    let poly_degree = 2;
    let num_mles = 2;
    let one: F = 1i32.map_to_field(config_ref);
    let const_evals = vec![one; 1 << num_vars];
    let poly_mles: Vec<DenseMultilinearExtension<F>> = (0..num_mles)
        .map(|_| {
            DenseMultilinearExtension::from_evaluations_vec(
                num_vars,
                const_evals.clone(),
                config_ref,
            )
        })
        .collect();

    let num_evals = 1 << num_vars;
    let sum = num_evals.map_to_field(config_ref);

    let comb_fn = |vals: &[F<'static>]| -> F<'static> { vals.iter().product() };

    let mut transcript = KeccakTranscript::default();
    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
        config_ref,
    );

    assert!(res.is_ok());
}

#[test]
fn sumcheck_with_single_variable() {
    let mut rng = ark_std::test_rng();
    let num_vars = 1;

    let config_ref = get_config();

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let comb_fn =
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
        config_ref,
    );

    assert!(res.is_ok());
}

#[test]
fn verifier_rejects_proof_if_transcript_is_tampered() {
    let mut rng = ark_std::test_rng();
    let num_vars = 3;

    let config_ref = get_config();

    let mut prover_transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let comb_fn =
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut prover_transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let mut clean_transcript = KeccakTranscript::default();
    let clean_res = MLSumcheck::verify_as_subprotocol(
        &mut clean_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
        config_ref,
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
        config_ref,
    );
    assert!(tampered_res.is_err());
}

#[test]
#[should_panic(expected = "Prover is not active")]
fn prover_panics_if_round_exceeds_num_vars() {
    let num_vars = 3;

    let config_ref = get_config();

    let mut prover_state = ProverState {
        randomness: vec![F::zero(); num_vars],
        mles: Vec::new(),
        num_vars,
        max_degree: 2,
        round: num_vars, // Set to the last valid round
    };

    let comb_fn = |_vals: &[F<'static>]| F::zero();

    let verifier_msg = Some(super::verifier::VerifierMsg {
        randomness: F::zero(),
    });

    IPForMLSumcheck::prove_round(&mut prover_state, &verifier_msg, comb_fn, config_ref);
}

#[test]
#[should_panic(expected = "insufficient rounds")]
fn verifier_panics_on_incomplete_proof() {
    let mut rng = ark_std::test_rng();
    let num_vars = 3;

    let config_ref = get_config();

    let mut transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(num_vars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let comb_fn =
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let mut incomplete_proof = proof.clone();
    incomplete_proof.0.pop(); // Remove the last prover message

    let mut verifier_transcript = KeccakTranscript::default();

    let _ = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &incomplete_proof,
        config_ref,
    );
}

#[test]
fn prover_handles_empty_mle_list() {
    let num_vars = 3;

    let config_ref = get_config();

    let poly_mles: Vec<DenseMultilinearExtension<F>> = Vec::new();
    let poly_degree = 0;
    let sum = F::zero();

    let comb_fn = |_vals: &[F<'static>]| -> F<'static> { F::zero() };

    let mut transcript = KeccakTranscript::default();
    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        num_vars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let mut verifier_transcript = KeccakTranscript::default();
    let res = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        num_vars,
        poly_degree,
        sum,
        &proof,
        config_ref,
    );

    assert!(res.is_ok());
}

#[test]
#[should_panic(expected = "Attempt to prove a constant.")]
fn prover_panics_with_zero_variables() {
    let num_vars = 0;
    let degree = 2;

    IPForMLSumcheck::<F>::prover_init(Vec::new(), num_vars, degree);
}

#[test]
#[should_panic(expected = "insufficient rounds")]
fn verifier_rejects_proof_with_mismatched_nvars() {
    let mut rng = ark_std::test_rng();
    let nvars_prover = 3;
    let nvars_verifier = 4;
    let config_ref = get_config();

    let (poly_degree, sum, proof) =
        generate_sumcheck_proof::<F>(nvars_prover, &mut rng, config_ref);

    let mut transcript = KeccakTranscript::default();
    let _ = MLSumcheck::verify_as_subprotocol(
        &mut transcript,
        nvars_verifier,
        poly_degree,
        sum,
        &proof,
        config_ref,
    );
}

#[test]
fn verifier_produces_correct_subclaim() {
    let mut rng = ark_std::test_rng();
    let nvars = 3;
    let config_ref = get_config();

    let mut prover_transcript = KeccakTranscript::default();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(nvars, (2, 5), 7, config_ref, &mut rng).unwrap();

    let original_mles = poly_mles.clone();
    let products_for_verification = products.clone();

    let comb_fn =
        move |vals: &[F<'static>]| -> F<'static> { rand_poly_comb_fn(vals, &products, config_ref) };

    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut prover_transcript,
        poly_mles,
        nvars,
        poly_degree,
        comb_fn,
        config_ref,
    );

    let mut verifier_transcript = KeccakTranscript::default();
    let subclaim = MLSumcheck::verify_as_subprotocol(
        &mut verifier_transcript,
        nvars,
        poly_degree,
        sum,
        &proof,
        config_ref,
    )
    .unwrap();

    let mle_evals_at_point: Vec<F> = original_mles
        .iter()
        .map(|mle| mle.evaluate(&subclaim.point, config_ref).unwrap())
        .collect();

    let manual_eval =
        rand_poly_comb_fn(&mle_evals_at_point, &products_for_verification, config_ref);

    assert_eq!(manual_eval, subclaim.expected_evaluation);
}
