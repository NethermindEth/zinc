use std::str::FromStr;

use crate::{
    biginteger::BigInt,
    ccs::{ccs_z::get_test_ccs_stuff_Z, test_utils::get_dummy_ccs_Z_from_z_length},
    field_config::FieldConfig,
    transcript::KeccakTranscript,
    zinc::{
        prover::SpartanProver,
        structs::{ZincProver, ZincVerifier},
        verifier::SpartanVerifier,
    },
    zip::code::ZipSpec1,
};
#[test]
fn test_dummy_spartan_prover() {
    const N: usize = 3;
    let n = 1 << 13;
    let mut rng = ark_std::test_rng();
    let config: *const FieldConfig<N> = &FieldConfig::new(
        BigInt::<N>::from_str("312829638388039969874974628075306023441").unwrap(),
    );
    let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<N, _> {
        data: std::marker::PhantomData::<ZipSpec1>,
    };

    let (z_ccs, z_mle, ccs_f, statement_f) =
        ZincProver::<N, ZipSpec1>::prepare_for_random_field_piop(&statement, &wit, &ccs, config)
            .expect("Failed to prepare for random field PIOP");

    let proof = SpartanProver::<N>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
        config,
    );

    assert!(proof.is_ok())
}

#[test]
fn test_spartan_verifier() {
    const N: usize = 3;
    let input = 3;
    let config: *const FieldConfig<N> = &FieldConfig::new(
        BigInt::<N>::from_str("312829638388039969874974628075306023441").unwrap(),
    );
    let (ccs, statement, wit, _) = get_test_ccs_stuff_Z(input);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<N, _> {
        data: std::marker::PhantomData::<ZipSpec1>,
    };

    let (z_ccs, z_mle, ccs_f, statement_f) =
        ZincProver::<N, ZipSpec1>::prepare_for_random_field_piop(&statement, &wit, &ccs, config)
            .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<N>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
        config,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier::<N, _> {
        data: std::marker::PhantomData::<ZipSpec1>,
    };
    let mut verifier_transcript = KeccakTranscript::new();

    let res =
        SpartanVerifier::<N>::verify(&verifier, &spartan_proof, &mut verifier_transcript, &ccs_f);

    assert!(res.is_ok())
}

#[test]
fn test_dummy_spartan_verifier() {
    const N: usize = 3;
    let n = 1 << 13;
    let mut rng = ark_std::test_rng();
    let config: *const FieldConfig<N> = &FieldConfig::new(
        BigInt::<N>::from_str("312829638388039969874974628075306023441").unwrap(),
    );
    let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<N, _> {
        data: std::marker::PhantomData::<ZipSpec1>,
    };

    let (z_ccs, z_mle, ccs_f, statement_f) =
        ZincProver::<N, ZipSpec1>::prepare_for_random_field_piop(&statement, &wit, &ccs, config)
            .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<N>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
        config,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier::<N, _> {
        data: std::marker::PhantomData::<ZipSpec1>,
    };
    let mut verifier_transcript = KeccakTranscript::new();

    let res =
        SpartanVerifier::<N>::verify(&verifier, &spartan_proof, &mut verifier_transcript, &ccs_f);

    assert!(res.is_ok())
}

#[test]
fn test_failing_spartan_verifier() {
    const N: usize = 3;
    let input = 3;
    let config: *const FieldConfig<N> = &FieldConfig::new(
        BigInt::<N>::from_str("312829638388039969874974628075306023441").unwrap(),
    );
    let (ccs, statement, mut wit, _) = get_test_ccs_stuff_Z(input);
    // Change the witness such that it is no longer valid
    assert!(wit.w_ccs[3] != 0);
    wit.w_ccs[3] = 0i64;
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<N, _> {
        data: std::marker::PhantomData::<ZipSpec1>,
    };

    let (z_ccs, z_mle, ccs_f, statement_f) =
        ZincProver::<N, ZipSpec1>::prepare_for_random_field_piop(&statement, &wit, &ccs, config)
            .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<N>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
        config,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier {
        data: std::marker::PhantomData::<ZipSpec1>,
    };
    let mut verifier_transcript = KeccakTranscript::new();

    let res =
        SpartanVerifier::<N>::verify(&verifier, &spartan_proof, &mut verifier_transcript, &ccs_f);

    assert!(res.is_err())
}
