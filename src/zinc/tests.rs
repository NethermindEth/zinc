use ark_std::marker::PhantomData;
use crypto_bigint::Zero;

use crate::{
    big_int,
    ccs::{ccs_z::get_test_ccs_stuff_Z, test_utils::get_dummy_ccs_Z_from_z_length},
    field::{ConfigRef, Int, RandomField},
    field_config,
    traits::ConfigReference,
    transcript::KeccakTranscript,
    zinc::{
        prover::SpartanProver,
        structs::{ZincProver, ZincVerifier},
        verifier::SpartanVerifier,
    },
    zip::code::ZipLinearCodeSpec1,
};

#[test]
fn test_dummy_spartan_prover() {
    const I: usize = 1;
    const N: usize = 3;
    let n = 1 << 13;
    let mut rng = ark_std::test_rng();
    let config = field_config!(312829638388039969874974628075306023441, N);
    let config = ConfigRef::from(&config);

    let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<Int<I>, RandomField<N>, ZipLinearCodeSpec1> { data: PhantomData };

    let (z_ccs, z_mle, ccs_f, statement_f) =
        ZincProver::<Int<I>, RandomField<N>, ZipLinearCodeSpec1>::prepare_for_random_field_piop(
            &statement, &wit, &ccs, config,
        )
        .expect("Failed to prepare for random field PIOP");

    let proof = SpartanProver::<Int<I>, RandomField<N>>::prove(
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
    const I: usize = 1;
    const N: usize = 3;
    let input = 3;
    let config = field_config!(312829638388039969874974628075306023441, N);
    let config = ConfigRef::from(&config);

    let (ccs, statement, wit, _) = get_test_ccs_stuff_Z(input);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<Int<I>, RandomField<N>, ZipLinearCodeSpec1> { data: PhantomData };

    let (z_ccs, z_mle, ccs_f, statement_f) =
        ZincProver::<Int<I>, RandomField<N>, ZipLinearCodeSpec1>::prepare_for_random_field_piop(
            &statement, &wit, &ccs, config,
        )
        .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<Int<I>, RandomField<N>>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
        config,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier::<Int<I>, RandomField<N>, ZipLinearCodeSpec1> { data: PhantomData };
    let mut verifier_transcript = KeccakTranscript::new();

    config.reference().expect("Field config cannot be none");

    let res = SpartanVerifier::<RandomField<N>>::verify(
        &verifier,
        &spartan_proof,
        &ccs_f,
        &mut verifier_transcript,
        config,
    );

    assert!(res.is_ok())
}

#[test]
fn test_dummy_spartan_verifier() {
    const I: usize = 1;
    const N: usize = 3;
    let n = 1 << 13;
    let mut rng = ark_std::test_rng();
    let config = field_config!(312829638388039969874974628075306023441, N);
    let config = ConfigRef::from(&config);
    let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<Int<I>, RandomField<N>, ZipLinearCodeSpec1> { data: PhantomData };

    let (z_ccs, z_mle, ccs_f, statement_f) =
        ZincProver::<Int<I>, RandomField<N>, ZipLinearCodeSpec1>::prepare_for_random_field_piop(
            &statement, &wit, &ccs, config,
        )
        .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<Int<I>, RandomField<N>>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
        config,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier::<Int<I>, RandomField<N>, ZipLinearCodeSpec1> { data: PhantomData };
    let mut verifier_transcript = KeccakTranscript::new();
    config.reference().expect("Field config cannot be none");
    let res = SpartanVerifier::<RandomField<N>>::verify(
        &verifier,
        &spartan_proof,
        &ccs_f,
        &mut verifier_transcript,
        config,
    );

    assert!(res.is_ok())
}

#[test]
fn test_failing_spartan_verifier() {
    const I: usize = 1;
    const N: usize = 3;
    let input = 3;
    let config = field_config!(312829638388039969874974628075306023441, N);
    let config = ConfigRef::from(&config);
    let (ccs, statement, mut wit, _) = get_test_ccs_stuff_Z(input);
    // Change the witness such that it is no longer valid
    assert!(wit.w_ccs[3] != Int::<I>::zero());
    wit.w_ccs[3] = Int::<I>::zero();
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver::<Int<I>, RandomField<N>, ZipLinearCodeSpec1> { data: PhantomData };

    let (z_ccs, z_mle, ccs_f, statement_f) =
        ZincProver::<Int<I>, RandomField<N>, ZipLinearCodeSpec1>::prepare_for_random_field_piop(
            &statement, &wit, &ccs, config,
        )
        .expect("Failed to prepare for random field PIOP");

    let (spartan_proof, _) = SpartanProver::<Int<I>, RandomField<N>>::prove(
        &prover,
        &statement_f,
        &z_ccs,
        &z_mle,
        &ccs_f,
        &mut prover_transcript,
        config,
    )
    .expect("Failed to generate Spartan proof");

    let verifier = ZincVerifier::<Int<I>, RandomField<N>, ZipLinearCodeSpec1> { data: PhantomData };
    let mut verifier_transcript = KeccakTranscript::new();

    config.reference().expect("Field config cannot be none");
    let res = SpartanVerifier::<RandomField<N>>::verify(
        &verifier,
        &spartan_proof,
        &ccs_f,
        &mut verifier_transcript,
        config,
    );

    assert!(res.is_err())
}
