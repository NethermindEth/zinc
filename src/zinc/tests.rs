use std::str::FromStr;

use crate::{
    biginteger::BigInt,
    brakedown::code::BrakedownSpec1,
    ccs::{ccs_f::get_test_ccs_stuff_F, test_utils::get_dummy_ccs_F_from_z_length},
    field::{conversion::FieldMap, RandomField},
    field_config::FieldConfig,
    transcript::KeccakTranscript,
    zinc::{
        prover::SpartanProver,
        structs::{ZincProver, ZincVerifier},
        utils::draw_random_field,
        verifier::SpartanVerifier,
    },
};
#[test]
fn test_spartan_prover() {
    let n = 1 << 13;
    let mut rng = ark_std::test_rng();
    let config =
        FieldConfig::new(BigInt::<3>::from_str("312829638388039969874974628075306023441").unwrap());
    let (_, ccs, statement, wit) = get_dummy_ccs_F_from_z_length(n, &mut rng, &config);
    let mut transcript = KeccakTranscript::new();

    let prover = ZincProver::<3, _> {
        // If we are keeping primes around 128 bits we should stay with N = 3 hardcoded
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };

    let proof = SpartanProver::<3>::prove(
        &prover,
        &statement,
        &wit,
        &mut transcript,
        &ccs,
        &config,
    );

    assert!(proof.is_ok())
}

#[test]
fn test_spartan_verifier() {
    let input = 3;
    let config =
        FieldConfig::new(BigInt::<3>::from_str("312829638388039969874974628075306023441").unwrap());
    let (ccs, statement, wit, _) = get_test_ccs_stuff_F(input, &config);

    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver {
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };

    let proof = SpartanProver::<3>::prove(
        &prover,
        &statement,
        &wit,
        &mut prover_transcript,
        &ccs,
        &config,
    )
    .unwrap();

    let verifier = ZincVerifier::<3, _> {
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };
    let mut verifier_transcript = KeccakTranscript::new();

    let res = verifier.verify(&statement, proof, &mut verifier_transcript, &ccs);

    assert!(res.is_ok())
}

#[test]
fn test_failing_spartan_verifier() {
    const N: usize = 3;
    let config: *const FieldConfig<N> = &FieldConfig::new(
        BigInt::<3>::from_str("312829638388039969874974628075306023441").unwrap(),
    );
    let input = 3;
    let (ccs, statement, mut wit, _) = get_test_ccs_stuff_F(input, config);
    // Change the witness such that it is no longer valid
    wit.w_ccs[3] = RandomField::from_bigint(config, 0u32.into()).unwrap();
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver {
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };

    let proof = SpartanProver::<N>::prove(
        &prover,
        &statement,
        &wit,
        &mut prover_transcript,
        &ccs,
        config,
    )
    .unwrap();

    let verifier = ZincVerifier {
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };
    let mut verifier_transcript = KeccakTranscript::new();

    let res = verifier.verify(&statement, proof, &mut verifier_transcript, &ccs);

    assert!(res.is_err())
}
