use std::str::FromStr;

use crate::{
    biginteger::BigInt,
    brakedown::code::BrakedownSpec1,
    ccs::test_utils::get_dummy_ccs_from_z_length,
    field_config::FieldConfig,
    spartan::{
        structs::{ZincProver, ZincVerifier},
        SpartanProver, SpartanVerifier,
    },
    transcript::KeccakTranscript,
};
#[test]
fn test_spartan_prover() {
    const N: usize = 3;
    let config: *const FieldConfig<N> = &FieldConfig::new(
        BigInt::<3>::from_str("312829638388039969874974628075306023441").unwrap(),
    );
    let mut rng = ark_std::test_rng();
    let n = 1 << 13;
    let (_, ccs, statement, wit) = get_dummy_ccs_from_z_length::<N>(n, &mut rng, config);
    let mut transcript = KeccakTranscript::new();

    let prover = ZincProver {
        config,
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };

    let proof = prover.prove(&statement, &wit, &mut transcript, &ccs);

    assert!(proof.is_ok())
}

#[test]
fn test_spartan_verifier() {
    const N: usize = 3;
    let config: *const FieldConfig<N> = &FieldConfig::new(
        BigInt::<3>::from_str("312829638388039969874974628075306023441").unwrap(),
    );
    let mut rng = ark_std::test_rng();
    let n = 1 << 13;
    let (_, ccs, statement, wit) = get_dummy_ccs_from_z_length::<N>(n, &mut rng, config);
    let mut prover_transcript = KeccakTranscript::new();

    let prover = ZincProver {
        config,
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };

    let proof = prover
        .prove(&statement, &wit, &mut prover_transcript, &ccs)
        .unwrap();

    let verifier = ZincVerifier {
        config,
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };
    let mut verifier_transcript = KeccakTranscript::new();

    let res = verifier.verify(&statement, &proof, &mut verifier_transcript, &ccs);

    assert!(res.is_ok())
}
