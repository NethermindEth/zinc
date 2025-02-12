use std::str::FromStr;

use crate::{
    biginteger::BigInt,
    brakedown::{code::BrakedownSpec1, pcs::MultilinearBrakedown, pcs_transcript::PcsTranscript},
    field::RandomField,
    field_config::FieldConfig,
    poly::mle::DenseMultilinearExtension,
};
const N: usize = 2;
#[test]
fn test_brakedown_commitment() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = BrakedownSpec1;

    let param: MultilinearBrakedown<N, S>::Param =
        MultilinearBrakedown::<N, S>::setup(8, 2, rng, config);

    let evaluations = [
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 1u32.into()).unwrap(),
        RandomField::from_bigint(config, 2u32.into()).unwrap(),
        RandomField::from_bigint(config, 3u32.into()).unwrap(),
        RandomField::from_bigint(config, 4u32.into()).unwrap(),
        RandomField::from_bigint(config, 5u32.into()).unwrap(),
        RandomField::from_bigint(config, 6u32.into()).unwrap(),
        RandomField::from_bigint(config, 7u32.into()).unwrap(),
    ];
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations, config);

    let res = MultilinearBrakedown::<N, BrakedownSpec1>::commit(&param, &mle);

    assert!(res.is_ok())
}

#[test]
fn test_failing_brakedown_commitment() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = BrakedownSpec1;

    let param: MultilinearBrakedown<N, S>::Param =
        MultilinearBrakedown::<N, S>::setup(8, 2, rng, config);

    let evaluations = [
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 1u32.into()).unwrap(),
        RandomField::from_bigint(config, 2u32.into()).unwrap(),
        RandomField::from_bigint(config, 3u32.into()).unwrap(),
        RandomField::from_bigint(config, 4u32.into()).unwrap(),
        RandomField::from_bigint(config, 5u32.into()).unwrap(),
        RandomField::from_bigint(config, 6u32.into()).unwrap(),
        RandomField::from_bigint(config, 7u32.into()).unwrap(),
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 1u32.into()).unwrap(),
        RandomField::from_bigint(config, 2u32.into()).unwrap(),
        RandomField::from_bigint(config, 3u32.into()).unwrap(),
        RandomField::from_bigint(config, 4u32.into()).unwrap(),
        RandomField::from_bigint(config, 5u32.into()).unwrap(),
        RandomField::from_bigint(config, 6u32.into()).unwrap(),
        RandomField::from_bigint(config, 7u32.into()).unwrap(),
    ];
    let n = 4;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations, config);

    let res = MultilinearBrakedown::<N, BrakedownSpec1>::commit(&param, &mle);

    assert!(res.is_err())
}

#[test]
fn test_brakedown_opening() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = BrakedownSpec1;
    let mut transcript = PcsTranscript::new();
    let param: MultilinearBrakedown<N, S>::Param =
        MultilinearBrakedown::<N, S>::setup(8, 2, rng, config);

    let evaluations = [
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 1u32.into()).unwrap(),
        RandomField::from_bigint(config, 2u32.into()).unwrap(),
        RandomField::from_bigint(config, 3u32.into()).unwrap(),
        RandomField::from_bigint(config, 4u32.into()).unwrap(),
        RandomField::from_bigint(config, 5u32.into()).unwrap(),
        RandomField::from_bigint(config, 6u32.into()).unwrap(),
        RandomField::from_bigint(config, 7u32.into()).unwrap(),
    ];
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations, config);

    let comm = MultilinearBrakedown::<N, BrakedownSpec1>::commit(&param, &mle).unwrap();

    let point = vec![
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
    ];
    let eval = RandomField::from_bigint(config, 0u32.into()).unwrap();

    let res =
        MultilinearBrakedown::<N, S>::open(&param, &mle, &comm, &point, &eval, &mut transcript);

    assert!(res.is_ok())
}

#[test]
fn test_brakedown_evaluation() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = BrakedownSpec1;
    let mut transcript = PcsTranscript::new();
    let param: MultilinearBrakedown<N, S>::Param =
        MultilinearBrakedown::<N, S>::setup(8, 2, rng, config);

    let evaluations = [
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 1u32.into()).unwrap(),
        RandomField::from_bigint(config, 2u32.into()).unwrap(),
        RandomField::from_bigint(config, 3u32.into()).unwrap(),
        RandomField::from_bigint(config, 4u32.into()).unwrap(),
        RandomField::from_bigint(config, 5u32.into()).unwrap(),
        RandomField::from_bigint(config, 6u32.into()).unwrap(),
        RandomField::from_bigint(config, 7u32.into()).unwrap(),
    ];
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations, config);

    let comm = MultilinearBrakedown::<N, BrakedownSpec1>::commit(&param, &mle).unwrap();

    let point = vec![
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
    ];
    let eval = RandomField::from_bigint(config, 0u32.into()).unwrap();
    let res = MultilinearBrakedown::<N, S>::verify(&param, &comm, &point, &eval, &mut transcript);

    assert!(res.is_ok())
}
