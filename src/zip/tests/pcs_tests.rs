use std::str::FromStr;

use crate::{
    biginteger::BigInt,
    field::{rand_with_config, RandomField},
    field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension,
    zip::{code::ZipSpec1, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};
const N: usize = 2;
#[test]
fn test_zip_commitment() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = ZipSpec1;

    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng, config);

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

    let res = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_commitment() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = ZipSpec1;

    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng, config);

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

    let res = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle);

    assert!(res.is_err())
}

#[test]
fn test_zip_opening() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = ZipSpec1;
    let mut transcript = PcsTranscript::new();
    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng, config);

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

    let comm = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle).unwrap();

    let point = vec![
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
        RandomField::from_bigint(config, 0u32.into()).unwrap(),
    ];
    let eval = RandomField::from_bigint(config, 0u32.into()).unwrap();

    let res = MultilinearZip::<N, S>::open(&param, &mle, &comm, &point, &eval, &mut transcript);

    assert!(res.is_ok())
}

#[test]
fn test_zip_evaluation() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let mut rng = ark_std::test_rng();
    type S = ZipSpec1;
    let n = 8;
    let param: MultilinearZip<N, S>::Param =
        MultilinearZip::<N, S>::setup(1 << 8, &mut rng, config);

    let evaluations: Vec<RandomField<N>> = (0..(1 << n))
        .map(|_| rand_with_config(&mut rng, config))
        .collect();

    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations, config);

    let comm = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle).unwrap();

    let point: Vec<RandomField<N>> = (0..n).map(|_| rand_with_config(&mut rng, config)).collect();

    let eval = mle.evaluate(&point, config).unwrap();

    let mut transcript = PcsTranscript::new();
    let _ = MultilinearZip::<N, S>::open(&param, &mle, &comm, &point, &eval, &mut transcript);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);

    let res = MultilinearZip::<N, S>::verify(&param, &comm, &point, &eval, &mut transcript);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_evaluation() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = ZipSpec1;

    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng, config);

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

    let comm = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle).unwrap();

    let point = vec![
        RandomField::from_bigint(config, 1u32.into()).unwrap(),
        RandomField::from_bigint(config, 1u32.into()).unwrap(),
        RandomField::from_bigint(config, 1u32.into()).unwrap(),
    ];
    let eval = RandomField::from_bigint(config, 0u32.into()).unwrap();

    let mut transcript = PcsTranscript::new();
    let _ = MultilinearZip::<N, S>::open(&param, &mle, &comm, &point, &eval, &mut transcript);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);

    let res = MultilinearZip::<N, S>::verify(&param, &comm, &point, &eval, &mut transcript);

    assert!(res.is_err())
}
