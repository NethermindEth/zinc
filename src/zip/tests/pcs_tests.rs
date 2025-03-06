use std::str::FromStr;

use crate::{
    biginteger::BigInt,
    field_config::FieldConfig,
    poly_z::mle::DenseMultilinearExtension,
    zip::{code::ZipSpec1, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};
use ark_ff::UniformRand;
const N: usize = 2;
#[test]
fn test_zip_commitment() {
    let rng = ark_std::test_rng();
    type S = ZipSpec1;

    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng);

    let evaluations = [0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64];
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let res = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_commitment() {
    let rng = ark_std::test_rng();
    type S = ZipSpec1;

    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng);

    let evaluations = [
        0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64, 0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64,
        7i64,
    ];
    let n = 4;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

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
    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng);

    let evaluations = [0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64];
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let comm = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle).unwrap();

    let point = vec![0i64, 0i64, 0i64];

    let res = MultilinearZip::<N, S>::open(&param, &mle, &comm, &point, config, &mut transcript);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_evaluation() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = ZipSpec1;

    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng);

    let evaluations = [0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64];
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let comm = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle).unwrap();

    let point = vec![0i64, 0i64, 0i64];
    let eval = 7i64;

    let mut transcript = PcsTranscript::new();
    let _ = MultilinearZip::<N, S>::open(&param, &mle, &comm, &point, config, &mut transcript);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);

    let res = MultilinearZip::<N, S>::verify(&param, &comm, &point, &eval, &mut transcript, config);

    assert!(res.is_err())
}

#[test]
fn test_zip_evaluation() {
    let config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let rng = ark_std::test_rng();
    type S = ZipSpec1;

    let param: MultilinearZip<N, S>::Param = MultilinearZip::<N, S>::setup(8, rng);

    let evaluations = [0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64];
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let comm = MultilinearZip::<N, ZipSpec1>::commit(&param, &mle).unwrap();

    let point = vec![0i64, 0i64, 0i64];
    let eval = 0i64;

    let mut transcript = PcsTranscript::new();
    let _ = MultilinearZip::<N, S>::open(&param, &mle, &comm, &point, config, &mut transcript);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);

    let res = MultilinearZip::<N, S>::verify(&param, &comm, &point, &eval, &mut transcript, config);

    assert!(res.is_ok())
}
