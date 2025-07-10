use ark_ff::UniformRand;
use ark_std::{str::FromStr, vec, vec::Vec};
use crypto_bigint::Int;

use crate::{
    biginteger::BigInt,
    field_config::{ConfigRef, FieldConfig},
    poly_z::mle::DenseMultilinearExtension,
    traits::{ConfigReference, FieldMap},
    transcript::KeccakTranscript,
    zip::{code::ZipSpec1, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};

const I: usize = 1;
const N: usize = 2;

type TestZip = MultilinearZip<I, { 2 * I }, { 4 * I }, { 8 * I }, ZipSpec1, KeccakTranscript>;

#[test]
fn test_zip_commitment() {
    let mut transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(8, &mut transcript);

    let evaluations: Vec<_> = (0..8).map(Int::<I>::from_i32).collect();

    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let res = TestZip::commit::<N>(&param, &mle);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_commitment() {
    let mut transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(8, &mut transcript);

    let evaluations: Vec<_> = (0..16).map(Int::<I>::from_i32).collect();
    let n = 4;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let res = TestZip::commit::<N>(&param, &mle);

    assert!(res.is_err())
}

#[test]
fn test_zip_opening() {
    let config = FieldConfig::<N>::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let config = ConfigRef::from(&config);

    let mut keccak_transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(8, &mut keccak_transcript);

    let mut transcript = PcsTranscript::new();

    let evaluations: Vec<_> = (0..8).map(Int::<I>::from_i32).collect();
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let (data, _) = TestZip::commit::<N>(&param, &mle).unwrap();

    let point = vec![0i64, 0i64, 0i64].map_to_field(config);

    let res = TestZip::open(&param, &mle, &data, &point, config, &mut transcript);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_evaluation() {
    let config = FieldConfig::<N>::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let config = ConfigRef::from(&config);

    let mut keccak_transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(8, &mut keccak_transcript);

    let evaluations: Vec<_> = (0..8).map(Int::<I>::from_i32).collect();
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let (data, comm) = TestZip::commit::<N>(&param, &mle).unwrap();

    let point = vec![0i64, 0i64, 0i64].map_to_field(config);
    let eval = 7i64.map_to_field(config);

    let mut transcript = PcsTranscript::new();
    let _ = TestZip::open(&param, &mle, &data, &point, config, &mut transcript);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);
    config.reference().expect("Field config cannot be none");
    let res = TestZip::verify(&param, &comm, &point, eval, &mut transcript, config);

    assert!(res.is_err())
}

#[test]
fn test_zip_evaluation() {
    let config = FieldConfig::<N>::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let config = ConfigRef::from(&config);
    let mut rng = ark_std::test_rng();

    let n = 8;
    let mut keccak_transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(1 << n, &mut keccak_transcript);
    let evaluations: Vec<_> = (0..(1 << n))
        .map(|_| Int::<I>::from_i8(i8::rand(&mut rng)))
        .collect();
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let (data, comm) = TestZip::commit::<N>(&param, &mle).unwrap();

    let point: Vec<_> = (0..n).map(|_| Int::<I>::from(i8::rand(&mut rng))).collect();
    let eval = mle.evaluate(&point).unwrap().map_to_field(config);

    let point = point.map_to_field(config);
    let mut transcript = PcsTranscript::new();
    let _ = TestZip::open(&param, &mle, &data, &point, config, &mut transcript);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);
    config.reference().expect("Field config cannot be none");
    TestZip::verify(&param, &comm, &point, eval, &mut transcript, config)
        .expect("Failed to verify");
}
#[test]
fn test_zip_batch_evaluation() {
    let config = FieldConfig::<N>::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let config = ConfigRef::from(&config);
    let mut rng = ark_std::test_rng();

    let n = 8;
    // the number of polynomials we will batch verify;
    let m = 10;
    let mut keccak_transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(1 << n, &mut keccak_transcript);
    let evaluations: Vec<Vec<Int<I>>> = (0..m)
        .map(|_| {
            (0..(1 << n))
                .map(|_| Int::<I>::from_i8(i8::rand(&mut rng)))
                .collect::<Vec<Int<I>>>()
        })
        .collect();

    let mles: Vec<_> = evaluations
        .iter()
        .map(|evaluations| DenseMultilinearExtension::from_evaluations_slice(n, evaluations))
        .collect();

    let commitments: Vec<_> = TestZip::batch_commit::<N>(&param, &mles).unwrap();
    let (data, commitments): (Vec<_>, Vec<_>) = commitments.into_iter().unzip();
    let point: Vec<_> = (0..n).map(|_| Int::<I>::from(i8::rand(&mut rng))).collect();
    let eval: Vec<_> = mles
        .iter()
        .map(|mle| mle.evaluate(&point).unwrap().map_to_field(config))
        .collect();

    let point = point.map_to_field(config);
    let points: Vec<_> = (0..m).map(|_| point.clone()).collect();
    let mut transcript = PcsTranscript::new();
    let _ = TestZip::batch_open(&param, &mles, &data, &points, &mut transcript, config);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);
    config.reference().expect("Field config cannot be none");
    TestZip::batch_verify_z(
        &param,
        &commitments,
        &points,
        &eval,
        &mut transcript,
        config,
    )
    .expect("Failed to verify");
}
