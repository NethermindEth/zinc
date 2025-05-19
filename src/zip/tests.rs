use crate::field_config::ConfigPtr;
use crate::{
    biginteger::BigInt,
    field::conversion::FieldMap,
    field_config::FieldConfig,
    poly_z::mle::DenseMultilinearExtension,
    transcript::KeccakTranscript,
    zip::{code::ZipSpec1, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};
use ark_ff::UniformRand;
use crypto_bigint::Int;
use std::str::FromStr;

const N: usize = 2;

type TestZip<'cfg> =
    MultilinearZip<N, { 2 * N }, { 4 * N }, { 8 * N }, ZipSpec1, KeccakTranscript<'cfg>>;

#[test]
fn test_zip_commitment() {
    let mut transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(8, &mut transcript);

    let evaluations: Vec<_> = (0..8).map(Int::<N>::from_i32).collect();

    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let res = TestZip::commit(&param, &mle);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_commitment() {
    let mut transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(8, &mut transcript);

    let evaluations: Vec<_> = (0..16).map(Int::<N>::from_i32).collect();
    let n = 4;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let res = TestZip::commit(&param, &mle);

    assert!(res.is_err())
}

#[test]
fn test_zip_opening() {
    let config = FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let config = ConfigPtr::from(&config);

    let mut keccak_transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(8, &mut keccak_transcript);

    let mut transcript = PcsTranscript::new();

    let evaluations: Vec<_> = (0..8).map(Int::<N>::from_i32).collect();
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let (data, _) = TestZip::commit(&param, &mle).unwrap();

    let point = vec![0i64, 0i64, 0i64].map_to_field(config);

    let res = TestZip::open(&param, &mle, &data, &point, config, &mut transcript);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_evaluation() {
    let config = FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let config = ConfigPtr::from(&config);

    let mut keccak_transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(8, &mut keccak_transcript);

    let evaluations: Vec<_> = (0..8).map(Int::<N>::from_i32).collect();
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let (data, comm) = TestZip::commit(&param, &mle).unwrap();

    let point = vec![0i64, 0i64, 0i64].map_to_field(config);
    let eval = 7i64.map_to_field(config);

    let mut transcript = PcsTranscript::new();
    let _ = TestZip::open(&param, &mle, &data, &point, config, &mut transcript);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);

    let res = TestZip::verify(
        &param,
        &comm,
        &point,
        eval,
        &mut transcript,
        config.reference().expect("Field config cannot be none"),
    );

    assert!(res.is_err())
}

#[test]
fn test_zip_evaluation() {
    let config = FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let config = ConfigPtr::from(&config);
    let mut rng = ark_std::test_rng();

    let n = 8;
    let mut keccak_transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(1 << n, &mut keccak_transcript);
    let evaluations: Vec<_> = (0..(1 << n))
        .map(|_| Int::<N>::from_i8(i8::rand(&mut rng)))
        .collect();
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let (data, comm) = TestZip::commit(&param, &mle).unwrap();

    let point: Vec<_> = (0..n).map(|_| Int::<N>::from(i8::rand(&mut rng))).collect();
    let eval = mle.evaluate(&point).unwrap().map_to_field(config);

    let point = point.map_to_field(config);
    let mut transcript = PcsTranscript::new();
    let _ = TestZip::open(&param, &mle, &data, &point, config, &mut transcript);

    let proof = transcript.into_proof();
    let mut transcript = PcsTranscript::from_proof(&proof);

    TestZip::verify(
        &param,
        &comm,
        &point,
        eval,
        &mut transcript,
        config.reference().expect("Field config cannot be none"),
    )
    .expect("Failed to verify");
}
#[test]
fn test_zip_batch_evaluation() {
    let config = FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
    let config = ConfigPtr::from(&config);
    let mut rng = ark_std::test_rng();

    let n = 8;
    // the number of polynomials we will batch verify;
    let m = 10;
    let mut keccak_transcript = KeccakTranscript::new();
    let param: TestZip::Param = TestZip::setup(1 << n, &mut keccak_transcript);
    let evaluations: Vec<Vec<Int<N>>> = (0..m)
        .map(|_| {
            (0..(1 << n))
                .map(|_| Int::<N>::from_i8(i8::rand(&mut rng)))
                .collect::<Vec<Int<N>>>()
        })
        .collect();

    let mles: Vec<_> = evaluations
        .iter()
        .map(|evaluations| DenseMultilinearExtension::from_evaluations_slice(n, evaluations))
        .collect();

    let commitments: Vec<_> = TestZip::batch_commit(&param, &mles).unwrap();
    let (data, commitments): (Vec<_>, Vec<_>) = commitments.into_iter().unzip();
    let point: Vec<_> = (0..n).map(|_| Int::<N>::from(i8::rand(&mut rng))).collect();
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

    TestZip::batch_verify_z(
        &param,
        &commitments,
        &points,
        &eval,
        &mut transcript,
        config.reference().expect("Field config cannot be none"),
    )
    .expect("Failed to verify");
}
