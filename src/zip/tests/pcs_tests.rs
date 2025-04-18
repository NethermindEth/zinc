use crate::{
    biginteger::BigInt,
    field::conversion::FieldMap,
    field_config::FieldConfig,
    poly_z::mle::DenseMultilinearExtension,
    transcript::KeccakTranscript,
    zip::{code::ZipSpec1, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};
use ark_ff::UniformRand;
use std::str::FromStr;

const N: usize = 2;
#[test]
fn test_zip_commitment() {
    type S = ZipSpec1;
    type T = KeccakTranscript;
    let mut transcript = KeccakTranscript::new();
    let param: MultilinearZip<N, S, T>::Param =
        MultilinearZip::<N, S, T>::setup(8, &mut transcript);

    let evaluations = [0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64];
    let n = 3;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let res = MultilinearZip::<N, S, T>::commit(&param, &mle);

    assert!(res.is_ok())
}

#[test]
fn test_failing_zip_commitment() {
    type S = ZipSpec1;
    type T = KeccakTranscript;
    let mut transcript = KeccakTranscript::new();
    let param: MultilinearZip<N, S, T>::Param =
        MultilinearZip::<N, S, T>::setup(8, &mut transcript);

    let evaluations = [
        0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64, 0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64,
        7i64,
    ];
    let n = 4;
    let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

    let res = MultilinearZip::<N, ZipSpec1, T>::commit(&param, &mle);

    assert!(res.is_err())
}

// #[test]
// fn test_zip_opening() {
//     let config: *const FieldConfig<N> =
//         &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());

//     type S = ZipSpec1;
//     type T = KeccakTranscript;
//     let mut keccak_transcript = KeccakTranscript::new();
//     let param: MultilinearZip<N, S, T>::Param =
//         MultilinearZip::<N, S, T>::setup(8, &mut keccak_transcript);

//     let mut transcript = PcsTranscript::new();

//     let evaluations = [0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64];
//     let n = 3;
//     let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

//     let (data, _) = MultilinearZip::<N, ZipSpec1, T>::commit(&param, &mle).unwrap();

//     let point = vec![0i64, 0i64, 0i64];

//     let res =
//         MultilinearZip::<N, S, T>::open_z(&param, &mle, &data, &point, config, &mut transcript);

//     assert!(res.is_ok())
// }

// #[test]
// fn test_failing_zip_evaluation() {
//     let config: *const FieldConfig<N> =
//         &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());

//     type S = ZipSpec1;
//     type T = KeccakTranscript;
//     let mut keccak_transcript = KeccakTranscript::new();
//     let param: MultilinearZip<N, S, T>::Param =
//         MultilinearZip::<N, S, T>::setup(8, &mut keccak_transcript);

//     let evaluations = [0i64, 1i64, 2i64, 3i64, 4i64, 5i64, 6i64, 7i64];
//     let n = 3;
//     let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

//     let (data, comm) = MultilinearZip::<N, ZipSpec1, T>::commit(&param, &mle).unwrap();

//     let point = vec![0i64, 0i64, 0i64];
//     let eval = 7i64;

//     let mut transcript = PcsTranscript::new();
//     let _ = MultilinearZip::<N, S, T>::open_z(&param, &mle, &data, &point, config, &mut transcript);

//     let proof = transcript.into_proof();
//     let mut transcript = PcsTranscript::from_proof(&proof);

//     let res =
//         MultilinearZip::<N, S, T>::verify_z(&param, &comm, &point, &eval, &mut transcript, config);

//     assert!(res.is_err())
// }

// #[test]
// fn test_zip_evaluation() {
//     let config: *const FieldConfig<N> =
//         &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
//     let mut rng = ark_std::test_rng();
//     type S = ZipSpec1;
//     type T = KeccakTranscript;
//     let n = 8;
//     let mut keccak_transcript = KeccakTranscript::new();
//     let param: MultilinearZip<N, S, T>::Param =
//         MultilinearZip::<N, S, T>::setup(1 << n, &mut keccak_transcript);
//     let evaluations: Vec<_> = (0..(1 << n))
//         .map(|_| i64::from(i8::rand(&mut rng)))
//         .collect();
//     let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

//     let (data, comm) = MultilinearZip::<N, ZipSpec1, T>::commit(&param, &mle).unwrap();

//     let point: Vec<_> = (0..n).map(|_| i64::from(i8::rand(&mut rng))).collect();
//     let eval = mle.evaluate(&point).unwrap();

//     let mut transcript = PcsTranscript::new();
//     let _ = MultilinearZip::<N, S, T>::open_z(&param, &mle, &data, &point, config, &mut transcript);

//     let proof = transcript.into_proof();
//     let mut transcript = PcsTranscript::from_proof(&proof);

//     MultilinearZip::<N, S, T>::verify_z(&param, &comm, &point, &eval, &mut transcript, config)
//         .expect("Failed to verify");
// }

// #[test]
// fn test_zip_evaluation_field() {
//     let config: *const FieldConfig<N> =
//         &FieldConfig::new(BigInt::from_str("57316695564490278656402085503").unwrap());
//     let mut rng = ark_std::test_rng();
//     type S = ZipSpec1;
//     type T = KeccakTranscript;
//     let n = 8;
//     let mut keccak_transcript = KeccakTranscript::new();
//     let param: MultilinearZip<N, S, T>::Param =
//         MultilinearZip::<N, S, T>::setup(1 << n, &mut keccak_transcript);
//     let evaluations: Vec<_> = (0..(1 << n)).map(|_| i64::rand(&mut rng)).collect();
//     let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

//     let (data, comm) = MultilinearZip::<N, ZipSpec1, T>::commit(&param, &mle).unwrap();

//     let point: Vec<_> = (0..n).map(|_| 1u32.map_to_field(config)).collect();
//     let eval = evaluations[(1 << n) - 1].map_to_field(config);

//     let mut transcript = PcsTranscript::new();
//     let _ = MultilinearZip::<N, S, T>::open_f(&param, &mle, &data, &point, config, &mut transcript);

//     let proof = transcript.into_proof();
//     let mut transcript = PcsTranscript::from_proof(&proof);

//     let res =
//         MultilinearZip::<N, S, T>::verify_f(&param, &comm, &point, &eval, &mut transcript, config);

//     assert!(res.is_ok())
// }
