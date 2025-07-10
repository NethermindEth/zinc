#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use std::{
    str::FromStr,
    time::{Duration, Instant},
};

use ark_std::test_rng;
use criterion::{
    criterion_group, criterion_main, measurement::WallTime, BenchmarkGroup, Criterion,
};
use zinc::{
    biginteger::BigInt,
    field_config::{ConfigRef, FieldConfig},
    poly_z::mle::{DenseMultilinearExtension, MultilinearExtension},
    traits::FieldMap,
    transcript::KeccakTranscript,
    zip::{code::ZipSpec1, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};
const INT_LIMBS: usize = 1;
const FIELD_LIMBS: usize = 4;
type BenchZip = MultilinearZip<
    INT_LIMBS,
    { 2 * INT_LIMBS },
    { 4 * INT_LIMBS },
    { 8 * INT_LIMBS },
    ZipSpec1,
    KeccakTranscript,
>;

fn commit<const P: usize>(group: &mut BenchmarkGroup<WallTime>, modulus: &str, spec: usize) {
    let mut rng = test_rng();
    type T = KeccakTranscript;
    let mut keccak_transcript = T::new();
    let params = BenchZip::setup(1 << P, &mut keccak_transcript);

    group.bench_function(
        format!("Commit: RandomField<{FIELD_LIMBS}>, poly_size = 2^{P}(Int limbs = {INT_LIMBS}), ZipSpec{spec}, modulus={modulus}"),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let poly = DenseMultilinearExtension::rand(P, &mut rng);
                    let timer = Instant::now();
                    let _ = BenchZip::commit::<FIELD_LIMBS>(&params, &poly).expect("Failed to commit");
                    total_duration += timer.elapsed()
                }

                total_duration / iters as u32
            })
        },
    );
}

fn open<const P: usize>(group: &mut BenchmarkGroup<WallTime>, modulus: &str, spec: usize) {
    let mut rng = test_rng();
    let config = FieldConfig::new(BigInt::<FIELD_LIMBS>::from_str(modulus).unwrap());
    let field_config = ConfigRef::from(&config);

    type T = KeccakTranscript;
    let mut keccak_transcript = T::new();
    let params = BenchZip::setup(1 << P, &mut keccak_transcript);

    let poly = DenseMultilinearExtension::rand(P, &mut rng);
    let (data, _) = BenchZip::commit::<FIELD_LIMBS>(&params, &poly).unwrap();
    let point = vec![1i64; P];

    group.bench_function(
        format!("Open: RandomField<{FIELD_LIMBS}>, poly_size = 2^{P}(Int limbs = {INT_LIMBS}), ZipSpec{spec}, modulus={modulus}"),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let mut transcript = PcsTranscript::new();
                    let timer = Instant::now();
                    BenchZip::open(
                        &params,
                        &poly,
                        &data,
                        &point.map_to_field(field_config),
                        field_config,
                        &mut transcript,
                    )
                    .expect("Failed to make opening");
                    total_duration += timer.elapsed();
                }
                total_duration / iters as u32
            })
        },
    );
}
fn verify<const P: usize>(group: &mut BenchmarkGroup<WallTime>, modulus: &str, spec: usize) {
    let mut rng = test_rng();
    let config = FieldConfig::new(BigInt::<FIELD_LIMBS>::from_str(modulus).unwrap());
    let field_config = ConfigRef::from(&config);

    type T = KeccakTranscript;
    let mut keccak_transcript = T::new();
    let params = BenchZip::setup(1 << P, &mut keccak_transcript);

    let poly = DenseMultilinearExtension::rand(P, &mut rng);
    let (data, commitment) = BenchZip::commit::<FIELD_LIMBS>(&params, &poly).unwrap();
    let point = vec![1i64; P];
    let eval = poly.evaluations.last().unwrap();
    let mut transcript = PcsTranscript::new();

    BenchZip::open(
        &params,
        &poly,
        &data,
        &point.map_to_field(field_config),
        field_config,
        &mut transcript,
    )
    .unwrap();

    let proof = transcript.into_proof();
    field_config
        .reference()
        .expect("Field config cannot be none");
    group.bench_function(
        format!("Verify: RandomField<{FIELD_LIMBS}>, poly_size = 2^{P}(Int limbs = {INT_LIMBS}), ZipSpec{spec}, modulus={modulus}"),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let mut transcript = PcsTranscript::from_proof(&proof);
                    let timer = Instant::now();
                    BenchZip::verify(
                        &params,
                        &commitment,
                        &point.map_to_field(field_config),
                        eval.map_to_field(field_config),
                        &mut transcript,
                        field_config,
                    )
                    .expect("Failed to verify");
                    total_duration += timer.elapsed();
                }
                total_duration / iters as u32
            })
        },
    );
}

fn zip_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("Zip");

    commit::<12>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );
    commit::<16>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );

    open::<12>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );
    open::<16>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );

    verify::<12>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );
    verify::<16>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );

    group.finish();
}

criterion_group!(benches, zip_benchmarks);
criterion_main!(benches);
