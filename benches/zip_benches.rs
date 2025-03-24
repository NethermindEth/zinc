#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion};
use std::str::FromStr;
use std::time::{Duration, Instant};
use zinc::transcript::KeccakTranscript;

use ark_std::test_rng;
use criterion::measurement::WallTime;
use zinc::biginteger::BigInt;

use zinc::field_config::FieldConfig;
use zinc::poly_z::mle::{DenseMultilinearExtension, MultilinearExtension};
use zinc::zip::code::{ZipSpec, ZipSpec1};
use zinc::zip::pcs::structs::MultilinearZip;
use zinc::zip::pcs_transcript::PcsTranscript;

fn commit<const N: usize, B: ZipSpec, const P: usize>(
    group: &mut BenchmarkGroup<WallTime>,
    modulus: &str,
    spec: usize,
) {
    let mut rng = test_rng();
    type T = KeccakTranscript;
    let mut keccak_transcript = T::new();
    let params = MultilinearZip::<N, B, T>::setup(1 << P, &mut keccak_transcript);

    group.bench_function(
        format!("Commit: RandomField<{N}>, poly_size = {P}, ZipSpec{spec}, modululus={modulus}"),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let poly = DenseMultilinearExtension::rand(P, &mut rng);
                    let timer = Instant::now();
                    let _ = MultilinearZip::<N, B, T>::commit(&params, &poly).unwrap();
                    total_duration += timer.elapsed()
                }

                total_duration / iters as u32
            })
        },
    );
}

fn open<const N: usize, B: ZipSpec, const P: usize>(
    group: &mut BenchmarkGroup<WallTime>,
    modulus: &str,
    spec: usize,
) {
    let mut rng = test_rng();
    let field_config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::<N>::from_str(modulus).unwrap());

    type T = KeccakTranscript;
    let mut keccak_transcript = T::new();
    let params = MultilinearZip::<N, B, T>::setup(1 << P, &mut keccak_transcript);

    let poly = DenseMultilinearExtension::rand(P, &mut rng);
    let commitment = MultilinearZip::<N, B, T>::commit(&params, &poly).unwrap();
    let point = vec![1i64; P];

    group.bench_function(
        format!("Open: RandomField<{N}>, poly_size = {P}, ZipSpec{spec}, modulus={modulus}"),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let mut transcript = PcsTranscript::new();
                    let timer = Instant::now();
                    let _ = MultilinearZip::<N, B, T>::open_z(
                        &params,
                        &poly,
                        &commitment,
                        &point,
                        field_config,
                        &mut transcript,
                    )
                    .unwrap();
                    total_duration += timer.elapsed();
                }
                total_duration / iters as u32
            })
        },
    );
}
fn verify<const N: usize, B: ZipSpec, const P: usize>(
    group: &mut BenchmarkGroup<WallTime>,
    modulus: &str,
    spec: usize,
) {
    let mut rng = test_rng();
    let field_config: *const FieldConfig<N> =
        &FieldConfig::new(BigInt::<N>::from_str(modulus).unwrap());

    type T = KeccakTranscript;
    let mut keccak_transcript = T::new();
    let params = MultilinearZip::<N, B, T>::setup(1 << P, &mut keccak_transcript);

    let poly = DenseMultilinearExtension::rand(P, &mut rng);
    let commitment = MultilinearZip::<N, B, T>::commit(&params, &poly).unwrap();
    let point = vec![1i64; P];
    let eval = poly.evaluations.last().unwrap();
    let mut transcript = PcsTranscript::new();

    let _ = MultilinearZip::<N, B, T>::open_z(
        &params,
        &poly,
        &commitment,
        &point,
        field_config,
        &mut transcript,
    )
    .unwrap();

    group.bench_function(
        format!("Verify: RandomField<{N}>, poly_size = {P}, ZipSpec{spec}, modulus={modulus}"),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let mut transcript = PcsTranscript::new();
                    let timer = Instant::now();
                    let _ = MultilinearZip::<N, B, T>::verify_z(
                        &params,
                        &commitment,
                        &point,
                        eval,
                        &mut transcript,
                        field_config,
                    );
                    total_duration += timer.elapsed();
                }
                total_duration / iters as u32
            })
        },
    );
}

fn zip_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("Zip");

    commit::<4, ZipSpec1, 12>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );
    commit::<4, ZipSpec1, 16>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );

    open::<4, ZipSpec1, 12>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );
    open::<4, ZipSpec1, 16>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );

    verify::<4, ZipSpec1, 12>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );
    verify::<4, ZipSpec1, 16>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );

    group.finish();
}

criterion_group!(benches, zip_benchmarks);
criterion_main!(benches);
