#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion};
use std::str::FromStr;
use std::time::{Duration, Instant};

use ark_std::test_rng;
use criterion::measurement::WallTime;
use zinc::biginteger::BigInt;
use zinc::brakedown::code::{BrakedownSpec, BrakedownSpec1};
use zinc::brakedown::pcs::MultilinearBrakedown;
use zinc::brakedown::pcs_transcript::PcsTranscript;
use zinc::field::RandomField;
use zinc::field_config::FieldConfig;
use zinc::poly::mle::{DenseMultilinearExtension, MultilinearExtension};

fn commit<const N: usize, B: BrakedownSpec, const P: usize>(
    group: &mut BenchmarkGroup<WallTime>,
    modulus: &str,
    spec: usize,
) {
    let mut rng = test_rng();
    let field_config = FieldConfig::new(BigInt::<N>::from_str(modulus).unwrap());

    let params = MultilinearBrakedown::<N, B>::setup(P, 0, &mut rng);

    group.bench_function(
        format!(
            "Commit: RandomField<{N}>, poly_size = {P}, BrakedownSpec{spec}, modululus={modulus}"
        ),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let poly = DenseMultilinearExtension::<N>::rand(8, &field_config, &mut rng);
                    let timer = Instant::now();
                    let _ = MultilinearBrakedown::<N, B>::commit(&params, &poly).unwrap();
                    total_duration += timer.elapsed()
                }

                total_duration / iters as u32
            })
        },
    );
}

fn batch_commit<const N: usize, B: BrakedownSpec, const P: usize>(
    group: &mut BenchmarkGroup<WallTime>,
    modulus: &str,
    spec: usize,
    batch_size: usize,
) {
    let mut rng = test_rng();
    let field_config = FieldConfig::new(BigInt::<N>::from_str(modulus).unwrap());

    let params = MultilinearBrakedown::<N, B>::setup(P, 0, &mut rng);

    group.bench_function(
        format!(
            "BatchCommit: RandomField<{N}>, poly_size = {P}, BrakedownSpec{spec}, batch_size = {batch_size}, modulus={modulus}"
        ),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let polys: Vec<_> = (0..batch_size)
                        .map(|_| DenseMultilinearExtension::<N>::rand(4, &field_config, &mut rng))
                        .collect();

                    let timer = Instant::now();
                    let _ = MultilinearBrakedown::<N, B>::batch_commit(&params, &polys).unwrap();
                    total_duration += timer.elapsed()
                }

                total_duration / iters as u32
            })
        }
    );
}

fn verify<const N: usize, B: BrakedownSpec, const P: usize>(
    group: &mut BenchmarkGroup<WallTime>,
    modulus: &str,
    spec: usize,
) {
    let mut rng = test_rng();
    let field_config = FieldConfig::new(BigInt::<N>::from_str(modulus).unwrap());

    let params = MultilinearBrakedown::<N, B>::setup(P, 0, &mut rng);

    let poly = DenseMultilinearExtension::<N>::rand(4, &field_config, &mut rng);
    let commitment = MultilinearBrakedown::<N, B>::commit(&params, &poly).unwrap();
    let point = vec![RandomField::<N>::from(2u32); 4];
    let eval = RandomField::<N>::from(10u32);
    let mut transcript = PcsTranscript::new();

    let _ = MultilinearBrakedown::<N, B>::open(
        &params,
        &poly,
        &commitment,
        &point,
        &eval,
        &mut transcript,
    )
    .unwrap();

    group.bench_function(
        format!(
            "Verify: RandomField<{N}>, poly_size = {P}, BrakedownSpec{spec}, modulus={modulus}"
        ),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let mut transcript = PcsTranscript::new();
                    let timer = Instant::now();
                    let _ = MultilinearBrakedown::<N, B>::verify(
                        &params,
                        &commitment,
                        &point,
                        &eval,
                        &mut transcript,
                    );
                    total_duration += timer.elapsed();
                }
                total_duration / iters as u32
            })
        },
    );
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("PCS");

    commit::<1, BrakedownSpec1, { 1 << 10 }>(&mut group, "18446744073709551557", 1);
    commit::<2, BrakedownSpec1, { 1 << 10 }>(
        &mut group,
        "292235773983524966147578932818224452679",
        1,
    );
    commit::<4, BrakedownSpec1, { 1 << 10 }>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );

    let batch_size = 8;

    batch_commit::<1, BrakedownSpec1, { 1 << 10 }>(
        &mut group,
        "18446744073709551557",
        1,
        batch_size,
    );
    batch_commit::<2, BrakedownSpec1, { 1 << 10 }>(
        &mut group,
        "292235773983524966147578932818224452679",
        1,
        batch_size,
    );
    batch_commit::<4, BrakedownSpec1, { 1 << 10 }>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
        batch_size,
    );

    verify::<1, BrakedownSpec1, { 1 << 10 }>(&mut group, "18446744073709551557", 1);
    verify::<2, BrakedownSpec1, { 1 << 10 }>(
        &mut group,
        "292235773983524966147578932818224452679",
        1,
    );
    verify::<4, BrakedownSpec1, { 1 << 10 }>(
        &mut group,
        "106319353542452952636349991594949358997917625194731877894581586278529202198383",
        1,
    );

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
