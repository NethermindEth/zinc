#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion};
use std::str::FromStr;
use std::time::{Duration, Instant};
use zinc::field::conversion::FieldMap;
use zinc::transcript::KeccakTranscript;

use ark_std::test_rng;
use criterion::measurement::WallTime;
use zinc::biginteger::BigInt;

use zinc::field_config::{ConfigRef, FieldConfig};
use zinc::poly_z::mle::{DenseMultilinearExtension, MultilinearExtension};
use zinc::zip::code::ZipSpec1;
use zinc::zip::pcs::structs::MultilinearZip;
use zinc::zip::pcs_transcript::PcsTranscript;
const N: usize = 4;
type BenchZip<'cfg> =
    MultilinearZip<N, { 2 * N }, { 4 * N }, { 8 * N }, ZipSpec1, KeccakTranscript<'cfg>>;

fn commit<const P: usize>(group: &mut BenchmarkGroup<WallTime>, modulus: &str, spec: usize) {
    let mut rng = test_rng();
    type T<'cfg> = KeccakTranscript<'cfg>;
    let mut keccak_transcript = T::new();
    let params = BenchZip::setup(1 << P, &mut keccak_transcript);

    group.bench_function(
        format!("Commit: RandomField<{N}>, poly_size = 2^{P}, ZipSpec{spec}, modulus={modulus}"),
        |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::ZERO;
                for _ in 0..iters {
                    let poly = DenseMultilinearExtension::rand(P, &mut rng);
                    let timer = Instant::now();
                    let _ = BenchZip::commit(&params, &poly).expect("Failed to commit");
                    total_duration += timer.elapsed()
                }

                total_duration / iters as u32
            })
        },
    );
}

fn open<const P: usize>(group: &mut BenchmarkGroup<WallTime>, modulus: &str, spec: usize) {
    let mut rng = test_rng();
    let config = FieldConfig::new(BigInt::<N>::from_str(modulus).unwrap());
    let field_config = ConfigRef::from(&config);

    type T<'cfg> = KeccakTranscript<'cfg>;
    let mut keccak_transcript = T::new();
    let params = BenchZip::setup(1 << P, &mut keccak_transcript);

    let poly = DenseMultilinearExtension::rand(P, &mut rng);
    let (data, _) = BenchZip::commit(&params, &poly).unwrap();
    let point = vec![1i64; P];

    group.bench_function(
        format!("Open: RandomField<{N}>, poly_size = 2^{P}, ZipSpec{spec}, modulus={modulus}"),
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
    let config = FieldConfig::new(BigInt::<N>::from_str(modulus).unwrap());
    let field_config = ConfigRef::from(&config);

    type T<'cfg> = KeccakTranscript<'cfg>;
    let mut keccak_transcript = T::new();
    let params = BenchZip::setup(1 << P, &mut keccak_transcript);

    let poly = DenseMultilinearExtension::rand(P, &mut rng);
    let (data, commitment) = BenchZip::commit(&params, &poly).unwrap();
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

    group.bench_function(
        format!("Verify: RandomField<{N}>, poly_size = 2^{P}, ZipSpec{spec}, modulus={modulus}"),
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
                        field_config
                            .reference()
                            .expect("Field config cannot be none"),
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
