#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]
use std::str::FromStr;

use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};

use zinc::biginteger::BigInt;
use zinc::sumcheck::utils::{rand_poly, rand_poly_comb_fn};
use zinc::sumcheck::MLSumcheck;
use zinc::transcript::KeccakTranscript;
use zinc::{field::RandomField, field_config::FieldConfig};

fn run_sumcheck<const N: usize>(config: &FieldConfig<N>) {
    let nvars = 20;
    let mut rng = ark_std::test_rng();
    let ((poly_mles, poly_degree), products, sum) =
        rand_poly(nvars, (2, 5), 7, config, &mut rng).unwrap();

    let comb_fn =
        |vals: &[RandomField<N>]| -> RandomField<N> { rand_poly_comb_fn(vals, &products, config) };

    let mut transcript = KeccakTranscript::new();
    let (proof, _) = MLSumcheck::prove_as_subprotocol(
        &mut transcript,
        poly_mles,
        nvars,
        poly_degree,
        comb_fn,
        config,
    );

    let mut transcript = KeccakTranscript::default();
    let _ =
        MLSumcheck::verify_as_subprotocol(&mut transcript, nvars, poly_degree, sum, &proof, config);
}

fn bench_sumcheck_1(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("312829638388039969874974628075306023441").unwrap());

    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 1"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}

fn bench_sumcheck_2(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("322403673970412282254454204757070554431").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 2"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}
fn bench_sumcheck_3(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("285377653308127403920843585487834553973").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 3"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}
fn bench_sumcheck_4(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("233164262138933757225139946152020066289").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 4"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}
fn bench_sumcheck_5(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("236731782032802149747299945609116943963").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 5"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}
fn bench_sumcheck_6(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("219683254296065967274427818235999335877").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 6"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}

fn bench_sumcheck_7(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("301671071065105344993631035132007692107").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 7"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}
fn bench_sumcheck_8(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("198610996558066700123786608463974256119").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 8"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}
fn bench_sumcheck_9(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("195037227084167124209519505146483700673").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 9"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}
fn bench_sumcheck_10(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("267885485539660112750009519228185107803").unwrap());
    group.bench_with_input(
        BenchmarkId::new("Sumcheck", "Prime 10"),
        &config,
        |b, config| {
            b.iter(|| {
                run_sumcheck::<3>(config);
            });
        },
    );
}
fn bench_sumchecks(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    bench_sumcheck_1(group);
    bench_sumcheck_2(group);
    bench_sumcheck_3(group);
    bench_sumcheck_4(group);
    bench_sumcheck_5(group);
    bench_sumcheck_6(group);
    bench_sumcheck_7(group);
    bench_sumcheck_8(group);
    bench_sumcheck_9(group);
    bench_sumcheck_10(group);
}
pub fn sumcheck_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);

    let mut group = c.benchmark_group("Sumcheck Benchmarks");
    group.sample_size(10);
    group.plot_config(plot_config);

    bench_sumchecks(&mut group);
    group.finish();
}

criterion_group!(benches, sumcheck_benchmarks);
criterion_main!(benches);
