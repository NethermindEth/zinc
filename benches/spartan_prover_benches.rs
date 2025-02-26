use std::str::FromStr;

use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use zinc::{
    biginteger::BigInt,
    brakedown::code::BrakedownSpec1,
    ccs::test_utils::get_dummy_ccs_from_z_length,
    field_config::FieldConfig,
    spartan::{prover::SpartanProver, structs::ZincProver},
    transcript::KeccakTranscript,
};

fn run_spartan_prover<const N: usize>(n: usize, config: &FieldConfig<N>) {
    let mut rng = ark_std::test_rng();

    let (_, ccs, statement, wit) = get_dummy_ccs_from_z_length::<N>(n, &mut rng, config);
    let mut transcript = KeccakTranscript::new();

    let prover = ZincProver {
        config,
        data: std::marker::PhantomData::<BrakedownSpec1>,
    };
    let _proof = prover.prove(&statement, &wit, &mut transcript, &ccs);
}

fn bench_spartan_prover_1(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInt::<3>::from_str("312829638388039969874974628075306023441").unwrap());

    group.bench_with_input(
        BenchmarkId::new("Spartan Prover", "Prime 1, wit len 2^10"),
        &config,
        |b, config| {
            b.iter(|| {
                run_spartan_prover::<3>(1 << 10, config);
            });
        },
    );
    group.bench_with_input(
        BenchmarkId::new("Spartan Prover", "Prime 1, wit len 2^11"),
        &config,
        |b, config| {
            b.iter(|| {
                run_spartan_prover::<3>(1 << 11, config);
            });
        },
    );
    group.bench_with_input(
        BenchmarkId::new("Spartan Prover", "Prime 1, wit len 2^12"),
        &config,
        |b, config| {
            b.iter(|| {
                run_spartan_prover::<3>(1 << 12, config);
            });
        },
    );
    group.bench_with_input(
        BenchmarkId::new("Spartan Prover", "Prime 1, wit len 2^13"),
        &config,
        |b, config| {
            b.iter(|| {
                run_spartan_prover::<3>(1 << 13, config);
            });
        },
    );
    group.bench_with_input(
        BenchmarkId::new("Spartan Prover", "Prime 1, wit len 2^14"),
        &config,
        |b, config| {
            b.iter(|| {
                run_spartan_prover::<3>(1 << 14, config);
            });
        },
    );
    group.bench_with_input(
        BenchmarkId::new("Spartan Prover", "Prime 1, wit len 2^15"),
        &config,
        |b, config| {
            b.iter(|| {
                run_spartan_prover::<3>(1 << 15, config);
            });
        },
    );
    group.bench_with_input(
        BenchmarkId::new("Spartan Prover", "Prime 1, wit len 2^16"),
        &config,
        |b, config| {
            b.iter(|| {
                run_spartan_prover::<3>(1 << 16, config);
            });
        },
    );
}

pub fn bench_spartan_prover(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
) {
    bench_spartan_prover_1(group);
}

pub fn spartan_prover_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);

    let mut group = c.benchmark_group("spartan_prover");
    group.plot_config(plot_config);

    bench_spartan_prover(&mut group);
    group.finish();
}

criterion_group!(spartan_prover_benches, spartan_prover_benchmarks);
criterion_main!(spartan_prover_benches);
