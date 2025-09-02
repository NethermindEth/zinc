#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]
use criterion::{
    AxisScale, BenchmarkId, Criterion, PlotConfiguration, criterion_group, criterion_main,
};
use crypto_bigint::{U192, const_monty_params, modular::ConstMontyParams};
use rand::rng;
use zinc::{
    field::{RandomField, WORD_FACTOR},
    sumcheck::{
        MLSumcheck,
        utils::{rand_poly, rand_poly_comb_fn},
    },
    transcript::KeccakTranscript,
};

const N: usize = 3 * WORD_FACTOR;

fn run_sumcheck<MOD: ConstMontyParams<LIMBS>, const LIMBS: usize>() {
    let nvars = 20;
    let mut rng = rng();
    let ((poly_mles, poly_degree), products, sum) = rand_poly(nvars, (2, 5), 7, &mut rng).unwrap();

    let comb_fn = |vals: &[RandomField<MOD, LIMBS>]| -> RandomField<MOD, LIMBS> {
        rand_poly_comb_fn(vals, &products)
    };

    let mut transcript = KeccakTranscript::new();
    let (proof, _) =
        MLSumcheck::prove_as_subprotocol(&mut transcript, poly_mles, nvars, poly_degree, comb_fn);

    let mut transcript = KeccakTranscript::default();
    let _ = MLSumcheck::verify_as_subprotocol(&mut transcript, nvars, poly_degree, sum, &proof);
}

fn bench_sumcheck_1(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000EB58CBFB80010B1661534393E6C4FA11"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 1"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}

fn bench_sumcheck_2(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000F28CB06184A041766F29E7EFADF3A53F"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 2"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}
fn bench_sumcheck_3(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000D6B1BC9EB22A63FF926620E0604D9E75"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 3"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}
fn bench_sumcheck_4(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000F28CB06184A041766F29E7EFADF3A53F"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 4"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}
fn bench_sumcheck_5(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000B218E1FDE7571D007C3343422146865B"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 5"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}
fn bench_sumcheck_6(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000A54574AD88C665349B6C30ABF994B5C5"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 6"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}

fn bench_sumcheck_7(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000E2F3BC2330C9C8393D5051C71700374B"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 7"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}
fn bench_sumcheck_8(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000956B1628819CE25955477B33550C4DF7"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 8"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}
fn bench_sumcheck_9(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "000000000000000092BACDB8268D524F3DE8070A0849F3C1"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 9"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
}
fn bench_sumcheck_10(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    const_monty_params!(
        ModP,
        U192,
        "0000000000000000C988DE29E4E94A1F0DB0CB829F7A5D5B"
    );

    group.bench_function(BenchmarkId::new("Sumcheck", "Prime 10"), |b| {
        b.iter(|| {
            run_sumcheck::<ModP, N>();
        });
    });
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
