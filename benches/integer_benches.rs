#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use ark_ff::UniformRand;
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};

use zinc::biginteger::BigInteger256;

fn bench_big_integer(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let mut rng = ark_std::test_rng();
    let biginteger = BigInteger256::rand(&mut rng);
    let small_integer = u32::rand(&mut rng);

    group.bench_with_input(
        BenchmarkId::new("Multiply", "BigInteger256"),
        &(biginteger, small_integer),
        |b, (biginteger, small_integer)| {
            b.iter(move || {
                for _ in 0..10000 {
                    biginteger.clone().muln(*small_integer);
                }
            });
        },
    );
}

pub fn field_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);

    let mut group = c.benchmark_group("BigIntegerArithmetic");
    group.plot_config(plot_config);

    bench_big_integer(&mut group);
    group.finish();
}

criterion_group!(benches, field_benchmarks);
criterion_main!(benches);
