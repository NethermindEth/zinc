#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use ark_ff::UniformRand;
use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration,
};

use rand::thread_rng;
use zinc::biginteger::BigInteger256;

fn bench_big_integer(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    group.bench_function(BenchmarkId::new("Multiply", "BigInteger256"), |b| {
        b.iter_batched(
            || {
                let mut rng = thread_rng();
                let biginteger = BigInteger256::rand(&mut rng);
                let small_integer = u32::rand(&mut rng);
                (biginteger, small_integer)
            },
            |(biginteger, small_integer)| {
                for _ in 0..100000 {
                    black_box(biginteger).muln(black_box(small_integer));
                }
            },
            criterion::BatchSize::SmallInput,
        );
    });
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
