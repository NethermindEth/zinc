#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use ark_std::{
    iter::{Product, Sum},
    iterable::Iterable,
};
use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration,
};
use zinc::{
    big_int,
    field::{ConfigRef, RandomField},
    field_config, random_field,
};

fn bench_random_field(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config = field_config!(695962179703626800597079116051991347);
    let field_config = ConfigRef::from(&config);

    let field_elem = random_field!(695962179703, 4, field_config);
    group.bench_with_input(
        BenchmarkId::new("Multiply", "Random128BitFieldElement"),
        &field_elem,
        |b, unop_elem| {
            b.iter(|| {
                for _ in 0..10000 {
                    let _ = black_box(*unop_elem * *unop_elem);
                }
            });
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Addition", "Random128BitFieldElement"),
        &field_elem,
        |b, unop_elem| {
            b.iter(|| {
                for _ in 0..10000 {
                    let _ = black_box(*unop_elem + *unop_elem);
                }
            });
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Division", "Random128BitFieldElement"),
        &field_elem,
        |b, unop_elem| {
            b.iter(|| {
                for _ in 0..10000 {
                    let _ = black_box(*unop_elem / *unop_elem);
                }
            });
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Negation", "Random128BitFieldElement"),
        &field_elem,
        |b, unop_elem| {
            b.iter(|| {
                for _ in 0..10000 {
                    let _ = black_box(-*unop_elem);
                }
            });
        },
    );

    let v = vec![field_elem; 10];

    group.bench_with_input(
        BenchmarkId::new("Sum", "Random128BitFieldElement"),
        &v,
        |b, v| {
            b.iter(|| {
                for _ in 0..10000 {
                    let _ = black_box(RandomField::sum(v.iter()));
                }
            });
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Product", "Random128BitFieldElement"),
        &v,
        |b, v| {
            b.iter(|| {
                for _ in 0..10000 {
                    let _ = black_box(RandomField::product(v.iter()));
                }
            });
        },
    );
}

pub fn field_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);

    let mut group = c.benchmark_group("RandomFieldArithmetic");
    group.plot_config(plot_config);

    bench_random_field(&mut group);
    group.finish();
}

criterion_group!(benches, field_benchmarks);
criterion_main!(benches);
