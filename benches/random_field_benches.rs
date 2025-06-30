#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use ark_std::{hint::black_box, iterable::Iterable};
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use std::{
    iter::{Product, Sum},
    str::FromStr,
};
use zinc::field::conversion::FieldMap;

use zinc::{
    biginteger::BigInteger256,
    field::RandomField,
    field_config::{ConfigRef, FieldConfig},
};

fn bench_random_field(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let config =
        FieldConfig::new(BigInteger256::from_str("695962179703626800597079116051991347").unwrap());
    let field_config = ConfigRef::from(&config);

    let bigint = BigInteger256::from_str("695962179703").unwrap();

    let field_elem = bigint.map_to_field(field_config);
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
