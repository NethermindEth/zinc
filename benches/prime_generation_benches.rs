#![allow(non_local_definitions)]
#![allow(clippy::eq_op)]

use ark_std::hint::black_box;
use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use zinc::{
    field::{BigInt, RandomField},
    prime_gen,
    transcript::KeccakTranscript,
};

fn bench_prime_generation(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>) {
    let hasher = KeccakTranscript::new();
    const N: usize = 3;
    group.bench_with_input(BenchmarkId::new("PrimeGen", "196bits"), &hasher, |b, _| {
        b.iter(|| {
            let _: BigInt<N> =
                black_box(prime_gen::get_prime::<RandomField<N>>(&mut hasher.clone()));
        });
    });
}

pub fn field_benchmarks(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);

    let mut group = c.benchmark_group("PrimeGen");
    group.plot_config(plot_config);

    bench_prime_generation(&mut group);
    group.finish();
}

criterion_group!(benches, field_benchmarks);
criterion_main!(benches);
