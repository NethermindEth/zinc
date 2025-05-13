use std::str::FromStr;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use zinc::field_config::ConfigPtr;
use zinc::{
    biginteger::BigInt,
    ccs::test_utils::get_dummy_ccs_Z_from_z_length,
    field_config::FieldConfig,
    transcript::KeccakTranscript,
    zinc::{
        prover::SpartanProver,
        structs::{ZincProver, ZincVerifier},
        verifier::SpartanVerifier,
    },
    zip::code::ZipSpec1,
};

fn benchmark_spartan_prover<const N: usize>(c: &mut Criterion, config: ConfigPtr<N>, prime: &str) {
    let mut group = c.benchmark_group(format!("spartan_prover for {} prime", prime));
    let mut rng = ark_std::test_rng();

    let prover = ZincProver::<N, _> {
        // If we are keeping primes around 128 bits we should stay with N = 3 hardcoded
        data: std::marker::PhantomData::<ZipSpec1>,
    };

    for size in [12, 13, 14, 15, 16] {
        let n = 1 << size;
        let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);

        let (z_ccs, z_mle, ccs_f, statement_f) =
            ZincProver::<N, ZipSpec1>::prepare_for_random_field_piop(
                &statement, &wit, &ccs, config,
            )
            .expect("Failed to prepare for random field PIOP");

        group.bench_function(format!("n={}", n), |b| {
            b.iter_batched(
                KeccakTranscript::new,
                |mut prover_transcript| {
                    black_box(
                        SpartanProver::<N>::prove(
                            &prover,
                            &statement_f,
                            &z_ccs,
                            &z_mle,
                            &ccs_f,
                            &mut prover_transcript,
                            config,
                        )
                        .expect("Proof generation failed"),
                    )
                },
                criterion::BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn benchmark_spartan_verifier<const N: usize>(
    c: &mut Criterion,
    config: ConfigPtr<N>,
    prime: &str,
) {
    let mut group = c.benchmark_group(format!("spartan_verifier for {} prime", prime));
    let mut rng = ark_std::test_rng();

    let prover = ZincProver::<N, _> {
        data: std::marker::PhantomData::<ZipSpec1>,
    };

    let verifier = ZincVerifier::<N, _> {
        data: std::marker::PhantomData::<ZipSpec1>,
    };

    for size in [12, 13, 14, 15, 16] {
        let n = 1 << size;
        let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);
        let mut prover_transcript = KeccakTranscript::new();

        let (z_ccs, z_mle, ccs_f, statement_f) =
            ZincProver::<N, ZipSpec1>::prepare_for_random_field_piop(
                &statement, &wit, &ccs, config,
            )
            .expect("Failed to prepare for random field PIOP");

        let (spartan_proof, _) = SpartanProver::<N>::prove(
            &prover,
            &statement_f,
            &z_ccs,
            &z_mle,
            &ccs_f,
            &mut prover_transcript,
            config,
        )
        .expect("Failed to generate Spartan proof");

        group.bench_function(format!("n={}", n), |b| {
            b.iter_batched(
                KeccakTranscript::new,
                |mut verifier_transcript| {
                    black_box(
                        SpartanVerifier::<N>::verify(
                            &verifier,
                            &spartan_proof,
                            &mut verifier_transcript,
                            &ccs_f,
                            config.reference().expect("Field config cannot be none"),
                        )
                        .expect("Proof verification failed"),
                    )
                },
                criterion::BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn run_benches(c: &mut Criterion) {
    // Using a 256-bit prime field
    let config = ConfigPtr::from(&FieldConfig::new(
        BigInt::<4>::from_str(
            "115792089237316195423570985008687907853269984665640564039457584007913129639747",
        )
        .unwrap(),
    ));

    benchmark_spartan_prover::<4>(c, config, "256");
    benchmark_spartan_verifier::<4>(c, config, "256");

    let stark_config = ConfigPtr::from(&FieldConfig::new(
        BigInt::<4>::from_str(
            "3618502788666131213697322783095070105623107215331596699973092056135872020481",
        )
        .unwrap(),
    ));

    benchmark_spartan_prover::<4>(c, stark_config, "stark");
    benchmark_spartan_verifier::<4>(c, stark_config, "stark");
}

criterion_group!(benches, run_benches);
criterion_main!(benches);
