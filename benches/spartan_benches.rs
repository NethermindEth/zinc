use criterion::{black_box, criterion_group, criterion_main, Criterion};
use zinc::{
    ccs::test_utils::get_dummy_ccs_Z_from_z_length,
    define_random_field_zip_types,
    field::{ConfigRef, RandomField},
    field_config, implement_random_field_zip_types,
    traits::{ConfigReference, Field, FieldMap, Integer, ZipTypes},
    transcript::KeccakTranscript,
    zinc::{
        prover::SpartanProver,
        structs::{ZincProver, ZincVerifier},
        verifier::SpartanVerifier,
    },
    zip::code::DefaultLinearCodeSpec,
};

const INT_LIMBS: usize = 1;
const FIELD_LIMBS: usize = 4;

define_random_field_zip_types!();
implement_random_field_zip_types!(INT_LIMBS);

fn benchmark_spartan_prover<ZT: ZipTypes, F: Field>(c: &mut Criterion, config: F::R, prime: &str)
where
    for<'a> ZT::N: From<&'a F::I>,
    for<'a> F::I: From<&'a <ZT::N as Integer>::I>,
    for<'a> F::I: From<&'a ZT::N>,
    ZT::N: FieldMap<F, Output = F>,
{
    let mut group = c.benchmark_group(format!("spartan_prover for {prime} prime"));
    let mut rng = ark_std::test_rng();

    // If we are keeping primes around 128 bits we should stay with N = 3 hardcoded
    let prover = ZincProver::<ZT, F, _>::new(DefaultLinearCodeSpec);

    for size in [12, 13, 14, 15, 16] {
        let n = 1 << size;
        let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);

        let (z_ccs, z_mle, ccs_f, statement_f) =
            ZincProver::<ZT, F, DefaultLinearCodeSpec>::prepare_for_random_field_piop(
                &statement, &wit, &ccs, config,
            )
            .expect("Failed to prepare for random field PIOP");

        group.bench_function(format!("n={n}"), |b| {
            b.iter_batched(
                KeccakTranscript::new,
                |mut prover_transcript| {
                    black_box(
                        SpartanProver::<_, _>::prove(
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

fn benchmark_spartan_verifier<ZT: ZipTypes, F: Field>(c: &mut Criterion, config: F::R, prime: &str)
where
    for<'a> ZT::N: From<&'a F::I>,
    for<'a> F::I: From<&'a <ZT::N as Integer>::I>,
    for<'a> F::I: From<&'a ZT::N>,
    ZT::N: FieldMap<F, Output = F>,
{
    let mut group = c.benchmark_group(format!("spartan_verifier for {prime} prime"));
    let mut rng = ark_std::test_rng();

    let prover = ZincProver::<ZT, F, _>::new(DefaultLinearCodeSpec);

    let verifier = ZincVerifier::<ZT, F, _>::new(DefaultLinearCodeSpec);

    for size in [12, 13, 14, 15, 16] {
        let n = 1 << size;
        let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);
        let mut prover_transcript = KeccakTranscript::new();

        let (z_ccs, z_mle, ccs_f, statement_f) =
            ZincProver::<ZT, F, DefaultLinearCodeSpec>::prepare_for_random_field_piop(
                &statement, &wit, &ccs, config,
            )
            .expect("Failed to prepare for random field PIOP");

        let (spartan_proof, _) = SpartanProver::<_, _>::prove(
            &prover,
            &statement_f,
            &z_ccs,
            &z_mle,
            &ccs_f,
            &mut prover_transcript,
            config,
        )
        .expect("Failed to generate Spartan proof");
        config.reference().expect("Field config cannot be none");
        group.bench_function(format!("n={n}"), |b| {
            b.iter_batched(
                KeccakTranscript::new,
                |mut verifier_transcript| {
                    black_box(
                        SpartanVerifier::<F>::verify(
                            &verifier,
                            &spartan_proof,
                            &ccs_f,
                            &mut verifier_transcript,
                            config,
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
    type ZT = RandomFieldZipTypes<INT_LIMBS>;
    type F<'a> = RandomField<'a, FIELD_LIMBS>;

    // Using a 256-bit prime field
    let config = field_config!(
        115792089237316195423570985008687907853269984665640564039457584007913129639747,
        4
    );
    let config = ConfigRef::from(&config);

    benchmark_spartan_prover::<ZT, F>(c, config, "256");
    benchmark_spartan_verifier::<ZT, F>(c, config, "256");

    let stark_config = field_config!(
        3618502788666131213697322783095070105623107215331596699973092056135872020481,
        4
    );
    let stark_config = ConfigRef::from(&stark_config);

    benchmark_spartan_prover::<ZT, F>(c, stark_config, "stark");
    benchmark_spartan_verifier::<ZT, F>(c, stark_config, "stark");
}

criterion_group!(benches, run_benches);
criterion_main!(benches);
