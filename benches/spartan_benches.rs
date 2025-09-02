use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use crypto_bigint::{U256, const_monty_params};
use rand::rng;
use zinc::{
    ccs::get_dummy_ccs_Z_from_z_length,
    field::{F256, WORD_FACTOR},
    traits::Field,
    transcript::KeccakTranscript,
    zinc::{
        prover::SpartanProver,
        structs::{ZincProver, ZincVerifier},
        verifier::SpartanVerifier,
    },
    zip::code::DefaultLinearCodeSpec,
};

const INT_LIMBS: usize = WORD_FACTOR;
const FIELD_LIMBS: usize = 4 * WORD_FACTOR;

const N: usize = INT_LIMBS;
const L: usize = INT_LIMBS * 2;
const K: usize = INT_LIMBS * 4;
const M: usize = INT_LIMBS * 8;

const_monty_params!(
    ModP,
    U256,
    "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF43"
);
const_monty_params!(
    ModQ,
    U256,
    "0800000000000011000000000000000000000000000000000000000000000001"
);

fn benchmark_spartan_prover<
    const N: usize,
    const L: usize,
    const K: usize,
    const M: usize,
    F: Field<FIELD_LIMBS>,
>(
    c: &mut Criterion,
    prime: &str,
) {
    let mut group = c.benchmark_group(format!("spartan_prover for {prime} prime"));
    let mut rng = rng();

    // If we are keeping primes around 128 bits we should stay with N = 3 hardcoded
    let prover = ZincProver::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);

    for size in [12, 13, 14, 15, 16] {
        let n = 1 << size;
        let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);

        let (z_ccs, z_mle, ccs_f, statement_f) = ZincProver::<
            N,
            L,
            K,
            M,
            F,
            FIELD_LIMBS,
            DefaultLinearCodeSpec,
        >::prepare_for_random_field_piop(
            &statement, &wit, &ccs
        )
        .expect("Failed to prepare for random field PIOP");

        group.bench_function(format!("n={n}"), |b| {
            b.iter_batched(
                KeccakTranscript::new,
                |mut prover_transcript| {
                    black_box(
                        SpartanProver::prove(
                            &prover,
                            &statement_f,
                            &z_ccs,
                            &z_mle,
                            &ccs_f,
                            &mut prover_transcript,
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

fn benchmark_spartan_verifier<
    const N: usize,
    const L: usize,
    const K: usize,
    const M: usize,
    F: Field<FIELD_LIMBS>,
>(
    c: &mut Criterion,
    prime: &str,
) {
    let mut group = c.benchmark_group(format!("spartan_verifier for {prime} prime"));
    let mut rng = rng();

    let prover = ZincProver::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);

    let verifier = ZincVerifier::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);

    for size in [12, 13, 14, 15, 16] {
        let n = 1 << size;
        let (_, ccs, statement, wit) = get_dummy_ccs_Z_from_z_length(n, &mut rng);
        let mut prover_transcript = KeccakTranscript::new();

        let (z_ccs, z_mle, ccs_f, statement_f) = ZincProver::<
            N,
            L,
            K,
            M,
            F,
            FIELD_LIMBS,
            DefaultLinearCodeSpec,
        >::prepare_for_random_field_piop(
            &statement, &wit, &ccs
        )
        .expect("Failed to prepare for random field PIOP");

        let (spartan_proof, _) = SpartanProver::prove(
            &prover,
            &statement_f,
            &z_ccs,
            &z_mle,
            &ccs_f,
            &mut prover_transcript,
        )
        .expect("Failed to generate Spartan proof");
        group.bench_function(format!("n={n}"), |b| {
            b.iter_batched(
                KeccakTranscript::new,
                |mut verifier_transcript| {
                    black_box(
                        SpartanVerifier::<F, FIELD_LIMBS>::verify(
                            &verifier,
                            &spartan_proof,
                            &ccs_f,
                            &mut verifier_transcript,
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
    benchmark_spartan_prover::<N, L, K, M, F256<ModP>>(c, "256");
    benchmark_spartan_verifier::<N, L, K, M, F256<ModP>>(c, "256");

    benchmark_spartan_prover::<N, L, K, M, F256<ModP>>(c, "stark");
    benchmark_spartan_verifier::<N, L, K, M, F256<ModP>>(c, "stark");
}

criterion_group!(benches, run_benches);
criterion_main!(benches);
