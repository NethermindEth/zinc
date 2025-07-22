use std::marker::PhantomData;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use zinc::{
    ccs::test_utils::get_dummy_ccs_Z_from_z_length, field_config::ConfigRef, zinc::prelude::*,
};

const INT_LIMBS: usize = 1;
const FIELD_LIMBS: usize = 4;

fn benchmark_zinc_prove(c: &mut Criterion) {
    c.bench_function("Zinc/Prove", |b| {
        b.iter(|| {
            let prover = ZincProver::<INT_LIMBS, FIELD_LIMBS, ZipSpec1> { data: PhantomData };
            let mut prover_transcript = KeccakTranscript::new();
            let mut rng = ark_std::test_rng();

            let (_, ccs, statement, witness) =
                get_dummy_ccs_Z_from_z_length::<INT_LIMBS>(8, &mut rng);

            let field_config = draw_random_field::<INT_LIMBS, FIELD_LIMBS>(
                &statement.public_input,
                &mut prover_transcript,
            );
            let config_ref = ConfigRef::from(&field_config);

            let _ = black_box(
                prover
                    .prove(
                        &statement,
                        &witness,
                        &mut prover_transcript,
                        &ccs,
                        config_ref,
                    )
                    .expect("Proof generation failed"),
            );
        });
    });
}

fn benchmark_zinc_verify(c: &mut Criterion) {
    c.bench_function("Zinc/Verify", |b| {
        let prover = ZincProver::<INT_LIMBS, FIELD_LIMBS, ZipSpec1> { data: PhantomData };
        let verifier = ZincVerifier::<INT_LIMBS, FIELD_LIMBS, ZipSpec1> { data: PhantomData };

        let mut prover_transcript = KeccakTranscript::new();
        let mut rng = ark_std::test_rng();
        let (_, ccs, statement, witness) = get_dummy_ccs_Z_from_z_length::<INT_LIMBS>(8, &mut rng);

        let field_config = draw_random_field::<INT_LIMBS, FIELD_LIMBS>(
            &statement.public_input,
            &mut prover_transcript,
        );
        let config_ref = ConfigRef::from(&field_config);

        let proof = prover
            .prove(
                &statement,
                &witness,
                &mut prover_transcript,
                &ccs,
                config_ref,
            )
            .expect("Proof generation failed");
        b.iter_batched(
            || proof.clone(),
            |proof| {
                let mut verifier_transcript = KeccakTranscript::new();

                verifier
                    .verify(
                        &statement,
                        proof,
                        &mut verifier_transcript,
                        &ccs,
                        config_ref,
                    )
                    .expect("Proof verification failed");
            },
            criterion::BatchSize::LargeInput,
        );
    });
}

criterion_group!(benches, benchmark_zinc_prove, benchmark_zinc_verify);
criterion_main!(benches);
