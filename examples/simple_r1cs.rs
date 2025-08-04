use ark_std::{log2, marker::PhantomData};
use zinc::{crypto_int::Int, field::RandomField, field_config::ConfigRef, zinc::prelude::*};

// R1CS for: x^3 + x + 5 = y (example from article
// https://www.vitalik.ca/general/2016/12/10/qap.html )
fn main() {
    // Example code goes here
    const FIELD_LIMBS: usize = 4;
    const INT_LIMBS: usize = 1;
    let prover =
        ZincProver::<Int<INT_LIMBS>, RandomField<FIELD_LIMBS>, ZipSpec1> { data: PhantomData };
    let mut prover_transcript = KeccakTranscript::new();

    let (ccs, statement, witness) = get_ccs_stuff(3);
    let field_config = draw_random_field::<Int<INT_LIMBS>, RandomField<FIELD_LIMBS>>(
        &statement.public_input,
        &mut prover_transcript,
    );

    let config_ref = ConfigRef::from(&field_config);

    let proof = prover
        .prove::<Int<{ INT_LIMBS * 2 }>, Int<{ INT_LIMBS * 4 }>, Int<{ INT_LIMBS * 8 }>>(
            &statement,
            &witness,
            &mut prover_transcript,
            &ccs,
            ConfigRef::from(&field_config),
        )
        .expect("Proof generation failed");

    let verifier =
        ZincVerifier::<Int<INT_LIMBS>, RandomField<FIELD_LIMBS>, ZipSpec1> { data: PhantomData };

    let mut verifier_transcript = KeccakTranscript::new();
    verifier
        .verify::<Int<{ INT_LIMBS * 2 }>, Int<{ INT_LIMBS * 4 }>, Int<{ INT_LIMBS * 8 }>>(
            &statement,
            proof,
            &mut verifier_transcript,
            &ccs,
            config_ref,
        )
        .expect("Proof verification failed");
}

fn get_ccs<const N: usize>() -> CCS_Z<Int<N>> {
    let m = 4;
    let n = 6;
    CCS_Z {
        m,
        n,
        l: 1,
        t: 3,
        q: 2,
        d: 2,
        s: log2(m) as usize,
        s_prime: log2(n) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![1, -1],
        _phantom: PhantomData,
    }
}

#[allow(non_snake_case)]
fn get_ccs_statement<const N: usize>(input: i64) -> Statement_Z<Int<N>> {
    let A = to_Z_matrix(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 1, 0, 0],
        vec![1, 0, 0, 0, 1, 0],
        vec![0, 5, 0, 0, 0, 1],
    ]);
    let B = to_Z_matrix(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
    ]);
    let C = to_Z_matrix(vec![
        vec![0, 0, 0, 1, 0, 0],
        vec![0, 0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 0, 1],
        vec![0, 0, 1, 0, 0, 0],
    ]);
    let constraints = vec![A, B, C];
    let public_input = vec![Int::<N>::from(input)];
    Statement_Z {
        constraints,
        public_input,
    }
}

fn get_witness<const N: usize>(input: i64) -> Witness_Z<Int<N>> {
    Witness_Z::new(
        [
            input.pow(3) + input + 5, // x^3 + x + 5
            input.pow(2),             // x^2
            input.pow(3) * input,     // x^2 * x
            input.pow(3) + input,     // x^3 + x
        ]
        .iter()
        .cloned()
        .map(Int::from)
        .collect(),
    )
}

fn get_ccs_stuff<const N: usize>(
    input: i64,
) -> (CCS_Z<Int<N>>, Statement_Z<Int<N>>, Witness_Z<Int<N>>) {
    let mut ccs = get_ccs();
    let mut statement = get_ccs_statement(input);
    let matrices = statement.constraints.clone();
    let witness = get_witness(input);
    let z = statement.get_z_vector(&witness.w_ccs);
    ccs.check_relation(&matrices, &z)
        .expect("Failed to check relation over Integer Ring");
    let len = usize::max(ccs.m.next_power_of_two(), ccs.n.next_power_of_two());
    ccs.pad(&mut statement, len);
    (ccs, statement, witness)
}
