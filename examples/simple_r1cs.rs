use std::marker::PhantomData;

use ark_std::log2;
use crypto_bigint::Int;
use zinc::zinc::prelude::*;

// R1CS for: x^3 + x + 5 = y (example from article
// https://www.vitalik.ca/general/2016/12/10/qap.html )
fn main() {
    // Example code goes here
    const N: usize = 4;
    let prover = ZincProver::<N, _> {
        data: PhantomData::<ZipSpec1>,
    };
    let mut prover_transcript = KeccakTranscript::new();

    let (ccs, statement, witness) = get_ccs_stuff(3);
    let field_config = draw_random_field::<N>(&statement.public_input, &mut prover_transcript);
    let proof = prover
        .prove(
            &statement,
            &witness,
            &mut prover_transcript,
            &ccs,
            &field_config,
        )
        .expect("Proof generation failed");

    let verifier = ZincVerifier::<N, _> {
        data: PhantomData::<ZipSpec1>,
    };

    let mut verifier_transcript = KeccakTranscript::new();
    verifier
        .verify(
            &statement,
            proof,
            &mut verifier_transcript,
            &ccs,
            &field_config,
        )
        .expect("Proof verification failed");
}

fn get_ccs<const N: usize>() -> CCS_Z<N> {
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
    }
}

#[allow(non_snake_case)]
fn get_ccs_statement<const N: usize>(input: i64) -> Statement_Z<N> {
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
    let public_input = vec![Int::<N>::from_i64(input)];
    Statement_Z {
        constraints,
        public_input,
    }
}

fn get_witness<const N: usize>(input: i64) -> Witness_Z<N> {
    Witness_Z::new(vec![
        Int::<N>::from_i64(input * input * input + input + 5), // x^3 + x + 5
        Int::<N>::from_i64(input * input),                     // x^2
        Int::<N>::from_i64(input * input * input),             // x^2 * x
        Int::<N>::from_i64(input * input * input + input),     // x^3 + x
    ])
}

fn get_ccs_stuff<const N: usize>(input: i64) -> (CCS_Z<N>, Statement_Z<N>, Witness_Z<N>) {
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
