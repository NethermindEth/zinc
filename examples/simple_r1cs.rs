use ark_std::{log2, marker::PhantomData};
use crypto_bigint::{Int, U256, const_monty_params};
use zinc::{
    ccs::{Arith, CcsZ, Instance, Statement, Witness},
    field::{F256, WORD_FACTOR},
    sparse_matrix::to_Z_matrix,
    transcript::KeccakTranscript,
    zinc::prelude::{DefaultLinearCodeSpec, Prover, Verifier, ZincProver, ZincVerifier},
};

// Word factor addresses different WORD size on 32 and 64bit architectures
const FIELD_LIMBS: usize = 4 * WORD_FACTOR;

const INT_LIMBS: usize = 1 * WORD_FACTOR;

const N: usize = INT_LIMBS;
const L: usize = INT_LIMBS * 2;
const K: usize = INT_LIMBS * 4;
const M: usize = INT_LIMBS * 4;

// R1CS for: x^3 + x + 5 = y (example from article
// https://www.vitalik.ca/general/2016/12/10/qap.html )
fn main() {
    // Example code goes here
    let prover = ZincProver::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);
    let mut prover_transcript = KeccakTranscript::new();

    let (ccs, statement, witness) = get_ccs_stuff(3);

    let proof = prover
        .prove(&statement, &witness, &mut prover_transcript, &ccs)
        .expect("Proof generation failed");

    let verifier = ZincVerifier::<N, L, K, M, F, FIELD_LIMBS, _>::new(DefaultLinearCodeSpec);

    let mut verifier_transcript = KeccakTranscript::new();
    verifier
        .verify(&statement, proof, &mut verifier_transcript, &ccs)
        .expect("Proof verification failed");
}

const_monty_params!(
    ModP,
    U256,
    "0800000000000011000000000000000000000000000000000000000000000001"
);

type F = F256<ModP>;

fn get_ccs<const N: usize>() -> CcsZ<Int<N>> {
    let m = 4;
    let n = 6;
    CcsZ {
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
fn get_ccs_statement<const N: usize>(input: i64) -> Statement<Int<N>> {
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
    Statement {
        constraints,
        public_input,
    }
}

fn get_witness<const N: usize>(input: i64) -> Witness<Int<N>> {
    Witness::new(
        [
            input.pow(3) + input + 5, // x^3 + x + 5
            input.pow(2),             // x^2
            input.pow(2) * input,     // x^2 * x
            input.pow(3) + input,     // x^3 + x
        ]
        .iter()
        .cloned()
        .map(Int::from)
        .collect(),
    )
}

fn get_ccs_stuff<const N: usize>(input: i64) -> (CcsZ<Int<N>>, Statement<Int<N>>, Witness<Int<N>>) {
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
