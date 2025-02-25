use ark_ff::Zero;
use ark_std::iterable::Iterable;
use sha3::{digest::Output, Digest, Keccak256};

use crate::{
    brakedown::{
        code::{BrakedownSpec, LinearCodes},
        utils::{div_ceil, num_threads, parallelize, parallelize_iter},
        Error,
    },
    field::RandomField as F,
    poly::mle::DenseMultilinearExtension,
};

use super::{
    structs::{MultilinearBrakedown, MultilinearBrakedownCommitment},
    utils::validate_input,
};

impl<const N: usize, S> MultilinearBrakedown<N, S>
where
    S: BrakedownSpec,
{
    pub fn commit(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, Error> {
        validate_input("commit", pp.num_vars(), [poly], None)?;

        let row_len = pp.brakedown().row_len();
        let codeword_len = pp.brakedown().codeword_len();
        let merkle_depth = codeword_len.next_power_of_two().ilog2() as usize;

        let mut rows = vec![F::zero(); pp.num_rows() * codeword_len];
        let mut hashes = vec![Output::<Keccak256>::default(); (2 << merkle_depth) - 1];

        Self::encode_rows(pp, codeword_len, row_len, &mut rows, poly);
        Self::compute_column_hashes(&mut hashes, codeword_len, &rows);

        Self::merklize_column_hashes(merkle_depth, &mut hashes);

        let (intermediate_hashes, root) = {
            let mut intermediate_hashes = hashes;
            let root = intermediate_hashes.pop().unwrap();
            (intermediate_hashes, root)
        };

        Ok(MultilinearBrakedownCommitment::new(
            rows,
            intermediate_hashes,
            root,
        ))
    }

    pub fn batch_commit<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension<N>>,
    ) -> Result<Vec<Self::Commitment>, Error> {
        polys.iter().map(|poly| Self::commit(pp, poly)).collect()
    }
    fn encode_rows(
        pp: &Self::ProverParam,
        codeword_len: usize,
        row_len: usize,
        rows: &mut [F<N>],
        poly: &Self::Polynomial,
    ) {
        let chunk_size = div_ceil(pp.num_rows(), num_threads());
        parallelize_iter(
            rows.chunks_exact_mut(chunk_size * codeword_len)
                .zip(poly.evaluations.chunks_exact(chunk_size * row_len)),
            |(rows, evals)| {
                for (row, evals) in rows
                    .chunks_exact_mut(codeword_len)
                    .zip(evals.chunks_exact(row_len))
                {
                    row[..evals.len()].copy_from_slice(evals);
                    pp.brakedown().encode(row);
                }
            },
        );
    }

    fn compute_column_hashes(hashes: &mut [Output<Keccak256>], codeword_len: usize, rows: &[F<N>]) {
        parallelize(&mut hashes[..codeword_len], |(hashes, start)| {
            let mut hasher = Keccak256::new();
            for (hash, column) in hashes.iter_mut().zip(start..) {
                rows.iter()
                    .skip(column)
                    .step_by(codeword_len)
                    .for_each(|item| {
                        <Keccak256 as sha3::digest::Update>::update(
                            &mut hasher,
                            &item.value().to_bytes_be(),
                        )
                    });
                hasher.finalize_into_reset(hash);
            }
        });
    }

    fn merklize_column_hashes(depth: usize, hashes: &mut [Output<Keccak256>]) {
        let mut offset = 0;
        for width in (1..=depth).rev().map(|depth| 1 << depth) {
            let (input, output) = hashes[offset..].split_at_mut(width);
            //	    let num_threads = env::var("RAYON_NUM_THREADS").unwrap();
            let chunk_size = div_ceil(output.len(), num_threads());
            parallelize_iter(
                input
                    .chunks(2 * chunk_size)
                    .zip(output.chunks_mut(chunk_size)),
                |(input, output)| {
                    let mut hasher = Keccak256::new();

                    for (input, output) in input.chunks_exact(2).zip(output.iter_mut()) {
                        <Keccak256 as sha3::digest::Update>::update(&mut hasher, &input[0]);
                        <Keccak256 as sha3::digest::Update>::update(&mut hasher, &input[1]);
                        hasher.finalize_into_reset(output);
                    }
                },
            );
            offset += width;
        }
    }
}
