use ark_std::iterable::Iterable;
use i256::I256;
use sha3::{digest::Output, Digest, Keccak256};

use crate::{
    poly_z::mle::DenseMultilinearExtension,
    zip::{
        code::{LinearCodes, ZipSpec},
        utils::{div_ceil, num_threads, parallelize, parallelize_iter},
        Error,
    },
};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment},
    utils::validate_input,
};

impl<const N: usize, S> MultilinearZip<N, S>
where
    S: ZipSpec,
{
    pub fn commit(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, Error> {
        validate_input("commit", pp.num_vars(), [poly], None)?;

        let row_len = pp.zip().row_len();
        let codeword_len = pp.zip().codeword_len();
        let merkle_depth = codeword_len.next_power_of_two().ilog2() as usize;

        let mut rows = vec![0i64; pp.num_rows() * codeword_len];
        let mut hashes = vec![Output::<Keccak256>::default(); (2 << merkle_depth) - 1];

        let rows = Self::encode_rows(pp, codeword_len, row_len, &mut rows, poly);
        Self::compute_column_hashes(&mut hashes, codeword_len, &rows);

        Self::merklize_column_hashes(merkle_depth, &mut hashes);

        let (intermediate_hashes, root) = {
            let mut intermediate_hashes = hashes;
            let root = intermediate_hashes.pop().unwrap();
            (intermediate_hashes, root)
        };

        Ok(MultilinearZipCommitment::new(
            rows,
            intermediate_hashes,
            root,
        ))
    }

    pub fn batch_commit<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension>,
    ) -> Result<Vec<Self::Commitment>, Error> {
        polys.iter().map(|poly| Self::commit(pp, poly)).collect()
    }

    fn encode_rows(
        pp: &Self::ProverParam,
        codeword_len: usize,
        row_len: usize,
        rows: &mut [i64],
        poly: &Self::Polynomial,
    ) -> Vec<I256> {
        let chunk_size = div_ceil(pp.num_rows(), num_threads());
        let mut encoded_rows = vec![I256::default(); rows.len()];

        parallelize_iter(
            encoded_rows
                .chunks_exact_mut(chunk_size * codeword_len)
                .zip(poly.evaluations.chunks_exact(chunk_size * row_len)),
            |(encoded_chunk, evals)| {
                for (row, evals) in encoded_chunk
                    .chunks_exact_mut(codeword_len)
                    .zip(evals.chunks_exact(row_len))
                {
                    let mut temp_row = vec![0i64; codeword_len];
                    temp_row[..evals.len()].copy_from_slice(evals);
                    pp.zip().encode(&temp_row);

                    for (i, val) in temp_row.iter().enumerate() {
                        row[i] = I256::from(*val);
                    }
                }
            },
        );

        encoded_rows
    }

    fn compute_column_hashes(hashes: &mut [Output<Keccak256>], codeword_len: usize, rows: &[I256]) {
        parallelize(&mut hashes[..codeword_len], |(hashes, start)| {
            let mut hasher = Keccak256::new();
            for (hash, column) in hashes.iter_mut().zip(start..) {
                rows.iter()
                    .skip(column)
                    .step_by(codeword_len)
                    .for_each(|item| {
                        <Keccak256 as sha3::digest::Update>::update(
                            &mut hasher,
                            &item.to_be_bytes(),
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
