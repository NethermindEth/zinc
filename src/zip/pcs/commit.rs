use ark_std::iterable::Iterable;
use i256::I512;
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
    structs::{MultilinearZip, MultilinearZipCommitment, MultilinearZipData, ZipTranscript},
    utils::validate_input,
};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub fn commit(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
    ) -> Result<(Self::Data, Self::Commitment), Error> {
        validate_input("commit", pp.num_vars(), [poly], None)?;

        let row_len = pp.zip().row_len();
        let num_rows = row_len; // since initially the coefficient matrix is square, row_len = num_rows
        let codeword_len = pp.zip().codeword_len();
        // We deviate from the paper to merkleize each column instead of each row
        let merkle_depth = num_rows.next_power_of_two().ilog2() as usize;

        let mut hashes =
            vec![Output::<Keccak256>::default(); codeword_len * ((2 << merkle_depth) - 1)];

        let rows = Self::encode_rows(pp, codeword_len, row_len, poly);

        // Transpose rows into columns
        let mut columns = vec![I512::default(); rows.len()];
        for i in 0..pp.num_rows() {
            for j in 0..codeword_len {
                columns[j * pp.num_rows() + i] = rows[i * codeword_len + j];
            }
        }

        // First compute all hashes
        let mut temp_hashes = vec![Output::<Keccak256>::default(); num_rows * codeword_len];
        Self::compute_column_hashes(&mut temp_hashes, &columns);

        // Redistribute the hashes with proper spacing
        for i in 0..num_rows {
            let start_idx = i * ((2 << merkle_depth) - 1);
            hashes[start_idx] = temp_hashes[i];
        }

        for i in 0..codeword_len {
            let start_idx = i * ((2 << merkle_depth) - 1);
            let end_idx = (i + 1) * ((2 << merkle_depth) - 1);
            Self::merklize_column_hashes(merkle_depth, &mut hashes[start_idx..end_idx]);
        }

        let (intermediate_hashes, root) = {
            let mut intermediate_hashes = hashes;
            let root = intermediate_hashes.pop().unwrap();
            (intermediate_hashes, root)
        };

        Ok((
            MultilinearZipData::new(rows, intermediate_hashes, root),
            MultilinearZipCommitment::new(root),
        ))
    }
    #[allow(clippy::type_complexity)]
    pub fn batch_commit<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension>,
    ) -> Result<Vec<(Self::Data, Self::Commitment)>, Error> {
        polys.iter().map(|poly| Self::commit(pp, poly)).collect()
    }

    fn encode_rows(
        pp: &Self::ProverParam,
        codeword_len: usize,
        row_len: usize,
        poly: &Self::Polynomial,
    ) -> Vec<I512> {
        assert_eq!(pp.num_rows(), poly.evaluations.len().isqrt());
        assert_eq!(codeword_len, row_len * 2);
        let rows_per_thread = div_ceil(pp.num_rows(), num_threads());
        let mut encoded_rows = vec![I512::default(); pp.num_rows() * codeword_len];

        parallelize_iter(
            encoded_rows
                .chunks_exact_mut(rows_per_thread * codeword_len)
                .zip(poly.evaluations.chunks_exact(rows_per_thread * row_len)),
            |(encoded_chunk, evals)| {
                for (row, evals) in encoded_chunk
                    .chunks_exact_mut(codeword_len)
                    .zip(evals.chunks_exact(row_len))
                {
                    row.copy_from_slice(pp.zip().encode_i64(evals).as_slice());
                }
            },
        );

        encoded_rows
    }

    fn compute_column_hashes(hashes: &mut [Output<Keccak256>], columns: &[I512]) {
        // TODO:improve this without transposing the rows into columns
        parallelize(hashes, |(hashes, start)| {
            let mut hasher = Keccak256::new();
            for (hash, column) in hashes.iter_mut().zip(start..) {
                // For each column, iterate through all rows at that column position
                <Keccak256 as sha3::digest::Update>::update(
                    &mut hasher,
                    &columns[column].to_be_bytes(),
                );
                hasher.finalize_into_reset(hash);
            }
        });
    }

    fn merklize_column_hashes(depth: usize, hashes: &mut [Output<Keccak256>]) {
        let mut offset = 0;
        for width in (1..=depth).rev().map(|depth| 1 << depth) {
            let (input, output) = hashes[offset..].split_at_mut(width);

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
