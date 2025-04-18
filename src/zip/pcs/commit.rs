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
        let codeword_len = pp.zip().codeword_len();
        // We deviate from the paper to merkleize each column instead of each row
        let merkle_depth = codeword_len.next_power_of_two().ilog2() as usize;

        let mut hashes =
            vec![Output::<Keccak256>::default(); codeword_len * ((2 << merkle_depth) - 1)];
        let rows = Self::encode_rows(pp, codeword_len, row_len, poly);
        let mut temp_hashes = vec![Output::<Keccak256>::default(); row_len * codeword_len];
        Self::compute_rows_hashes(&mut temp_hashes, &rows);

        // Process each column's hashes
        for i in 0..codeword_len {
            let merkle_tree_size = (2 << merkle_depth) - 1;
            let start_idx = i * merkle_tree_size;

            // Copy hashes for this column
            for j in 0..row_len {
                hashes[start_idx + j] = temp_hashes[i * row_len + j];
            }

            // Merklize this column's hashes
            let end_idx = start_idx + merkle_tree_size;
            Self::merklize_rows_hashes(merkle_depth, &mut hashes[start_idx..end_idx]);
        }

        // Split hashes into chunks of size (2 << merkle_depth) - 1
        let mut split_hashes = Vec::with_capacity(codeword_len);
        let chunk_size = (2 << merkle_depth) - 1;
        let mut roots = Vec::with_capacity(codeword_len);

        for chunk in hashes.chunks(chunk_size) {
            let mut chunk_vec = chunk.to_vec();
            let root = chunk_vec.pop().unwrap();
            roots.push(root);
            split_hashes.push(chunk_vec);
        }

        Ok((
            MultilinearZipData::new(rows, split_hashes, roots.clone()),
            MultilinearZipCommitment::new(roots),
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
        // assert_eq!(pp.num_rows(), poly.evaluations.len().isqrt());
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

    fn compute_rows_hashes(hashes: &mut [Output<Keccak256>], rows: &[I512]) {
        // TODO:improve this without transposing the rows into columns
        parallelize(hashes, |(hashes, start)| {
            let mut hasher = Keccak256::new();
            for (hash, row) in hashes.iter_mut().zip(start..) {
                // For each row, iterate through all columns at that row position
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, &rows[row].to_be_bytes());
                hasher.finalize_into_reset(hash);
            }
        });
    }

    fn merklize_rows_hashes(depth: usize, hashes: &mut [Output<Keccak256>]) {
        let mut offset = 0;
        for width in (1..=depth).rev().map(|depth| 1 << depth) {
            let (current_layer, next_layer) = hashes[offset..].split_at_mut(width);

            let chunk_size = div_ceil(next_layer.len(), num_threads());
            parallelize_iter(
                current_layer
                    .chunks(2 * chunk_size)
                    .zip(next_layer.chunks_mut(chunk_size)),
                |(input, output)| {
                    let mut hasher = Keccak256::new();
                    for (input, output) in input.chunks_exact(2).zip(output.iter_mut()) {
                        hasher.update(&input[0]);
                        hasher.update(&input[1]);
                        hasher.finalize_into_reset(output);
                    }
                },
            );
            offset += width;
        }
    }
}
