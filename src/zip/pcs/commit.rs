use ark_std::iterable::Iterable;
use crypto_bigint::Int;

use crate::{
    poly_z::mle::DenseMultilinearExtension,
    zip::{
        code::{LinearCodes, Zip, ZipSpec},
        utils::{div_ceil, num_threads, parallelize_iter},
        Error,
    },
};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment, MultilinearZipData, ZipTranscript},
    utils::{validate_input, MerkleTree},
};

impl<const N: usize, const L: usize, const K: usize, const M: usize, S, T>
    MultilinearZip<N, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
{
    /// Commits to a polynomial over integers using the Zip, polynomial commitment scheme.
    ///
    /// # Arguments
    ///
    /// * `pp` - A reference to the prover parameters. i.e. number of variables in polynomial,
    ///     and the number of rows we split them into
    /// * `poly` - A reference to the polynomial that is to be committed.
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing a tuple `(data, commitment)`:
    /// * `data` - The encoding matrix, merkle tree and merkle root
    /// * `commitment` - The public commitment to the polynomial.
    ///     i.e just the merkle root.
    ///
    pub fn commit(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
    ) -> Result<(Self::Data, Self::Commitment), Error> {
        validate_input::<N>("commit", pp.num_vars(), [poly], None)?;

        let row_len = <Zip<N, L> as LinearCodes<N, M>>::row_len(pp.zip());
        let codeword_len = <Zip<N, L> as LinearCodes<N, M>>::codeword_len(pp.zip());
        let merkle_depth: usize = codeword_len.next_power_of_two().ilog2() as usize;

        let rows = Self::encode_rows(pp, codeword_len, row_len, poly);

        let rows_merkle_trees = rows
            .chunks_exact(codeword_len)
            .map(|row| MerkleTree::new(merkle_depth, row))
            .collect::<Vec<_>>();

        assert_eq!(rows_merkle_trees.len(), pp.num_rows());

        let roots = rows_merkle_trees
            .iter()
            .map(|tree| tree.root)
            .collect::<Vec<_>>();

        Ok((
            MultilinearZipData::<N, K>::new(rows, rows_merkle_trees),
            MultilinearZipCommitment::new(roots),
        ))
    }
    /// Allows for comitting to multiple polynomials at the same time
    #[allow(clippy::type_complexity)]
    pub fn batch_commit<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension<N>>,
    ) -> Result<Vec<(Self::Data, Self::Commitment)>, Error> {
        polys.iter().map(|poly| Self::commit(pp, poly)).collect()
    }

    /// Encodes the rows of the polynomial concatenating each encoded row
    fn encode_rows(
        pp: &Self::ProverParam,
        codeword_len: usize,
        row_len: usize,
        poly: &Self::Polynomial,
    ) -> Vec<Int<K>> {
        // assert_eq!(pp.num_rows(), poly.evaluations.len().isqrt());
        assert_eq!(codeword_len, row_len * 2);
        let rows_per_thread = div_ceil(pp.num_rows(), num_threads());
        let mut encoded_rows = vec![Int::<K>::default(); pp.num_rows() * codeword_len];

        parallelize_iter(
            encoded_rows
                .chunks_exact_mut(rows_per_thread * codeword_len)
                .zip(poly.evaluations.chunks_exact(rows_per_thread * row_len)),
            |(encoded_chunk, evals)| {
                for (row, evals) in encoded_chunk
                    .chunks_exact_mut(codeword_len)
                    .zip(evals.chunks_exact(row_len))
                {
                    row.copy_from_slice(pp.zip().encode(evals).as_slice());
                }
            },
        );

        encoded_rows
    }
}
