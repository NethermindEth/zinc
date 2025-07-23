use ark_std::{iterable::Iterable, vec, vec::Vec};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment, MultilinearZipData, ZipTranscript},
    utils::{validate_input, MerkleTree},
};
use crate::{
    poly_z::mle::{
        DenseMultilinearExtension as DenseMultilinearExtensionZ, DenseMultilinearExtension,
    },
    traits::{Field, Integer},
    zip::{
        code::{LinearCodes, Zip, ZipSpec},
        pcs::{structs::MultilinearZipParams, utils::ToBytes},
        utils::{div_ceil, num_threads, parallelize_iter},
        Error,
    },
};

impl<I: Integer, L: Integer, K: Integer, M: Integer, S, T> MultilinearZip<I, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
    M: for<'a> From<&'a I>,
    K: for<'a> From<&'a I> + ToBytes,
    Zip<I, L>: LinearCodes<I, M> + LinearCodes<I, K>,
{
    /// TODO: validate_input method requires a parameter points which is an iterable of type F
    pub fn commit<F: Field>(
        pp: &MultilinearZipParams<I, L>,
        poly: &DenseMultilinearExtensionZ<I>,
    ) -> Result<(MultilinearZipData<I, K>, MultilinearZipCommitment<I>), Error> {
        validate_input::<I, F>("commit", pp.num_vars(), [poly], None)?;

        let row_len = <Zip<I, L> as LinearCodes<I, M>>::row_len(pp.zip());
        let codeword_len = <Zip<I, L> as LinearCodes<I, M>>::codeword_len(pp.zip());
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
            MultilinearZipData::<I, K>::new(rows, rows_merkle_trees),
            MultilinearZipCommitment::new(roots),
        ))
    }
    #[allow(clippy::type_complexity)]
    pub fn batch_commit<'a, F: Field>(
        pp: &MultilinearZipParams<I, L>,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension<I>>,
    ) -> Result<Vec<(MultilinearZipData<I, K>, MultilinearZipCommitment<I>)>, Error>
    where
        I: 'a,
    {
        polys
            .iter()
            .map(|poly| Self::commit::<F>(pp, poly))
            .collect()
    }

    /// Encodes the rows of the polynomial concatenating each encoded row
    pub fn encode_rows(
        pp: &MultilinearZipParams<I, L>,
        codeword_len: usize,
        row_len: usize,
        poly: &DenseMultilinearExtensionZ<I>,
    ) -> Vec<K> {
        // assert_eq!(pp.num_rows(), poly.evaluations.len().isqrt());
        assert_eq!(codeword_len, row_len * 2);
        let rows_per_thread = div_ceil(pp.num_rows(), num_threads());
        let mut encoded_rows = vec![K::default(); pp.num_rows() * codeword_len];

        parallelize_iter(
            encoded_rows
                .chunks_exact_mut(rows_per_thread * codeword_len)
                .zip(poly.evaluations.chunks_exact(rows_per_thread * row_len)),
            |(encoded_chunk, evals)| {
                for (row, evals) in encoded_chunk
                    .chunks_exact_mut(codeword_len)
                    .zip(evals.chunks_exact(row_len))
                {
                    row.clone_from_slice(pp.zip().encode(evals).as_slice());
                }
            },
        );

        encoded_rows
    }
}
