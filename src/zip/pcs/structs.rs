use ark_std::vec::Vec;
use ark_std::{collections::BTreeSet, marker::PhantomData};

use crypto_bigint::Int;

use sha3::{digest::Output, Keccak256};

use crate::{
    poly_z::mle::DenseMultilinearExtension as DenseMultilinearExtensionZ,
    zip::code::{LinearCodes, Zip, ZipSpec},
};

use super::utils::MerkleTree;

// N is the width of elements in witness/ polynomial evaluations on hypercube
// L is the width of elements in the encoding matrices
// K is the width of elements in the code
// M is the width of elements in linear combination of code rows
#[derive(Debug)]
pub struct MultilinearZip<
    const N: usize,
    const L: usize,
    const K: usize,
    const M: usize,
    S: ZipSpec,
    T: ZipTranscript<L>,
>(PhantomData<S>, PhantomData<T>);

impl<
        const N: usize,
        const L: usize,
        const K: usize,
        const M: usize,
        S: ZipSpec,
        T: ZipTranscript<L>,
    > Clone for MultilinearZip<N, L, K, M, S, T>
{
    fn clone(&self) -> Self {
        Self(PhantomData, PhantomData)
    }
}

#[derive(Clone, Debug)]
pub struct MultilinearZipParams<const N: usize, const L: usize> {
    num_vars: usize,
    num_rows: usize,
    zip: Zip<N, L>,
}

impl<const N: usize, const L: usize> MultilinearZipParams<N, L> {
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    pub fn zip(&self) -> &Zip<N, L> {
        &self.zip
    }
}

/// Representantation of a zip commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearZipData<const N: usize, const K: usize> {
    /// The encoded rows of the polynomial matrix representation
    rows: Vec<Int<K>>,
    /// Merkle trees of each row
    rows_merkle_trees: Vec<MerkleTree>,
}
/// Representantation of a zip commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearZipCommitment<const N: usize> {
    /// Roots of the merkle tree of each row
    roots: Vec<Output<Keccak256>>,
}
impl<const N: usize> MultilinearZipCommitment<N> {
    pub fn new(roots: Vec<Output<Keccak256>>) -> MultilinearZipCommitment<N> {
        MultilinearZipCommitment { roots }
    }
    pub fn roots(&self) -> &[Output<Keccak256>] {
        &self.roots
    }
}

impl<const N: usize, const K: usize> MultilinearZipData<N, K> {
    pub fn new(rows: Vec<Int<K>>, rows_merkle_trees: Vec<MerkleTree>) -> MultilinearZipData<N, K> {
        MultilinearZipData {
            rows,
            rows_merkle_trees,
        }
    }

    pub fn rows(&self) -> &[Int<K>] {
        &self.rows
    }

    pub fn rows_merkle_trees(&self) -> &[MerkleTree] {
        &self.rows_merkle_trees
    }

    pub fn roots(&self) -> Vec<Output<Keccak256>> {
        self.rows_merkle_trees
            .iter()
            .map(|tree| tree.root)
            .collect::<Vec<_>>()
    }

    pub fn root_at_index(&self, index: usize) -> Output<Keccak256> {
        self.rows_merkle_trees[index].root
    }
}

pub trait ZipTranscript<const L: usize> {
    fn get_encoding_element(&mut self) -> Int<L>;
    fn sample_unique_columns(
        &mut self,
        range: ark_std::ops::Range<usize>,
        columns: &mut BTreeSet<usize>,
        count: usize,
    ) -> usize;
}
impl<const I: usize, const L: usize, const K: usize, const M: usize, S, T>
    MultilinearZip<I, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
{
    pub type Param = MultilinearZipParams<I, L>;
    pub type ProverParam = MultilinearZipParams<I, L>;
    pub type VerifierParam = MultilinearZipParams<I, L>;
    pub type Polynomial = DenseMultilinearExtensionZ<I>;
    pub type Data = MultilinearZipData<I, K>;
    pub type Commitment = MultilinearZipCommitment<I>;
    pub type CommitmentChunk = Output<Keccak256>;

    pub fn setup(poly_size: usize, transcript: &mut T) -> Self::Param {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        let zip = Zip::new_multilinear::<S, T>(num_vars, 20.min((1 << num_vars) - 1), transcript);

        MultilinearZipParams {
            num_vars,
            num_rows: ((1 << num_vars) / <Zip<I, L> as LinearCodes<I, M>>::row_len(&zip))
                .next_power_of_two(),
            zip,
        }
    }
}
