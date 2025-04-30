use std::{collections::BTreeSet, marker::PhantomData};

use i256::I512;
use sha3::{digest::Output, Keccak256};

use crate::{
    poly_z::mle::DenseMultilinearExtension,
    zip::code::{LinearCodes, Zip, ZipSpec},
};

use super::utils::MerkleTree;

#[derive(Debug)]
pub struct MultilinearZip<const N: usize, S: ZipSpec, T: ZipTranscript>(
    PhantomData<S>,
    PhantomData<T>,
);

impl<const N: usize, S: ZipSpec, T: ZipTranscript> Clone for MultilinearZip<N, S, T> {
    fn clone(&self) -> Self {
        Self(PhantomData, PhantomData)
    }
}

#[derive(Clone, Debug)]
pub struct MultilinearZipParams<const N: usize> {
    num_vars: usize,
    num_rows: usize,
    zip: Zip<N>,
}

impl<const N: usize> MultilinearZipParams<N> {
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    pub fn zip(&self) -> &Zip<N> {
        &self.zip
    }
}

/// Representantation of a zip commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearZipData<const N: usize> {
    /// The encoded rows of the polynomial matrix representation
    rows: Vec<I512>,
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

impl<const N: usize> MultilinearZipData<N> {
    pub fn new(rows: Vec<I512>, rows_merkle_trees: Vec<MerkleTree>) -> MultilinearZipData<N> {
        MultilinearZipData {
            rows,
            rows_merkle_trees,
        }
    }

    pub fn rows(&self) -> &[I512] {
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

pub trait ZipTranscript {
    fn get_encoding_element(&mut self) -> i128;
    fn sample_unique_columns(
        &mut self,
        range: std::ops::Range<usize>,
        columns: &mut BTreeSet<usize>,
        count: usize,
    ) -> usize;
}
impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub type Param = MultilinearZipParams<N>;
    pub type ProverParam = MultilinearZipParams<N>;
    pub type VerifierParam = MultilinearZipParams<N>;
    pub type Polynomial = DenseMultilinearExtension<N>;
    pub type Data = MultilinearZipData<N>;
    pub type Commitment = MultilinearZipCommitment<N>;
    pub type CommitmentChunk = Output<Keccak256>;

    pub fn setup(poly_size: usize, transcript: &mut T) -> Self::Param {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        let zip = Zip::new_multilinear::<S, T>(num_vars, 20.min((1 << num_vars) - 1), transcript);

        MultilinearZipParams {
            num_vars,
            num_rows: ((1 << num_vars) / zip.row_len()).next_power_of_two(),
            zip,
        }
    }
}
