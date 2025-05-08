//! Module that defines main structures for PCS as used in the zinc protocol
use std::{collections::BTreeSet, marker::PhantomData};

use crypto_bigint::Int;

use sha3::{digest::Output, Keccak256};

use crate::{
    poly_z::mle::DenseMultilinearExtension,
    zip::code::{LinearCodes, Zip, ZipSpec},
};

use super::utils::MerkleTree;
/// The main Zip PCS struct
/// N is the width of elements in witness/ polynomial evaluations on hypercube
/// L is the width of elements in the encoding matrices
/// K is the width of elements in the code
/// M is the width of elements in linear combination of code rows
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
impl<const N: usize, const L: usize, const K: usize, const M: usize, S, T>
    MultilinearZip<N, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
{
    /// Zip Paramaters
    pub type Param = MultilinearZipParams<N, L>;
    /// The kind of polyomial we can commit to
    pub type Polynomial = DenseMultilinearExtension<N>;
    /// The data that the prover stores in state
    pub type Data = MultilinearZipData<N, K>;
    /// The commitment that is given over to the prover
    pub type Commitment = MultilinearZipCommitment<N>;

    /// Gets the Zip parameters
    pub fn setup(poly_size: usize, transcript: &mut T) -> Self::Param {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        let zip = Zip::new_multilinear::<S, T>(num_vars, 20.min((1 << num_vars) - 1), transcript);

        MultilinearZipParams {
            num_vars,
            num_rows: ((1 << num_vars) / <Zip<N, L> as LinearCodes<N, M>>::row_len(&zip))
                .next_power_of_two(),
            zip,
        }
    }
}

/// Defines some parameters for the ML zip scheme
/// includes the number of variables in polynomials,
/// the number of rows we split the evaluations matrix into,
/// and the encoding scheme.
#[derive(Clone, Debug)]
pub struct MultilinearZipParams<const N: usize, const L: usize> {
    num_vars: usize,
    num_rows: usize,
    zip: Zip<N, L>,
}

impl<const N: usize, const L: usize> MultilinearZipParams<N, L> {
    pub(crate) fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub(crate) fn num_rows(&self) -> usize {
        self.num_rows
    }

    pub(crate) fn zip(&self) -> &Zip<N, L> {
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
    pub(crate) fn new(roots: Vec<Output<Keccak256>>) -> MultilinearZipCommitment<N> {
        MultilinearZipCommitment { roots }
    }
    pub(crate) fn roots(&self) -> &[Output<Keccak256>] {
        &self.roots
    }
}

impl<const N: usize, const K: usize> MultilinearZipData<N, K> {
    pub(crate) fn new(
        rows: Vec<Int<K>>,
        rows_merkle_trees: Vec<MerkleTree>,
    ) -> MultilinearZipData<N, K> {
        MultilinearZipData {
            rows,
            rows_merkle_trees,
        }
    }

    pub(crate) fn rows(&self) -> &[Int<K>] {
        &self.rows
    }

    pub(crate) fn rows_merkle_trees(&self) -> &[MerkleTree] {
        &self.rows_merkle_trees
    }
}

/// Some functionality we need a Fiat Shamir transcript to implement
/// so that we can use it in the PCS. i.e. the PCS can be non-interactive
/// if the prover and verifier have access to the following pseudorandom challenges.
pub trait ZipTranscript<const L: usize> {
    /// Get an element of an encoding matrix for the linear code
    fn get_encoding_element(&mut self) -> Int<L>;
    /// Decide which columns we should challenge
    fn sample_unique_columns(
        &mut self,
        range: std::ops::Range<usize>,
        columns: &mut BTreeSet<usize>,
        count: usize,
    ) -> usize;
}
