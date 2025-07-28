use ark_std::{collections::BTreeSet, marker::PhantomData, vec::Vec};
use sha3::{digest::Output, Keccak256};

use super::utils::MerkleTree;
use crate::{
    traits::Integer,
    zip::code::{LinearCodes, Zip, ZipSpec},
};

// N is the width of elements in witness/ polynomial evaluations on hypercube
// L is the width of elements in the encoding matrices
// K is the width of elements in the code
// M is the width of elements in linear combination of code rows
#[derive(Debug, Clone)]
pub struct MultilinearZip<
    N: Integer,
    L: Integer,
    K: Integer,
    M: Integer,
    S: ZipSpec,
    T: ZipTranscript<L>,
>(PhantomData<(N, L, K, M, S, T)>);

#[derive(Clone, Debug)]
pub struct MultilinearZipParams<N: Integer, L: Integer> {
    pub num_vars: usize,
    pub num_rows: usize,
    pub zip: Zip<N, L>,
}

/// Representantation of a zip commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearZipData<K: Integer> {
    /// The encoded rows of the polynomial matrix representation
    rows: Vec<K>,
    /// Merkle trees of each row
    rows_merkle_trees: Vec<MerkleTree>,
}
/// Representantation of a zip commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearZipCommitment<I: Integer> {
    /// Roots of the merkle tree of each row
    roots: Vec<Output<Keccak256>>,
    phantom: PhantomData<I>,
}
impl<I: Integer> MultilinearZipCommitment<I> {
    pub fn new(roots: Vec<Output<Keccak256>>) -> MultilinearZipCommitment<I> {
        MultilinearZipCommitment {
            roots,
            phantom: PhantomData,
        }
    }
    pub fn roots(&self) -> &[Output<Keccak256>] {
        &self.roots
    }
}

impl<K: Integer> MultilinearZipData<K> {
    pub fn new(rows: Vec<K>, rows_merkle_trees: Vec<MerkleTree>) -> MultilinearZipData<K> {
        MultilinearZipData {
            rows,
            rows_merkle_trees,
        }
    }

    pub fn rows(&self) -> &[K] {
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

pub trait ZipTranscript<I: Integer> {
    fn get_encoding_element(&mut self) -> I;
    fn sample_unique_columns(
        &mut self,
        range: ark_std::ops::Range<usize>,
        columns: &mut BTreeSet<usize>,
        count: usize,
    ) -> usize;
}
impl<I: Integer, L: Integer, K: Integer, M: Integer, S, T> MultilinearZip<I, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
    Zip<I, L>: LinearCodes<I, M>,
{
    pub fn setup(poly_size: usize, transcript: &mut T) -> MultilinearZipParams<I, L> {
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
