use ark_std::{collections::BTreeSet, marker::PhantomData, vec::Vec};
use sha3::{digest::Output, Keccak256};

use super::utils::MerkleTree;
use crate::{
    traits::Integer,
    zip::code::{LinearCode, ZipLinearCode, ZipLinearCodeSpec},
};

// N is the width of elements in witness/ polynomial evaluations on hypercube
// L is the width of elements in the encoding matrices
// K is the width of elements in the code
// M is the width of elements in linear combination of code rows
#[derive(Debug, Clone)]
pub struct MultilinearZip<N: Integer, L: Integer, K: Integer, M: Integer>(
    PhantomData<(N, L, K, M)>,
);

/// Parameters for the Zip PCS.
///
/// # Type Parameters
/// - `N`: Width of elements in witness/polynomial evaluations on hypercube.
/// - `L`: Width of elements in the encoding matrices.
#[derive(Clone, Debug)]
pub struct MultilinearZipParams<N: Integer, L: Integer> {
    pub num_vars: usize,
    pub num_rows: usize,
    pub linear_code: ZipLinearCode<N, L>,
}

/// Representantation of a zip commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearZipData<K: Integer> {
    /// The encoded rows of the polynomial matrix representation, referred to as "u-hat" in the Zinc paper
    pub rows: Vec<K>,
    /// Merkle trees of each row
    pub rows_merkle_trees: Vec<MerkleTree>,
}

/// Representantation of a zip commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearZipCommitment {
    /// Roots of the merkle tree of each row
    pub roots: Vec<Output<Keccak256>>,
}

impl<K: Integer> MultilinearZipData<K> {
    pub fn new(rows: Vec<K>, rows_merkle_trees: Vec<MerkleTree>) -> MultilinearZipData<K> {
        MultilinearZipData {
            rows,
            rows_merkle_trees,
        }
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

impl<I: Integer, L: Integer, K: Integer, M: Integer> MultilinearZip<I, L, K, M>
where
    ZipLinearCode<I, L>: LinearCode<I, M>,
{
    pub fn setup<S: ZipLinearCodeSpec, T: ZipTranscript<L>>(
        poly_size: usize,
        transcript: &mut T,
    ) -> MultilinearZipParams<I, L> {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        let linear_code = ZipLinearCode::new_multilinear::<S, T>(
            num_vars,
            20.min((1 << num_vars) - 1),
            transcript,
        );
        let num_rows = ((1 << num_vars) / linear_code.row_len()).next_power_of_two();

        MultilinearZipParams {
            num_vars,
            num_rows,
            linear_code,
        }
    }
}
