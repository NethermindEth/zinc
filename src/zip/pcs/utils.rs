use ark_ff::Zero;
use ark_std::{
    cmp::Ordering, fmt::Display, format, hash::Hasher, io::Write, iterable::Iterable, vec::Vec,
};
use itertools::Itertools;
use merkletree::{hash::Algorithm, merkle::Element, proof::Proof, store::VecStore};

use super::{error::MerkleError, structs::MultilinearZipData};
use crate::{
    poly_f::mle::DenseMultilinearExtension as MLE_F,
    poly_z::mle::DenseMultilinearExtension as MLE_Z,
    sumcheck::utils::build_eq_x_r as build_eq_x_r_f,
    traits::{Field, Integer},
    zip::{Error, pcs_transcript::PcsTranscript},
};

fn err_too_many_variates(function: &str, upto: usize, got: usize) -> Error {
    Error::InvalidPcsParam(format!(
        "Too many variates of poly to {function} (param supports variates up to {upto} but got {got})"
    ))
}

// Ensures that polynomials and evaluation points are of appropriate size
pub(super) fn validate_input<'a, I: Integer + 'a, F: Field + 'a>(
    function: &str,
    param_num_vars: usize,
    polys: impl Iterable<Item = &'a MLE_Z<I>>,
    points: impl Iterable<Item = &'a [F]>,
) -> Result<(), Error> {
    // Ensure all the number of variables in the polynomials don't exceed the limit
    for poly in polys.iter() {
        if param_num_vars < poly.num_vars {
            return Err(err_too_many_variates(
                function,
                param_num_vars,
                poly.num_vars,
            ));
        }
    }

    // Ensure all the points are of correct length
    let input_num_vars = polys
        .iter()
        .map(|poly| poly.num_vars)
        .chain(points.iter().map(|point| point.len()))
        .next()
        .expect("To have at least 1 poly or point");

    for point in points.iter() {
        if point.len() != input_num_vars {
            return Err(Error::InvalidPcsParam(format!(
                "Invalid point (expect point to have {input_num_vars} variates but got {})",
                point.len()
            )));
        }
    }
    Ok(())
}

// Define a new trait for converting to bytes
pub trait ToBytes {
    fn to_bytes(&self) -> Vec<u8>;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MtHash(pub blake3::Hash);

impl PartialOrd for MtHash {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MtHash {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.as_bytes().cmp(other.0.as_bytes())
    }
}

impl Default for MtHash {
    fn default() -> Self {
        MtHash(blake3::Hash::from_bytes([0; blake3::OUT_LEN]))
    }
}

impl AsRef<[u8]> for MtHash {
    fn as_ref(&self) -> &[u8] {
        self.0.as_bytes()
    }
}

impl Element for MtHash {
    fn byte_len() -> usize {
        blake3::OUT_LEN
    }

    fn from_slice(bytes: &[u8]) -> Self {
        assert_eq!(bytes.len(), blake3::OUT_LEN);
        MtHash(blake3::Hash::from_slice(bytes).expect("Invalid hash length"))
    }

    fn copy_to_slice(&self, bytes: &mut [u8]) {
        assert_eq!(bytes.len(), blake3::OUT_LEN);
        bytes.copy_from_slice(self.as_ref());
    }
}

#[derive(Debug, Default)]
pub struct MtHasher {
    hasher: blake3::Hasher,
}

impl Hasher for MtHasher {
    fn finish(&self) -> u64 {
        let hashed = self.hasher.finalize();
        hashed
            .as_bytes()
            .chunks(u64::BITS as usize / 8)
            .fold(0_u64, |acc, chunk| {
                acc ^ u64::from_le_bytes(chunk.try_into().unwrap())
            })
    }

    fn write(&mut self, bytes: &[u8]) {
        self.hasher.write_all(bytes).unwrap();
    }
}

impl Algorithm<MtHash> for MtHasher {
    fn hash(&mut self) -> MtHash {
        MtHash(self.hasher.finalize())
    }

    fn reset(&mut self) {
        self.hasher.reset();
    }
}

#[derive(Debug, Default)]
pub struct MerkleTree {
    inner: Option<merkletree::merkle::MerkleTree<MtHash, MtHasher, VecStore<MtHash>>>,
}

impl MerkleTree {
    pub fn new<T: ToBytes + Send + Sync>(leaves: &[T]) -> Self {
        assert!(leaves.len().is_power_of_two());
        let depth = leaves.len().ilog2() as usize;
        assert_eq!(leaves.len(), 1 << depth);

        let leaf_hashes = leaves
            .iter()
            .map(|leaf| MtHash(blake3::hash(&leaf.to_bytes())))
            .collect::<Vec<_>>();
        Self {
            inner: Some(merkletree::merkle::MerkleTree::new(leaf_hashes).unwrap()),
        }
    }

    pub fn root(&self) -> MtHash {
        self.inner
            .as_ref()
            .expect("Merkle tree not initialized")
            .root()
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MerkleProof {
    inner: Option<Proof<MtHash>>,
}

impl MerkleProof {
    pub fn new(proof: Proof<MtHash>) -> Self {
        MerkleProof { inner: Some(proof) }
    }

    pub fn create_proof(merkle_tree: &MerkleTree, leaf: usize) -> Result<Self, MerkleError> {
        let mt = merkle_tree
            .inner
            .as_ref()
            .ok_or(MerkleError::InvalidRootHash)?;
        let proof = mt.gen_proof(leaf).map_err(|e| {
            MerkleError::InvalidMerkleProof(format!("Failed to create Merkle proof: {}", e))
        })?;
        Ok(MerkleProof { inner: Some(proof) })
    }

    pub fn inner(&self) -> Option<&Proof<MtHash>> {
        self.inner.as_ref()
    }

    pub fn verify<T: ToBytes>(
        &self,
        root: &MtHash,
        leaf_value: &T,
        leaf_index: usize,
    ) -> Result<(), MerkleError> {
        let Some(proof) = self.inner.as_ref() else {
            return Err(MerkleError::InvalidMerkleProof(
                "Merkle proof is None".to_string(),
            ));
        };

        let binary_string_path = proof.path().iter().rev().join("");
        let usize_path = usize::from_str_radix(&binary_string_path, 2).map_err(|_| {
            MerkleError::InvalidMerkleProof(format!(
                "Failed to parse binary string path: {binary_string_path}"
            ))
        })?;

        if leaf_index != usize_path {
            return Err(MerkleError::InvalidMerkleProof(format!(
                "Leaf index to path mismatch: expected {leaf_index}, got {binary_string_path} (reversed {usize_path})",
            )));
        }

        if root != &proof.root() {
            return Err(MerkleError::InvalidMerkleProof(format!(
                "Root hash mismatch: expected {}, got {}",
                root.0,
                proof.root().0
            )));
        }

        let leaf_hash = hash_leaf(leaf_value);
        if leaf_hash != proof.item() {
            return Err(MerkleError::InvalidMerkleProof(format!(
                "Leaf hash mismatch: expected {}, got {}",
                leaf_hash.0,
                proof.item().0
            )));
        }

        proof.validate::<MtHasher>().map_err(|e| {
            MerkleError::InvalidMerkleProof(format!("Failed to validate Merkle proof: {}", e))
        })?;
        Ok(())
    }
}

impl Display for MerkleProof {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        match &self.inner {
            None => write!(f, "Merkle Proof: None")?,
            Some(proof) => {
                writeln!(f, "Merkle Path: {}", proof.path().iter().rev().join(""))?;
                writeln!(
                    f,
                    "Merkle Lemmas: {}",
                    proof.lemma().iter().map(|h| h.0).join(", ")
                )?;
            }
        };
        Ok(())
    }
}

/// This is a helper struct to open a column in a multilinear polynomial
/// Opening a column `j` in an `n x m` matrix `u_hat` requires opening `m` Merkle trees,
/// one for each row at position j
/// Note that the proof is written to the transcript and the order of the proofs is the same as the order of the columns
#[derive(Clone)]
pub struct ColumnOpening {}

impl ColumnOpening {
    pub fn open_at_column<F: Field, M: Integer>(
        column: usize,
        commit_data: &MultilinearZipData<M>,
        transcript: &mut PcsTranscript<F>,
    ) -> Result<(), MerkleError> {
        for row_merkle_tree in commit_data.rows_merkle_trees.iter() {
            let merkle_path = MerkleProof::create_proof(row_merkle_tree, column)?;
            transcript
                .write_merkle_proof(&merkle_path)
                .map_err(|_| MerkleError::FailedMerkleProofWriting)?;
        }
        Ok(())
    }

    pub fn verify_column<F: Field, T: ToBytes>(
        rows_roots: &[MtHash],
        column: &[T],
        column_index: usize,
        transcript: &mut PcsTranscript<F>,
    ) -> Result<(), MerkleError> {
        for (root, leaf) in rows_roots.iter().zip(column) {
            let proof = transcript
                .read_merkle_proof()
                .map_err(|_| MerkleError::FailedMerkleProofReading)?;
            proof.verify(root, leaf, column_index)?;
        }
        Ok(())
    }
}

fn hash_leaf<T: ToBytes>(data: &T) -> MtHash {
    let mut hasher = MtHasher::default();
    hasher.write(&data.to_bytes());
    let hash = hasher.hash();
    hasher.reset();
    hasher.leaf(hash)
}

/// For a polynomial arranged in matrix form, this splits the evaluation point into
/// two vectors, `q_0` multiplying on the left and `q_1` multiplying on the right
pub(super) fn point_to_tensor<F: Field>(
    num_rows: usize,
    point: &[F],
    config: F::R,
) -> Result<(Vec<F>, Vec<F>), Error> {
    assert!(num_rows.is_power_of_two());
    let (hi, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let q_0 = if !lo.is_empty() {
        build_eq_x_r_f(lo, config).unwrap()
    } else {
        MLE_F::zero()
    };

    let q_1 = if !hi.is_empty() {
        build_eq_x_r_f(hi, config).unwrap()
    } else {
        MLE_F::zero()
    };

    Ok((q_0.evaluations, q_1.evaluations))
}

/// For a polynomial arranged in matrix form, this splits the evaluation point into
/// two vectors, `q_0` multiplying on the left and `q_1` multiplying on the right
/// and returns the left vector only
pub(super) fn left_point_to_tensor<F: Field>(
    num_rows: usize,
    point: &[F],
    config: F::R,
) -> Result<Vec<F>, Error> {
    let (_, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let q_0 = if !lo.is_empty() {
        build_eq_x_r_f(lo, config).unwrap()
    } else {
        MLE_F::<F>::zero()
    };
    Ok(q_0.evaluations)
}

#[cfg(test)]
mod tests {
    use crypto_bigint::Random;

    use super::*;
    use crate::{field::Int, zip::utils::combine_rows};

    #[test]
    fn test_basic_combination() {
        let coeffs = vec![1, 2];
        let evaluations = vec![3, 4, 5, 6];
        let row_len = 2;

        let result = combine_rows(coeffs, evaluations, row_len);

        assert_eq!(result, vec![(3 + 2 * 5), (4 + 2 * 6)]);
    }

    #[test]
    fn test_second_combination() {
        let coeffs = vec![3, 4];
        let evaluations = vec![2, 4, 6, 8];
        let row_len = 2;

        let result = combine_rows(coeffs, evaluations, row_len);

        assert_eq!(result, vec![(3 * 2 + 4 * 6), (3 * 4 + 4 * 8)]);
    }
    #[test]
    fn test_large_values() {
        let coeffs = vec![1000, -500];
        let evaluations = vec![2000, -3000, 4000, -5000];
        let row_len = 2;

        let result = combine_rows(coeffs, evaluations, row_len);

        assert_eq!(
            result,
            vec![
                (1000 * 2000 + (-500) * 4000),
                (1000 * -3000 + (-500) * -5000)
            ]
        );
    }

    #[test]
    fn test_merkle_proof() {
        const N: usize = 3;
        let leaves_len = 1024;
        let mut rng = ark_std::test_rng();
        let leaves_data = (0..leaves_len)
            .map(|_| Int::random(&mut rng))
            .collect::<Vec<Int<N>>>();

        let merkle_tree = MerkleTree::new(&leaves_data);

        // Print tree structure after merklizing
        let root = merkle_tree.root();
        // Create a proof for the first leaf
        for (i, leaf) in leaves_data.iter().enumerate() {
            let proof =
                MerkleProof::create_proof(&merkle_tree, i).expect("Merkle proof creation failed");

            // Verify the proof
            proof
                .verify(&root, leaf, i)
                .expect("Merkle proof verification failed");
        }
    }
}
