use ark_ff::Zero;
use ark_std::{
    cfg_chunks, cfg_chunks_mut, cfg_iter_mut, format, iterable::Iterable, vec, vec::Vec,
};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::{error::MerkleError, structs::MultilinearZipData};
use crate::{
    poly_f::mle::DenseMultilinearExtension as MLE_F,
    poly_z::mle::DenseMultilinearExtension as MLE_Z,
    sumcheck::utils::build_eq_x_r as build_eq_x_r_f,
    traits::{Field, Integer},
    zip::{
        Error,
        pcs_transcript::PcsTranscript,
        utils::{div_ceil, num_threads},
    },
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

/// A merkle tree in which its layers are concatenated together in a single vector
#[derive(Clone, Debug)]
pub struct MerkleTree {
    pub root: blake3::Hash,
    pub depth: usize,
    pub layers: Vec<blake3::Hash>,
}

impl MerkleTree {
    pub fn new<T: ToBytes + Send + Sync>(depth: usize, leaves: &[T]) -> Self {
        assert!(leaves.len().is_power_of_two());
        assert_eq!(leaves.len(), 1 << depth);
        let mut layers = vec![blake3::Hash::from_bytes([0; blake3::OUT_LEN]); (2 << depth) - 1];
        Self::compute_leaves_hashes(&mut layers[..leaves.len()], leaves);
        Self::merklize_leaves_hashes(depth, &mut layers);
        Self {
            root: layers.pop().unwrap(),
            depth,
            layers,
        }
    }

    fn compute_leaves_hashes<T: ToBytes + Send + Sync>(hashes: &mut [blake3::Hash], leaves: &[T]) {
        cfg_iter_mut!(hashes)
            .enumerate()
            .for_each(|(row, hash)| *hash = blake3::hash(&leaves[row].to_bytes()));
    }

    fn merklize_leaves_hashes(depth: usize, hashes: &mut [blake3::Hash]) {
        assert_eq!(hashes.len(), (2 << depth) - 1);
        let mut offset = 0;
        for width in (1..=depth).rev().map(|depth| 1 << depth) {
            let (current_layer, next_layer) = hashes[offset..].split_at_mut(width);

            let chunk_size = div_ceil(next_layer.len(), num_threads());

            cfg_chunks!(current_layer, 2 * chunk_size) // Use chunks twice the size for the input layer
                .zip(cfg_chunks_mut!(next_layer, chunk_size))
                .for_each(|(input, output)| {
                    let mut hasher = blake3::Hasher::new();
                    for (input, output) in input.chunks_exact(2).zip(output.iter_mut()) {
                        hasher.update(input[0].as_bytes());
                        hasher.update(input[1].as_bytes());
                        *output = hasher.finalize();
                        hasher.reset();
                    }
                });

            offset += width;
        }
    }
}

impl Default for MerkleTree {
    fn default() -> Self {
        MerkleTree {
            root: blake3::Hash::from_bytes([0; blake3::OUT_LEN]),
            depth: 0,
            layers: vec![],
        }
    }
}

#[derive(Clone, Debug)]
pub struct MerkleProof {
    pub merkle_path: Vec<blake3::Hash>,
}

impl ark_std::fmt::Display for MerkleProof {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        writeln!(f, "Merkle Path:")?;
        for (i, hash) in self.merkle_path.iter().enumerate() {
            writeln!(f, "Level {i}: {hash:?}")?;
        }
        Ok(())
    }
}

impl Default for MerkleProof {
    fn default() -> Self {
        Self::new()
    }
}

impl MerkleProof {
    pub fn new() -> Self {
        Self {
            merkle_path: vec![],
        }
    }

    pub fn from_vec(vec: Vec<blake3::Hash>) -> Self {
        Self { merkle_path: vec }
    }

    pub fn create_proof(merkle_tree: &MerkleTree, leaf: usize) -> Result<Self, MerkleError> {
        let mut offset = 0;
        let path: Vec<blake3::Hash> = (1..=merkle_tree.depth)
            .rev()
            .map(|depth| {
                let width = 1 << depth;
                let idx = (leaf >> (merkle_tree.depth - depth)) ^ 1;
                let hash = merkle_tree.layers[offset + idx];
                offset += width;
                hash
            })
            .collect();
        Ok(MerkleProof::from_vec(path))
    }

    pub fn verify<T: ToBytes>(
        &self,
        root: blake3::Hash,
        leaf_value: &T,
        leaf_index: usize,
    ) -> Result<(), MerkleError> {
        let mut hasher = blake3::Hasher::new();
        let bytes = leaf_value.to_bytes();
        hasher.update(&bytes);
        let mut current = hasher.finalize();
        hasher.reset();

        let mut index = leaf_index;
        for path_hash in &self.merkle_path {
            if (index & 1) == 0 {
                hasher.update(current.as_bytes());
                hasher.update(path_hash.as_bytes());
            } else {
                hasher.update(path_hash.as_bytes());
                hasher.update(current.as_bytes());
            }

            current = hasher.finalize();
            hasher.reset();
            index /= 2;
        }
        if current != root {
            return Err(MerkleError::InvalidMerkleProof(
                "Merkle proof verification failed".into(),
            ));
        }
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
        rows_roots: &[blake3::Hash],
        column: &[T],
        column_index: usize,
        transcript: &mut PcsTranscript<F>,
    ) -> Result<(), MerkleError> {
        for (root, leaf) in rows_roots.iter().zip(column) {
            let proof = transcript
                .read_merkle_proof()
                .map_err(|_| MerkleError::FailedMerkleProofReading)?;
            proof.verify(*root, leaf, column_index)?;
        }
        Ok(())
    }
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

        let result = combine_rows(&coeffs, &evaluations, row_len);

        assert_eq!(result, vec![(3 + 2 * 5), (4 + 2 * 6)]);
    }

    #[test]
    fn test_second_combination() {
        let coeffs = vec![3, 4];
        let evaluations = vec![2, 4, 6, 8];
        let row_len = 2;

        let result = combine_rows(&coeffs, &evaluations, row_len);

        assert_eq!(result, vec![(3 * 2 + 4 * 6), (3 * 4 + 4 * 8)]);
    }
    #[test]
    fn test_large_values() {
        let coeffs = vec![1000, -500];
        let evaluations = vec![2000, -3000, 4000, -5000];
        let row_len = 2;

        let result = combine_rows(&coeffs, &evaluations, row_len);

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

        let merkle_depth = leaves_data.len().next_power_of_two().ilog2() as usize;
        let merkle_tree = MerkleTree::new(merkle_depth, &leaves_data);

        // Print tree structure after merklizing
        let root = merkle_tree.root;
        // Create a proof for the first leaf
        for (i, leaf) in leaves_data.iter().enumerate() {
            let proof =
                MerkleProof::create_proof(&merkle_tree, i).expect("Merkle proof creation failed");

            // Verify the proof
            proof
                .verify(root, leaf, i)
                .expect("Merkle proof verification failed");
        }
    }
}
