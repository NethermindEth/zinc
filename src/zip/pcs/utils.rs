use crate::field_config::ConfigRef;
use ark_ff::Zero;
use ark_std::format;
use ark_std::iterable::Iterable;
use ark_std::vec;
use ark_std::vec::Vec;

use crypto_bigint::Int;
use sha3::{digest::Output, Digest, Keccak256};

use super::error::MerkleError;
use super::structs::MultilinearZipData;
use crate::field::RandomField;
use crate::primitives::Signed;
use crate::{
    field::RandomField as F,
    poly_f::mle::DenseMultilinearExtension as MLE_F,
    poly_z::mle::DenseMultilinearExtension as MLE_Z,
    sumcheck::utils::build_eq_x_r as build_eq_x_r_f,
    zip::{
        pcs_transcript::PcsTranscript,
        utils::{div_ceil, num_threads, parallelize, parallelize_iter},
        Error,
    },
};

fn err_too_many_variates(function: &str, upto: usize, got: usize) -> Error {
    Error::InvalidPcsParam(
        format!(
            "Too many variates of poly to {function} (param supports variates up to {upto} but got {got})"
        )
    )
}

// Ensures that polynomials and evaluation points are of appropriate size
pub(super) fn validate_input<'cfg, 'a, const I: usize, const N: usize>(
    function: &str,
    param_num_vars: usize,
    polys: impl Iterable<Item = &'a MLE_Z<I>>,
    points: impl Iterable<Item = &'a [F<'cfg, N>]>,
) -> Result<(), Error>
where
    'cfg: 'a,
{
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

// Macro to implement ToBytes for integer types
// macro_rules! impl_to_bytes {
//     ($($t:ty),*) => {
//         $(
//             impl ToBytes for $t {
//                 fn to_bytes(&self) -> Vec<u8> {
//                     self.to_be_bytes().to_vec()
//                 }
//             }
//         )*
//     }
// }

impl<T: Signed> ToBytes for T {
    fn to_bytes(&self) -> Vec<u8> {
        <Self as Signed>::to_be_bytes(self)
    }
}

// Manual impl for generic type
impl<const N: usize> ToBytes for Int<N> {
    fn to_bytes(&self) -> Vec<u8> {
        self.to_words()
            .iter()
            .flat_map(|word| word.to_be_bytes())
            .collect()
    }
}
/// A merkle tree in which its layers are concatenated together in a single vector
#[derive(Clone, Debug, Default)]
pub struct MerkleTree {
    pub root: Output<Keccak256>,
    pub depth: usize,
    pub layers: Vec<Output<Keccak256>>,
}

impl MerkleTree {
    pub fn new<T: ToBytes + Send + Sync>(depth: usize, leaves: &[T]) -> Self {
        assert!(leaves.len().is_power_of_two());
        assert_eq!(leaves.len(), 1 << depth);
        let mut layers = vec![Output::<Keccak256>::default(); (2 << depth) - 1];
        Self::compute_leaves_hashes(&mut layers[..leaves.len()], leaves);
        Self::merklize_leaves_hashes(depth, &mut layers);
        Self {
            root: layers.pop().unwrap(),
            depth,
            layers,
        }
    }

    fn compute_leaves_hashes<T: ToBytes + Send + Sync>(
        hashes: &mut [Output<Keccak256>],
        leaves: &[T],
    ) {
        parallelize(hashes, |(hashes, start)| {
            let mut hasher = Keccak256::new();
            for (hash, row) in hashes.iter_mut().zip(start..) {
                let bytes = leaves[row].to_bytes();
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, &bytes);
                hasher.finalize_into_reset(hash);
            }
        });
    }

    fn merklize_leaves_hashes(depth: usize, hashes: &mut [Output<Keccak256>]) {
        assert_eq!(hashes.len(), (2 << depth) - 1);
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
                        hasher.update(input[0]);
                        hasher.update(input[1]);
                        hasher.finalize_into_reset(output);
                    }
                },
            );
            offset += width;
        }
    }
}

#[derive(Clone, Debug)]
pub struct MerkleProof {
    pub merkle_path: Vec<Output<Keccak256>>,
}

impl ark_std::fmt::Display for MerkleProof {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        writeln!(f, "Merkle Path:")?;
        for (i, hash) in self.merkle_path.iter().enumerate() {
            writeln!(f, "Level {}: {:?}", i, hash)?;
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

    pub fn from_vec(vec: Vec<Output<Keccak256>>) -> Self {
        Self { merkle_path: vec }
    }

    pub fn create_proof(merkle_tree: &MerkleTree, leaf: usize) -> Result<Self, MerkleError> {
        let mut offset = 0;
        let path: Vec<Output<Keccak256>> = (1..=merkle_tree.depth)
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
        root: Output<Keccak256>,
        leaf_value: &T,
        leaf_index: usize,
    ) -> Result<(), MerkleError> {
        let mut hasher = Keccak256::new();
        let bytes = leaf_value.to_bytes();
        hasher.update(&bytes);
        let mut current = hasher.finalize_reset();

        let mut index = leaf_index;
        for path_hash in &self.merkle_path {
            if (index & 1) == 0 {
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, &current);
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, path_hash);
            } else {
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, path_hash);
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, &current);
            }

            hasher.finalize_into_reset(&mut current);
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
    pub fn open_at_column<const N: usize, const I: usize, const M: usize>(
        column: usize,
        commit_data: &MultilinearZipData<I, M>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), MerkleError> {
        for row_merkle_tree in commit_data.rows_merkle_trees() {
            let merkle_path = MerkleProof::create_proof(row_merkle_tree, column)?;
            transcript
                .write_merkle_proof(&merkle_path)
                .map_err(|_| MerkleError::FailedMerkleProofWriting)?;
        }
        Ok(())
    }

    pub fn verify_column<const N: usize, T: ToBytes>(
        rows_roots: &[Output<Keccak256>],
        column: &[T],
        column_index: usize,
        transcript: &mut PcsTranscript<N>,
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
pub(super) fn point_to_tensor<'cfg, const N: usize>(
    num_rows: usize,
    point: &[F<'cfg, N>],
    config: ConfigRef<'cfg, N>,
) -> Result<(Vec<F<'cfg, N>>, Vec<F<'cfg, N>>), Error> {
    assert!(num_rows.is_power_of_two());
    let (hi, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let q_0 = if !lo.is_empty() {
        build_eq_x_r_f(lo, config).unwrap()
    } else {
        MLE_F::<RandomField<N>, ConfigRef<N>>::zero()
    };

    let q_1 = if !hi.is_empty() {
        build_eq_x_r_f(hi, config).unwrap()
    } else {
        MLE_F::<RandomField<N>, ConfigRef<N>>::zero()
    };

    Ok((q_0.evaluations, q_1.evaluations))
}

/// For a polynomial arranged in matrix form, this splits the evaluation point into
/// two vectors, `q_0` multiplying on the left and `q_1` multiplying on the right
/// and returns the left vector only
pub(super) fn left_point_to_tensor<'cfg, const N: usize>(
    num_rows: usize,
    point: &[F<'cfg, N>],
    config: ConfigRef<'cfg, N>,
) -> Result<Vec<F<'cfg, N>>, Error> {
    let (_, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let q_0 = if !lo.is_empty() {
        build_eq_x_r_f(lo, config).unwrap()
    } else {
        MLE_F::<RandomField<N>, ConfigRef<N>>::zero()
    };
    Ok(q_0.evaluations)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zip::utils::combine_rows;
    use crypto_bigint::Random;

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
        let leaves_data: Vec<Int<N>> = (0..leaves_len).map(|_| Int::random(&mut rng)).collect();

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
