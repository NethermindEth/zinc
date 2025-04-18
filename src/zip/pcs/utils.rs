use ark_ff::Zero;
use ark_std::iterable::Iterable;
use i256::I512;
use sha3::{digest::Output, Digest, Keccak256};

use crate::{
    field::RandomField as F,
    field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension as MLE_F,
    poly_z::mle::{build_eq_x_r as build_eq_x_r_z, DenseMultilinearExtension as MLE_Z},
    sumcheck::utils::build_eq_x_r as build_eq_x_r_f,
    zip::{pcs_transcript::PcsTranscript, Error},
};

use super::structs::MultilinearZipData;

fn err_too_many_variates(function: &str, upto: usize, got: usize) -> Error {
    Error::InvalidPcsParam(
        format!(
            "Too many variates of poly to {function} (param supports variates up to {upto} but got {got})"
        )
    )
}

// Ensures that polynomials and evaluation points are of appropriate size
pub(super) fn validate_input<'a>(
    function: &str,
    param_num_vars: usize,
    polys: impl Iterable<Item = &'a MLE_Z>,
    points: impl Iterable<Item = &'a Vec<i64>>,
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

#[derive(Clone)]
pub struct MerkleProof {
    pub merkle_path: Vec<Output<Keccak256>>,
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

    pub fn create_proof(
        row_hashes: &[Output<Keccak256>],
        leaf: usize,
        merkle_depth: usize,
    ) -> Self {
        let mut offset = 0;
        let path = (1..=merkle_depth)
            .rev()
            .map(|depth| {
                let width = 1 << depth;
                let idx = (leaf >> (merkle_depth - depth)) ^ 1;
                let hash = row_hashes[offset + idx];
                offset += width;
                hash
            })
            .collect();
        MerkleProof::from_vec(path)
    }

    pub fn verify(&self, root: Output<Keccak256>, leaf_value: I512) -> bool {
        let mut hasher = Keccak256::new();
        hasher.update(&leaf_value.to_be_bytes());
        let mut current = hasher.finalize_reset();

        for i in 0..self.merkle_path.len() {
            <Keccak256 as sha3::digest::Update>::update(&mut hasher, &current);
            <Keccak256 as sha3::digest::Update>::update(&mut hasher, &self.merkle_path[i]);
            hasher.finalize_into_reset(&mut current);
        }
        current == root
    }
}

/// This is a helper struct to open a column in a multilinear polynomial
/// Opening a column `j` in an `n x m` matrix `u_hat` requires opening `m` Merkle trees,
/// one for each row at position j
/// Note that the proof is written to the transcript and the order of the proofs is the same as the order of the columns
#[derive(Clone)]
pub struct ColumnOpening {}

impl ColumnOpening {
    pub fn open_at_column<const N: usize>(
        merkle_depth: usize,
        column: usize,
        comm: &MultilinearZipData<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> () {
        for row_hashes in comm.intermediate_rows_hashes() {
            let merkle_path = MerkleProof::create_proof(row_hashes, column, merkle_depth);
            transcript.write_merkle_proof(&merkle_path).unwrap();
        }
    }

    pub fn verify_column<const N: usize>(
        rows_roots: &[Output<Keccak256>],
        column: &[I512],
        transcript: &mut PcsTranscript<N>,
    ) -> bool {
        let mut valid = true;
        for (root, leaf) in rows_roots.iter().zip(column) {
            let proof = transcript.read_merkle_proof().unwrap();
            valid &= proof.verify(*root, *leaf);
        }
        valid
    }
}

pub(super) fn point_to_tensor_z(
    num_rows: usize,
    point: &[i64],
) -> Result<(Vec<i64>, Vec<i64>), Error> {
    assert!(num_rows.is_power_of_two());
    let (hi, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let t_0 = if !lo.is_empty() {
        build_eq_x_r_z(lo).unwrap()
    } else {
        MLE_Z::zero()
    };

    let t_1 = if !hi.is_empty() {
        build_eq_x_r_z(hi).unwrap()
    } else {
        MLE_Z::zero()
    };

    Ok((t_0.evaluations, t_1.evaluations))
}

pub(super) fn point_to_tensor_f<const N: usize>(
    num_rows: usize,
    point: &[F<N>],
    config: *const FieldConfig<N>,
) -> Result<(Vec<F<N>>, Vec<F<N>>), Error> {
    assert!(num_rows.is_power_of_two());
    let (hi, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let t_0 = if !lo.is_empty() {
        build_eq_x_r_f(lo, config).unwrap()
    } else {
        MLE_F::<N>::zero()
    };

    let t_1 = if !hi.is_empty() {
        build_eq_x_r_f(hi, config).unwrap()
    } else {
        MLE_F::<N>::zero()
    };

    Ok((t_0.evaluations, t_1.evaluations))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zip::utils::combine_rows;
    use i256::I512;

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

    fn hash_data(data: &[u8]) -> Output<Keccak256> {
        let mut hasher = Keccak256::new();
        hasher.update(data);
        hasher.finalize()
    }

    fn merklize_hashes(depth: usize, hashes: &mut [Output<Keccak256>]) {
        let mut offset = 0;
        for width in (1..=depth).rev().map(|depth| 1 << depth) {
            let (current_layer, next_layer) = hashes[offset..].split_at_mut(width);
            let mut hasher = Keccak256::new();

            for (input, output) in current_layer.chunks_exact(2).zip(next_layer.iter_mut()) {
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, &input[0]);
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, &input[1]);
                hasher.finalize_into_reset(output);
            }
            offset += width;
        }
    }

    #[test]
    fn test_merkle_proof() {
        // Example data to create a Merkle tree
        let leaves_data = vec![I512::from(1), I512::from(2), I512::from(3), I512::from(4)];

        // Hash the leaves
        let mut temp_hashes = vec![Output::<Keccak256>::default(); 4];
        let mut hasher = Keccak256::new();
        for (i, hash) in temp_hashes.iter_mut().enumerate() {
            // For each row, iterate through all columns at that row position
            <Keccak256 as sha3::digest::Update>::update(&mut hasher, &leaves_data[i].to_be_bytes());
            hasher.finalize_into_reset(hash);
        }

        let mut merkle_tree = vec![Output::<Keccak256>::default(); 2 * leaves_data.len() - 1];
        merkle_tree[..leaves_data.len()].copy_from_slice(&temp_hashes);
        merklize_hashes(2, &mut merkle_tree);
        let root = merkle_tree.pop().unwrap();
        // Create a proof for the first leaf
        let merkle_depth = 2; // Example depth for 4 leaves
        let proof = MerkleProof::create_proof(&merkle_tree, 0, merkle_depth);

        // Verify the proof
        let is_valid = proof.verify(root, leaves_data[0]);

        // Assert the verification result
        assert!(is_valid, "The Merkle proof verification failed");
    }
}
