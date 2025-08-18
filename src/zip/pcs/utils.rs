use ark_ff::Zero;
use ark_std::{
    cmp::Ordering,
    fmt::{Display, Formatter},
    format,
    io::Write,
    iterable::Iterable,
    vec::Vec,
};
use itertools::Itertools;
use p3_commit::{BatchOpeningRef, Mmcs};
use p3_field::Packable;
use p3_matrix::{Dimensions, dense::RowMajorMatrix};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{CryptographicHasher, PseudoCompressionFunction};

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
#[repr(transparent)]
pub struct MtHash(pub(crate) [u8; blake3::OUT_LEN]);

impl PartialOrd for MtHash {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MtHash {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl Default for MtHash {
    fn default() -> Self {
        MtHash([0; blake3::OUT_LEN])
    }
}

impl AsRef<[u8]> for MtHash {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl Display for MtHash {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let blake3_hash: blake3::Hash = self.0.into();
        <blake3::Hash as Display>::fmt(&blake3_hash, f)
    }
}

#[derive(Debug, Default, Clone)]
pub struct MtHasher;

impl<T: ToBytes + Clone> CryptographicHasher<T, [u8; 32]> for MtHasher {
    fn hash_iter<I>(&self, input: I) -> [u8; 32]
    where
        I: IntoIterator<Item = T>,
    {
        let mut hasher = blake3::Hasher::new();
        for item in input {
            hasher
                .write_all(&item.to_bytes())
                .expect("Failed to write to hasher");
        }
        let hash = hasher.finalize();
        hash.into()
    }
}

#[derive(Debug, Default, Clone)]
pub struct MtPerm;

impl PseudoCompressionFunction<[u8; 32], 2> for MtPerm {
    fn compress(&self, input: [[u8; 32]; 2]) -> [u8; 32] {
        let mut hasher = blake3::Hasher::new();
        for ref item in input {
            hasher.write_all(item).expect("Failed to write to hasher");
        }
        let hash = hasher.finalize();
        hash.into()
    }
}

// cannot reference blake3::OUT_LEN directly
type Matrix<T> = RowMajorMatrix<T>;
type MtMmcs<T> = MerkleTreeMmcs<T, u8, MtHasher, MtPerm, 32>;
type P3MerkleTree<T> = p3_merkle_tree::MerkleTree<T, u8, Matrix<T>, 32>;

#[derive(Debug, Default)]
pub struct MerkleTree<T>
where
    T: Packable + ToBytes + Clone + Send + Sync,
{
    inner: Option<MerkleTreeInner<T>>,
}

#[derive(Debug)]
struct MerkleTreeInner<T> {
    prover_data: P3MerkleTree<T>,
    num_leaves: usize,
}

impl<T> MerkleTree<T>
where
    T: Packable + ToBytes + Clone + Send + Sync,
{
    pub fn new(leaves: Vec<T>) -> Self {
        assert!(leaves.len().is_power_of_two());
        let depth = leaves.len().ilog2() as usize;
        assert_eq!(leaves.len(), 1 << depth);

        let num_leaves = leaves.len();
        let matrix = RowMajorMatrix::new_col(leaves);

        let prover_data = P3MerkleTree::new::<T, _, _, _>(&MtHasher, &MtPerm, vec![matrix]);

        Self {
            inner: Some(MerkleTreeInner {
                prover_data,
                num_leaves,
            }),
        }
    }

    pub fn root(&self) -> MtHash {
        // TODO: Avoid cloning
        MtHash(
            *self
                .inner
                .as_ref()
                .expect("Merkle tree not initialized")
                .prover_data
                .root()
                .as_ref(),
        )
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MerkleProof {
    path: Option<Vec<MtHash>>,
    pub num_leaves: usize,
}

impl MerkleProof {
    pub fn new(path: Vec<MtHash>, num_leaves: usize) -> Self {
        MerkleProof {
            path: Some(path),
            num_leaves,
        }
    }

    pub fn create_proof<T>(merkle_tree: &MerkleTree<T>, leaf: usize) -> Result<Self, MerkleError>
    where
        T: Packable + ToBytes + Clone,
    {
        let mt = merkle_tree
            .inner
            .as_ref()
            .ok_or(MerkleError::InvalidRootHash)?;
        let prover = MtMmcs::<T>::new(MtHasher, MtPerm);
        let opening = prover.open_batch(leaf, &mt.prover_data);
        let path = opening.opening_proof.into_iter().map(MtHash).collect();
        Ok(Self::new(path, mt.num_leaves))
    }

    pub fn path(&self) -> Option<&[MtHash]> {
        self.path.as_deref()
    }

    pub fn verify<T>(
        &self,
        root: &MtHash,
        leaf_value: &T,
        leaf_index: usize,
    ) -> Result<(), MerkleError>
    where
        T: Packable + ToBytes + Clone,
    {
        let Some(path) = self.path.as_ref() else {
            return Err(MerkleError::InvalidMerkleProof(
                "Merkle proof is None".to_string(),
            ));
        };
        let prover = MtMmcs::<T>::new(MtHasher, MtPerm);

        let values = vec![vec![*leaf_value]];
        let proof = path.iter().map(|h| h.0).collect_vec();
        let proof = BatchOpeningRef::new(&values, &proof);
        prover
            .verify_batch(
                &root.0.into(),
                &[Dimensions {
                    width: 0,
                    height: self.num_leaves,
                }],
                leaf_index,
                proof,
            )
            .map_err(|e| {
                MerkleError::InvalidMerkleProof(format!("Failed to validate Merkle proof: {:?}", e))
            })
    }
}

impl Display for MerkleProof {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        match &self.path {
            None => write!(f, "Merkle Proof: None")?,
            Some(path) => {
                writeln!(f, "Merkle Path: {}", path.iter().join(", "))?;
                writeln!(f, "# Leaves: {}", self.num_leaves)?;
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

    pub fn verify_column<F: Field, T: Packable + ToBytes + Clone>(
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

        let merkle_tree = MerkleTree::new(leaves_data.clone());

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
