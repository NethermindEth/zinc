use ark_std::{vec, vec::Vec};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment, MultilinearZipData},
    utils::{MerkleTree, validate_input},
};
use crate::{
    poly_z::mle::DenseMultilinearExtension,
    traits::{Field, ZipTypes},
    zip::{
        Error,
        code::LinearCode,
        pcs::structs::MultilinearZipParams,
        utils::{div_ceil, num_threads, parallelize_iter},
    },
};

impl<ZT: ZipTypes, LC: LinearCode<ZT>> MultilinearZip<ZT, LC> {
    /// Creates a commitment to a multilinear polynomial using the ZIP PCS scheme.
    ///
    /// This function implements the commitment phase of the ZIP polynomial commitment scheme.
    /// It encodes the polynomial's evaluations using a linear error-correcting code and
    /// then creates Merkle tree commitments to the encoded rows.
    ///
    /// # Algorithm
    /// 1.  Validates that the polynomial's number of variables matches the parameters.
    /// 2.  Arranges the polynomial's evaluations into a matrix with `pp.num_rows` rows.
    /// 3.  Encodes each row using the specified linear code, expanding its length from `row_len` to `codeword_len`.
    /// 4.  Constructs a Merkle tree for each encoded row to produce a root hash.
    /// 5.  Returns the full commitment data (for the prover) and a compact commitment (for the verifier).
    ///
    /// # Type Parameters
    /// - `F`: A generic [`Field`] type used for validation purposes, though not for the commitment itself, which operates on integers.
    ///
    /// # Parameters
    /// - `pp`: Public parameters (`MultilinearZipParams`) containing the configuration for the commitment scheme.
    /// - `poly`: The multilinear polynomial (`DenseMultilinearExtension`) to be committed to.
    ///
    /// # Returns
    /// A `Result` containing a tuple of:
    /// - `MultilinearZipData`: Full data including encoded rows and Merkle trees, kept by the prover for the opening phase.
    /// - `MultilinearZipCommitment`: The compact commitment, consisting of only the Merkle roots, to be sent to the verifier.
    ///
    /// # Errors
    /// - Returns `Error::InvalidPcsParam` if the polynomial has more variables than the parameters support.
    ///
    /// # Panics
    /// - Panics if the number of polynomial evaluations does not perfectly match the expected matrix size (`pp.num_rows * pp.linear_code.row_len()`).
    /// - Panics if the number of generated Merkle trees does not match `pp.num_rows`, indicating an internal logic error.
    pub fn commit<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        poly: &DenseMultilinearExtension<ZT::N>,
    ) -> Result<(MultilinearZipData<ZT::K>, MultilinearZipCommitment), Error> {
        validate_input("commit", pp.num_vars, [poly], None::<&[F]>)?;

        let expected_num_evals = pp.num_rows * pp.linear_code.row_len();
        assert_eq!(
            poly.evaluations.len(),
            expected_num_evals,
            "Polynomial has an incorrect number of evaluations ({}) for the expected matrix size ({})",
            poly.evaluations.len(),
            expected_num_evals
        );

        let row_len = pp.linear_code.row_len();
        let codeword_len = pp.linear_code.codeword_len();
        let merkle_depth: usize = codeword_len.next_power_of_two().ilog2() as usize;

        let rows = Self::encode_rows(pp, codeword_len, row_len, &poly.evaluations);

        let rows_merkle_trees = rows
            .chunks_exact(codeword_len)
            .map(|row| MerkleTree::new(merkle_depth, row))
            .collect::<Vec<_>>();

        assert_eq!(rows_merkle_trees.len(), pp.num_rows);

        let roots = rows_merkle_trees
            .iter()
            .map(|tree| tree.root)
            .collect::<Vec<_>>();

        Ok((
            MultilinearZipData::new(rows, rows_merkle_trees),
            MultilinearZipCommitment { roots },
        ))
    }

    /// Creates a commitment without constructing Merkle trees.
    ///
    /// This function performs the encoding step of the commitment phase but deliberately
    /// skips the computationally intensive step of building Merkle trees. It is intended
    /// **for testing and benchmarking purposes only**, where the full commitment
    /// structure is not required.
    ///
    /// # Parameters
    /// - `pp`: Public parameters (`MultilinearZipParams`).
    /// - `poly`: The multilinear polynomial to commit to.
    ///
    /// # Returns
    /// A `Result` containing `MultilinearZipData` with the encoded rows but empty Merkle trees,
    /// and a `MultilinearZipCommitment` with an empty vector of roots.
    #[allow(dead_code)]
    pub fn commit_no_merkle<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        poly: &DenseMultilinearExtension<ZT::N>,
    ) -> Result<(MultilinearZipData<ZT::K>, MultilinearZipCommitment), Error> {
        validate_input("commit", pp.num_vars, [poly], None::<&[F]>)?;

        let row_len = pp.linear_code.row_len();
        let codeword_len = pp.linear_code.codeword_len();

        let rows = Self::encode_rows(pp, codeword_len, row_len, &poly.evaluations);

        Ok((
            MultilinearZipData::new(rows, vec![]),
            MultilinearZipCommitment { roots: vec![] },
        ))
    }

    /// Commits to a batch of multilinear polynomials.
    ///
    /// This function iterates through a slice of polynomials and calls [`commit`] on each one,
    /// collecting the results into a vector.
    ///
    /// # Parameters
    /// - `pp`: Public parameters (`MultilinearZipParams`) shared across all commitments.
    /// - `polys`: A slice of `DenseMultilinearExtension` polynomials to commit to.
    ///
    /// # Returns
    /// A `Result` containing a `Vec` of commitment tuples, where each tuple corresponds
    /// to a polynomial in the input slice.
    #[allow(clippy::type_complexity)]
    pub fn batch_commit<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        polys: &[DenseMultilinearExtension<ZT::N>],
    ) -> Result<Vec<(MultilinearZipData<ZT::K>, MultilinearZipCommitment)>, Error> {
        polys
            .iter()
            .map(|poly| Self::commit::<F>(pp, poly))
            .collect()
    }

    /// Encodes the evaluations of a polynomial by arranging them into rows and applying a linear code.
    ///
    /// This function treats the polynomial's flat evaluation vector as a matrix with `pp.num_rows`
    /// and encodes each row individually. The resulting encoded rows are concatenated into a
    /// single flat vector. This operation can be parallelized if the `parallel` feature is enabled.
    ///
    /// # Parameters
    /// - `pp`: Public parameters containing matrix dimensions and the linear code.
    /// - `codeword_len`: The length of an encoded row.
    /// - `row_len`: The length of a row before encoding.
    /// - `poly`: The polynomial whose evaluations are to be encoded.
    ///
    /// # Returns
    /// A `Vec<ZT::K>` containing all the encoded rows concatenated together.
    pub fn encode_rows(
        pp: &MultilinearZipParams<ZT, LC>,
        codeword_len: usize,
        row_len: usize,
        evals: &[ZT::N],
    ) -> Vec<ZT::K> {
        let rows_per_thread = div_ceil(pp.num_rows, num_threads());
        let mut encoded_rows = vec![ZT::K::default(); pp.num_rows * codeword_len];

        parallelize_iter(
            encoded_rows
                .chunks_exact_mut(rows_per_thread * codeword_len)
                .zip(evals.chunks_exact(rows_per_thread * row_len)),
            |(encoded_chunk, evals)| {
                for (row, evals) in encoded_chunk
                    .chunks_exact_mut(codeword_len)
                    .zip(evals.chunks_exact(row_len))
                {
                    let encoded: Vec<ZT::K> = pp.linear_code.encode_wide(evals);
                    row.clone_from_slice(encoded.as_slice());
                }
            },
        );

        encoded_rows
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{UniformRand, mem::size_of, slice::from_ref, vec, vec::Vec};
    use crypto_bigint::Random;
    use sha3::{Digest, Keccak256};

    use crate::{
        field::{BigInt, ConfigRef, Int, RandomField},
        field_config,
        poly_z::mle::DenseMultilinearExtension,
        traits::{FieldMap, ZipTypes},
        transcript::KeccakTranscript,
        zip::{
            code::{DefaultLinearCodeSpec, LinearCode, ZipLinearCode},
            code_raa::RaaCode,
            pcs::{
                MerkleTree,
                structs::{MultilinearZip, MultilinearZipParams},
                tests::{MockTranscript, RandomFieldZipTypes},
            },
            pcs_transcript::PcsTranscript,
            utils::div_ceil,
        },
    };

    const INT_LIMBS: usize = 1;
    const FIELD_LIMBS: usize = 4;

    type ZT = RandomFieldZipTypes<INT_LIMBS>;
    type LC = RaaCode<ZT>;
    type F<'cfg> = RandomField<'cfg, FIELD_LIMBS>;
    type TestZip = MultilinearZip<ZT, LC>;

    /// Helper function to set up common parameters for tests.
    fn setup_test_params(
        num_vars: usize,
    ) -> (
        MultilinearZipParams<ZT, ZipLinearCode<ZT>>,
        DenseMultilinearExtension<Int<INT_LIMBS>>,
    ) {
        let poly_size = 1 << num_vars;
        // Correctly calculate num_rows for both even and odd num_vars
        let num_rows = 1 << div_ceil(num_vars, 2);

        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, poly_size, &mut transcript);
        let pp = MultilinearZipParams::new(num_vars, num_rows, code);

        let evaluations: Vec<_> = (1..=poly_size).map(|v| Int::from(v as i32)).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(num_vars, evaluations);

        (pp, poly)
    }

    #[test]
    fn commit_rejects_too_many_variables() {
        let (pp, _) = setup_test_params(3); // Setup for 3 variables

        // Create polynomial with 4 variables (which is > 3)
        let evaluations = (1..=16).map(Int::from).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(4, evaluations);

        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_err());
    }

    #[test]
    fn commit_is_deterministic() {
        let (pp, poly) = setup_test_params(3);

        let result1 = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly).unwrap();
        let result2 = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly).unwrap();

        assert_eq!(result1.1.roots, result2.1.roots);
    }

    #[test]
    fn different_polynomials_produce_different_commitments() {
        let (pp, _) = setup_test_params(3);

        let poly1 = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(1); 8]);
        let poly2 = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(2); 8]);

        let (_, commitment1) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly1).unwrap();
        let (_, commitment2) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly2).unwrap();

        assert_ne!(commitment1.roots, commitment2.roots);
    }

    #[test]
    fn commit_succeeds_for_small_polynomial() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 16, &mut transcript);
        let pp = MultilinearZipParams::new(4, 4, code);

        let evaluations = vec![Int::from(42); 16];
        let poly = DenseMultilinearExtension::from_evaluations_vec(4, evaluations);

        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_ok());
    }

    #[test]
    fn commit_succeeds_for_two_variables() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 4, &mut transcript);
        let pp = MultilinearZipParams::new(2, 2, code);

        let evaluations = vec![Int::from(1), Int::from(2), Int::from(3), Int::from(4)];
        let poly = DenseMultilinearExtension::from_evaluations_vec(2, evaluations);

        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_ok());
    }

    #[test]
    fn merkle_tree_has_correct_depth() {
        let (pp, poly) = setup_test_params(3);
        let (data, _) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly).unwrap();

        let expected_depth = pp.linear_code.codeword_len().next_power_of_two().ilog2() as usize;
        for tree in &data.rows_merkle_trees {
            assert_eq!(tree.depth, expected_depth);
        }
    }

    #[test]
    fn commit_no_merkle_produces_empty_trees() {
        let (pp, poly) = setup_test_params(3);
        let result = MultilinearZip::<ZT, _>::commit_no_merkle::<F>(&pp, &poly);
        assert!(result.is_ok());

        let (data, commitment) = result.unwrap();
        assert_eq!(data.rows.len(), pp.num_rows * pp.linear_code.codeword_len());
        assert!(data.rows_merkle_trees.is_empty());
        assert!(commitment.roots.is_empty());
    }

    #[test]
    fn batch_commit_succeeds() {
        let (pp, _) = setup_test_params(3);

        let polys = vec![
            DenseMultilinearExtension::from_evaluations_vec(3, (1..=8).map(Int::from).collect()),
            DenseMultilinearExtension::from_evaluations_vec(3, (9..=16).map(Int::from).collect()),
        ];

        let results = MultilinearZip::<ZT, _>::batch_commit::<F>(&pp, &polys);
        assert!(results.is_ok());

        let outputs = results.unwrap();
        assert_eq!(outputs.len(), 2);
        assert_ne!(outputs[0].1.roots, outputs[1].1.roots);
    }

    #[test]
    fn encode_rows_produces_correct_size() {
        let (pp, poly) = setup_test_params(3);
        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly.evaluations,
        );

        assert_eq!(encoded.len(), pp.num_rows * pp.linear_code.codeword_len());
    }

    /// Verifies that the output of `encode_rows` is semantically correct by
    /// comparing it to a direct, row-by-row encoding.
    #[test]
    fn encoded_rows_match_linear_code_definition() {
        let (pp, poly) = setup_test_params(3);
        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly.evaluations,
        );

        for (i, row_chunk) in encoded
            .chunks_exact(pp.linear_code.codeword_len())
            .enumerate()
        {
            let start = i * pp.linear_code.row_len();
            let end = start + pp.linear_code.row_len();
            let row_evals = &poly.evaluations[start..end];
            let expected_encoding = pp.linear_code.encode_wide(row_evals);
            assert_eq!(
                row_chunk,
                expected_encoding.as_slice(),
                "Row {i} encoding mismatch",
            );
        }
    }

    /// Verifies that corrupting the encoded data after commitment results in a different Merkle root.
    #[test]
    fn corrupted_encoding_changes_merkle_root() {
        let (pp, poly) = setup_test_params(3);
        let (mut data, commitment) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly).unwrap();

        if !data.rows.is_empty() {
            data.rows[0] = Int::from(999999);
            let codeword_len = pp.linear_code.codeword_len();
            let corrupted_row = &data.rows[0..codeword_len];
            let new_tree = MerkleTree::new(data.rows_merkle_trees[0].depth, corrupted_row);
            assert_ne!(
                new_tree.root, commitment.roots[0],
                "Corruption should change Merkle root"
            );
        }
    }

    #[test]
    fn batch_commit_with_single_polynomial_is_consistent() {
        let (pp, poly) = setup_test_params(3);

        let batch_result = MultilinearZip::<ZT, _>::batch_commit::<F>(&pp, from_ref(&poly));
        let mut batch_outputs = batch_result.unwrap();
        let (batch_data, batch_commitment) = batch_outputs.remove(0);

        let single_result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        let (single_data, single_commitment) = single_result.unwrap();

        assert_eq!(batch_commitment.roots, single_commitment.roots);
        assert_eq!(batch_data.rows, single_data.rows);
    }

    #[test]
    fn encoded_rows_are_nonzero_for_nonzero_input() {
        let (pp, poly) = setup_test_params(3);
        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly.evaluations,
        );

        assert_eq!(encoded.len(), pp.num_rows * pp.linear_code.codeword_len());
        let non_zero_count = encoded.iter().filter(|&&x| x != Int::from(0)).count();
        assert!(non_zero_count > 0);
    }

    #[test]
    fn commit_produces_correct_merkle_tree_count() {
        let (pp, poly) = setup_test_params(3);
        let (data, _) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly).unwrap();

        assert_eq!(data.rows_merkle_trees.len(), pp.num_rows);
        assert_eq!(data.rows.len(), pp.num_rows * pp.linear_code.codeword_len());
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn encoding_is_consistent_across_threads() {
        use rayon::prelude::*;

        let num_vars = 6;
        let poly_size = 1 << num_vars;
        let evaluations = (1..=poly_size).map(|v| Int::from(v as i32)).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(num_vars, evaluations);

        let results: Vec<Vec<Int<4>>> = (0..10)
            .into_par_iter()
            .map(|_| {
                let mut transcript = MockTranscript::default();
                let code =
                    ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, poly_size, &mut transcript);
                let pp = MultilinearZipParams::new(num_vars, 8, code);

                MultilinearZip::<ZT, _>::encode_rows(
                    &pp,
                    pp.linear_code.codeword_len(),
                    pp.linear_code.row_len(),
                    &poly.evaluations,
                )
            })
            .collect();

        assert!(
            results.windows(2).all(|w| w[0] == w[1]),
            "Parallel encoding runs produced inconsistent results"
        );
    }

    #[test]
    fn commit_succeeds_for_zero_polynomial() {
        let (pp, _) = setup_test_params(3);
        let zero_poly = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(0); 8]);
        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &zero_poly);
        assert!(result.is_ok());

        let (data, commitment) = result.unwrap();
        assert_eq!(commitment.roots.len(), pp.num_rows);
        assert_eq!(data.rows_merkle_trees.len(), pp.num_rows);
    }

    #[test]
    fn commit_succeeds_for_alternating_values() {
        let (pp, _) = setup_test_params(3);
        let alternating = (0..8)
            .map(|i| Int::from(if i % 2 == 0 { 1 } else { -1 }))
            .collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, alternating);
        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_ok());
    }

    #[test]
    fn batch_commit_on_empty_slice_is_ok() {
        let (pp, _) = setup_test_params(3);
        let empty_polys: Vec<DenseMultilinearExtension<Int<INT_LIMBS>>> = vec![];
        let results = MultilinearZip::<ZT, _>::batch_commit::<F>(&pp, &empty_polys);
        assert!(results.is_ok());
        assert!(results.unwrap().is_empty());
    }

    #[test]
    fn encode_rows_succeeds_for_single_row() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 4, &mut transcript);
        let pp = MultilinearZipParams::new(2, 1, code);

        let poly = DenseMultilinearExtension::from_evaluations_vec(2, vec![Int::from(5); 4]);
        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly.evaluations,
        );
        assert_eq!(encoded.len(), pp.linear_code.codeword_len());
    }

    #[test]
    fn merkle_root_integrity_is_maintained() {
        let (pp, _) = setup_test_params(3);
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(42); 8]);
        let (data, commitment) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly).unwrap();

        let codeword_len = pp.linear_code.codeword_len();
        for (i, tree) in data.rows_merkle_trees.iter().enumerate() {
            let start = i * codeword_len;
            let end = start + codeword_len;
            let row_data = &data.rows[start..end];
            let independent_tree = MerkleTree::new(tree.depth, row_data);
            assert_eq!(tree.root, independent_tree.root);
            assert_eq!(commitment.roots[i], independent_tree.root);
        }
    }

    #[test]
    fn matrix_dimensions_are_invariant() {
        let test_cases = vec![(2, 2), (4, 4), (6, 8)];
        for (num_vars, expected_rows) in test_cases {
            assert_eq!(1 << (num_vars / 2), expected_rows);

            let (pp, poly) = setup_test_params(num_vars);
            assert_eq!(pp.num_rows, expected_rows);
            let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
            assert!(result.is_ok());
        }
    }

    #[test]
    #[should_panic]
    fn reject_incompatible_dimensions() {
        let (pp, poly) = setup_test_params(3);
        let incompatible_pp = MultilinearZipParams::new(3, 3, pp.linear_code);
        let _ = MultilinearZip::<ZT, _>::commit::<F>(&incompatible_pp, &poly);
    }

    #[test]
    fn linear_code_preserves_linearity() {
        let (pp, poly) = setup_test_params(4);
        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly.evaluations,
        );
        let row_len = pp.linear_code.row_len();
        let codeword_len = pp.linear_code.codeword_len();
        let row1_evals = &poly.evaluations[0..row_len];
        let row2_evals = &poly.evaluations[row_len..2 * row_len];
        let a = Int::from(3);
        let b = Int::from(5);
        let combined_evals: Vec<_> = (0..row_len)
            .map(|i| a * row1_evals[i] + b * row2_evals[i])
            .collect();
        let combined_encoded = pp.linear_code.encode_wide::<_, Int<4>>(&combined_evals);
        let row1_encoded = &encoded[0..codeword_len];
        let row2_encoded = &encoded[codeword_len..2 * codeword_len];
        let expected_combined: Vec<_> = (0..codeword_len)
            .map(|i| Int::<4>::from(&a) * row1_encoded[i] + Int::<4>::from(&b) * row2_encoded[i])
            .collect();
        assert_eq!(combined_encoded, expected_combined);
    }

    #[test]
    #[should_panic]
    fn commit_panics_if_evaluations_not_multiple_of_row_len() {
        let (pp, mut poly) = setup_test_params(4);
        poly.evaluations.truncate(15);
        assert_eq!(poly.evaluations.len(), 15);
        let _ = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
    }

    #[test]
    fn commit_with_many_variables() {
        let num_vars = 16;
        let (pp, poly) = setup_test_params(num_vars);
        assert_eq!(pp.num_vars, num_vars);
        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_ok());
        let (_, commitment) = result.unwrap();
        assert_eq!(commitment.roots.len(), pp.num_rows);
    }

    #[test]
    fn commit_with_smallest_matrix_arrangement() {
        let (pp, poly) = setup_test_params(2);
        assert_eq!(pp.num_rows, 2);
        assert_eq!(pp.linear_code.row_len(), 2);
        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_ok());
        let (data, commitment) = result.unwrap();
        assert_eq!(commitment.roots.len(), 2);
        assert_eq!(data.rows_merkle_trees.len(), 2);
    }

    #[test]
    fn encode_rows_handles_large_integer_values() {
        let (pp, _) = setup_test_params(3);
        let max_val = Int::<INT_LIMBS>::from(i64::MAX);
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, vec![max_val; 8]);
        let encoded_rows = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly.evaluations,
        );
        assert_eq!(
            encoded_rows.len(),
            pp.num_rows * pp.linear_code.codeword_len()
        );
    }

    #[test]
    #[should_panic(expected = "leaves.len().is_power_of_two()")]
    fn merkle_tree_new_panics_on_non_power_of_two_leaves() {
        let leaves_data: Vec<Int<INT_LIMBS>> = (0..7).map(Int::from).collect();
        let merkle_depth = 3;
        let _ = MerkleTree::new(merkle_depth, &leaves_data);
    }

    #[test]
    fn verifier_rejects_commitment_with_bad_proximity() {
        fn evaluate_in_field<'cfg>(
            evaluations: &[Int<INT_LIMBS>],
            point: &[RandomField<'cfg, FIELD_LIMBS>],
            config: ConfigRef<'cfg, FIELD_LIMBS>,
        ) -> F<'cfg> {
            let num_vars = point.len();
            assert_eq!(evaluations.len(), 1 << num_vars);
            let mut current_evals: Vec<F> = evaluations.map_to_field(config);
            for p in point.iter().take(num_vars) {
                let one_minus_p_i = FieldMap::<F>::map_to_field(&1i32, config) - p;
                let mut next_evals = Vec::with_capacity(current_evals.len() / 2);
                for j in (0..current_evals.len()).step_by(2) {
                    let val = current_evals[j] * one_minus_p_i + current_evals[j + 1] * p;
                    next_evals.push(val);
                }
                current_evals = next_evals;
            }
            current_evals[0]
        }

        type F<'cfg> = RandomField<'cfg, FIELD_LIMBS>;
        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();
        let n = 3;
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);
        let evaluations: Vec<_> = (0..poly_size)
            .map(|_| Int::<INT_LIMBS>::from(i8::rand(&mut rng)))
            .collect();
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);
        let point_int: Vec<_> = (0..n).map(|_| Int::<INT_LIMBS>::random(&mut rng)).collect();
        let point_f = point_int.map_to_field(config);

        let (mut data, comm) = TestZip::commit::<F>(&param, &mle).unwrap();
        if !data.rows.is_empty() {
            data.rows[0] += Int::<{ 4 * INT_LIMBS }>::from(1);
        }

        let mut prover_transcript = PcsTranscript::new();
        TestZip::open(
            &param,
            &mle,
            &data,
            &point_f,
            config,
            &mut prover_transcript,
        )
        .unwrap();
        let proof = prover_transcript.into_proof();

        let mut verifier_transcript = PcsTranscript::from_proof(&proof);
        let eval = evaluate_in_field(&mle.evaluations, &point_f, config);
        let verification_result = TestZip::verify(
            &param,
            &comm,
            &point_f,
            eval,
            &mut verifier_transcript,
            config,
        );

        assert!(verification_result.is_err());
    }

    #[test]
    fn proof_size_is_correct_for_parameters() {
        fn calculate_expected_proof_size_bits<ZT: ZipTypes, LC: LinearCode<ZT>>(
            pp: &MultilinearZipParams<ZT, LC>,
        ) -> usize {
            let size_of_zt_k = size_of::<ZT::K>();
            let size_of_zt_m = size_of::<ZT::M>();
            let size_of_f_b = size_of::<BigInt<FIELD_LIMBS>>();
            let size_of_hash = Keccak256::output_size();
            let size_of_path_len = size_of::<u64>();

            let codeword_len = pp.linear_code.codeword_len();
            let merkle_depth = codeword_len.next_power_of_two().ilog2() as usize;

            let proximity_phase_size =
                pp.linear_code.num_proximity_testing() * pp.linear_code.row_len() * size_of_zt_m;

            let column_values_size = pp.num_rows * size_of_zt_k;
            let single_merkle_proof_size = size_of_path_len + merkle_depth * size_of_hash;
            let all_merkle_proofs_size = pp.num_rows * single_merkle_proof_size;
            let size_per_column_opening = column_values_size + all_merkle_proofs_size;
            let column_opening_phase_size =
                pp.linear_code.num_column_opening() * size_per_column_opening;

            let evaluation_phase_size = pp.linear_code.row_len() * size_of_f_b;

            (proximity_phase_size + column_opening_phase_size + evaluation_phase_size) * 8
        }

        type F<'cfg> = RandomField<'cfg, FIELD_LIMBS>;
        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();
        let num_vars = 4;
        let poly_size = 1 << num_vars;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);
        let evaluations: Vec<_> = (0..poly_size)
            .map(|_| Int::<INT_LIMBS>::from(i8::rand(&mut rng)))
            .collect();
        let mle = DenseMultilinearExtension::from_evaluations_slice(num_vars, &evaluations);
        let point_int: Vec<_> = (0..num_vars)
            .map(|_| Int::<INT_LIMBS>::random(&mut rng))
            .collect();
        let point_f: Vec<F> = point_int.map_to_field(config);

        let (data, _) = TestZip::commit::<F>(&param, &mle).unwrap();
        let mut prover_transcript = PcsTranscript::new();
        TestZip::open(
            &param,
            &mle,
            &data,
            &point_f,
            config,
            &mut prover_transcript,
        )
        .unwrap();
        let proof = prover_transcript.into_proof();

        let actual_proof_size_bits = proof.len() * 8;
        let expected_proof_size_bits = calculate_expected_proof_size_bits(&param);

        assert_eq!(actual_proof_size_bits, expected_proof_size_bits);
    }
}
