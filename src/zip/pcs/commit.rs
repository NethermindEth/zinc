use ark_std::{vec, vec::Vec};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment, MultilinearZipData},
    utils::{validate_input, MerkleTree},
};
use crate::{
    poly_z::mle::DenseMultilinearExtension,
    traits::{Field, ZipTypes},
    zip::{
        code::LinearCode,
        pcs::structs::MultilinearZipParams,
        utils::{div_ceil, num_threads, parallelize_iter},
        Error,
    },
};

impl<ZT: ZipTypes, LC: LinearCode<ZT>> MultilinearZip<ZT, LC> {
    /// Creates a commitment to a multilinear polynomial using the ZIP PCS scheme.
    ///
    /// This function implements the commitment phase of the ZIP polynomial commitment scheme,
    /// which encodes the polynomial evaluations using a linear error-correcting code and
    /// creates Merkle tree commitments to the encoded rows.
    ///
    /// # Algorithm
    /// 1. Validates that the polynomial dimensions match the parameters
    /// 2. Arranges polynomial evaluations into a matrix with `pp.num_rows` rows
    /// 3. Encodes each row using the linear code (expanding row_len → codeword_len)
    /// 4. Creates a Merkle tree for each encoded row
    /// 5. Returns both the full data (for proving) and the commitment (Merkle roots)
    ///
    /// # Type Parameters
    /// - `F`: Field type for validation (not directly used in commitment)
    ///
    /// # Parameters
    /// - `pp`: Public parameters containing:
    ///   - `num_vars`: Number of variables the scheme supports (log of evaluation domain size)
    ///   - `num_rows`: Number of rows in matrix arrangement
    ///   - `linear_code`: Error-correcting code for proximity testing
    /// - `poly`: Multilinear polynomial to commit to, with evaluations over boolean hypercube
    ///
    /// # Returns
    /// - `MultilinearZipData`: Full commitment data including encoded rows and Merkle trees (kept by prover)
    /// - `MultilinearZipCommitment`: Commitment containing only Merkle roots (sent to verifier)
    ///
    /// # Errors
    /// - `InvalidPcsParam`: If polynomial has more variables than `pp.num_vars` supports
    ///
    /// # Complexity
    /// - Time: O(n × expansion_factor) where n = 2^{num_vars} (dominated by encoding)
    /// - Space: O(n × expansion_factor) for encoded rows
    ///
    /// # Panics
    /// Panics if the number of Merkle trees doesn't match `pp.num_rows` (internal logic error)
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

        let rows = Self::encode_rows(pp, codeword_len, row_len, poly);

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

    /// Execute a commit phase without constructing Merkle treesю
    /// Needed for tests only.
    #[allow(dead_code)]
    pub fn commit_no_merkle<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        poly: &DenseMultilinearExtension<ZT::N>,
    ) -> Result<(MultilinearZipData<ZT::K>, MultilinearZipCommitment), Error> {
        validate_input("commit", pp.num_vars, [poly], None::<&[F]>)?;

        let row_len = pp.linear_code.row_len();
        let codeword_len = pp.linear_code.codeword_len();

        let rows = Self::encode_rows(pp, codeword_len, row_len, poly);

        Ok((
            MultilinearZipData::new(rows, vec![]),
            MultilinearZipCommitment { roots: vec![] },
        ))
    }

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

    /// Encodes the rows of the polynomial concatenating each encoded row
    pub fn encode_rows(
        pp: &MultilinearZipParams<ZT, LC>,
        codeword_len: usize,
        row_len: usize,
        poly: &DenseMultilinearExtension<ZT::N>,
    ) -> Vec<ZT::K> {
        let rows_per_thread = div_ceil(pp.num_rows, num_threads());
        let mut encoded_rows = vec![ZT::K::default(); pp.num_rows * codeword_len];

        parallelize_iter(
            encoded_rows
                .chunks_exact_mut(rows_per_thread * codeword_len)
                .zip(poly.evaluations.chunks_exact(rows_per_thread * row_len)),
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
    use crate::field::{BigInt, ConfigRef};
    use crate::traits::{FieldMap, ZipTypes};
    use crate::transcript::KeccakTranscript;
    use crate::zip::code_raa::RaaCode;
    use crate::zip::pcs::structs::{MultilinearZipCommitment, MultilinearZipData};
    use crate::zip::pcs_transcript::PcsTranscript;
    use crate::zip::Error;
    use crate::{
        define_random_field_zip_types,
        field::{Int, RandomField},
        field_config, implement_random_field_zip_types,
        poly_z::mle::DenseMultilinearExtension,
        traits::Integer,
        zip::{
            code::{DefaultLinearCodeSpec, LinearCode, ZipLinearCode},
            pcs::{
                structs::{MultilinearZip, MultilinearZipParams, ZipTranscript},
                MerkleTree,
            },
        },
    };
    use ark_std::println;
    use ark_std::vec::Vec;
    use ark_std::{collections::BTreeSet, ops::Range, vec, UniformRand};
    use crypto_bigint::Random;
    use sha3::{Digest, Keccak256};

    const INT_LIMBS: usize = 1;
    const FIELD_LIMBS: usize = 4;

    define_random_field_zip_types!();
    implement_random_field_zip_types!(INT_LIMBS);

    type ZT = RandomFieldZipTypes<INT_LIMBS>;

    type F<'cfg> = RandomField<'cfg, FIELD_LIMBS>;

    type TestZip<LC> = MultilinearZip<ZT, LC>;

    type LC = RaaCode<ZT>;

    #[derive(Default)]
    pub struct MockTranscript {
        pub counter: i64,
    }

    impl<L: Integer> ZipTranscript<L> for MockTranscript {
        fn get_encoding_element(&mut self) -> L {
            self.counter += 1;
            L::from(self.counter)
        }

        fn get_u64(&mut self) -> u64 {
            self.counter += 1;
            self.counter as u64
        }

        fn sample_unique_columns(
            &mut self,
            range: Range<usize>,
            columns: &mut BTreeSet<usize>,
            count: usize,
        ) -> usize {
            self.counter += 1;

            let mut inserted = 0;
            for i in range.clone() {
                if columns.insert(i) {
                    inserted += 1;
                    if inserted == count {
                        break;
                    }
                }
            }

            inserted
        }
    }

    #[allow(clippy::type_complexity)]
    fn setup_zip_params(
        num_vars: usize,

        num_rows: usize,

        poly_size: usize,

        poly_input: Range<i32>,
        poly_num_vars: usize,
    ) -> (
        Result<
            (
                MultilinearZipData<<ZT as ZipTypes>::K>,
                MultilinearZipCommitment,
            ),
            Error,
        >,
        MultilinearZipParams<ZT, LC>,
    ) {
        let mut transcript = MockTranscript::default();

        let code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut transcript);

        let pp = MultilinearZipParams::new(num_vars, num_rows, code);

        let evaluations = poly_input.map(Int::from).collect();

        let poly = DenseMultilinearExtension::from_evaluations_vec(poly_num_vars, evaluations);

        (MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly), pp)
    }

    #[test]
    fn commit_rejects_too_many_variables() {
        let (result, _) = setup_zip_params(3, 4, 8, 1..17, 4);
        assert!(result.is_err());
    }

    #[test]
    fn commit_deterministic() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);

        let pp = MultilinearZipParams::new(3, 4, code);

        let evaluations = (1..=8).map(Int::from).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, evaluations);

        let result1 = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly).unwrap();
        let result2 = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly).unwrap();

        assert_eq!(result1.1.roots, result2.1.roots);
    }

    #[test]
    fn commit_different_polynomials_different_commitments() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);

        let pp = MultilinearZipParams::new(3, 4, code);

        let poly1 = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(1); 8]);
        let poly2 = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(2); 8]);

        let (_, commitment1) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly1).unwrap();
        let (_, commitment2) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly2).unwrap();

        assert_ne!(commitment1.roots, commitment2.roots);
    }

    #[test]
    fn commit_small_polynomial() {
        let (result, _) = setup_zip_params(4, 4, 16, 1..17, 4);
        assert!(result.is_ok());
    }

    #[test]
    fn commit_edge_case_two_variables() {
        let (result, _) = setup_zip_params(2, 2, 4, 1..5, 2);
        assert!(result.is_ok());
    }

    #[test]
    fn merkle_tree_depth_correct() {
        let (result, pp) = setup_zip_params(3, 4, 8, 1..9, 3);
        let expected_depth = pp.linear_code.codeword_len().next_power_of_two().ilog2() as usize;

        for tree in &result.unwrap().0.rows_merkle_trees {
            assert_eq!(tree.depth, expected_depth);
        }
    }

    #[test]
    fn commit_no_merkle_behaves_correctly() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);

        let pp = MultilinearZipParams::new(3, 4, code);

        let evaluations = (1..=8).map(Int::from).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, evaluations);

        let result = MultilinearZip::<ZT, _>::commit_no_merkle::<F>(&pp, &poly);
        assert!(result.is_ok());

        let (data, commitment) = result.unwrap();
        assert_eq!(data.rows.len(), pp.num_rows * pp.linear_code.codeword_len());
        assert!(data.rows_merkle_trees.is_empty());
        assert!(commitment.roots.is_empty());
    }

    #[test]
    fn batch_commit_behaves_correctly() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);

        let pp = MultilinearZipParams::new(3, 4, code);

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
    fn encode_rows_behaves_correctly() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);

        let pp = MultilinearZipParams::new(3, 4, code);

        let evaluations = (1..=8).map(Int::from).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, evaluations);

        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly,
        );

        assert_eq!(encoded.len(), pp.num_rows * pp.linear_code.codeword_len());
    }

    #[test]
    fn test_encode_rows_output_size() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);
        let pp = MultilinearZipParams::new(3, 4, code);

        let poly = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(1); 8]);
        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly,
        );

        assert_eq!(encoded.len(), pp.num_rows * pp.linear_code.codeword_len());

        let non_zero_count = encoded.iter().filter(|&&x| x != Int::from(0)).count();
        assert!(non_zero_count > 0);
    }

    #[test]
    fn test_commit_merkle_tree_count() {
        let (result, pp) = setup_zip_params(3, 4, 8, 1..9, 3);
        let (data, _) = result.unwrap();
        assert_eq!(data.rows_merkle_trees.len(), pp.num_rows);
        assert_eq!(data.rows.len(), pp.num_rows * pp.linear_code.codeword_len());
    }

    #[test]
    fn test_encode_rows_parallelization_consistency() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 16, &mut transcript);
        let pp = MultilinearZipParams::new(4, 4, code);

        let poly1 = DenseMultilinearExtension::from_evaluations_vec(4, vec![Int::from(7); 16]);
        let poly2 = DenseMultilinearExtension::from_evaluations_vec(4, vec![Int::from(7); 16]);

        let encoded1 = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly1,
        );
        let encoded2 = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly2,
        );

        assert_eq!(encoded1, encoded2);
    }

    #[test]
    fn test_commit_with_zero_polynomial() {
        let (result, pp) = setup_zip_params(3, 4, 8, 1..9, 3);

        assert!(result.is_ok());
        let (data, commitment) = result.unwrap();

        assert_eq!(commitment.roots.len(), pp.num_rows);
        assert_eq!(data.rows_merkle_trees.len(), pp.num_rows);
    }

    #[test]
    fn test_commit_with_alternating_values() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);
        let pp = MultilinearZipParams::new(3, 4, code);

        let alternating = (0..8)
            .map(|i| Int::from(if i % 2 == 0 { 1 } else { -1 }))
            .collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, alternating);

        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_ok());
    }

    #[test]
    fn test_batch_commit_empty_slice() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);
        let pp = MultilinearZipParams::new(3, 4, code);

        let empty_polys: Vec<DenseMultilinearExtension<Int<INT_LIMBS>>> = vec![];
        let results = MultilinearZip::<ZT, _>::batch_commit::<F>(&pp, &empty_polys);

        assert!(results.is_ok());
        assert_eq!(results.unwrap().len(), 0);
    }

    #[test]
    fn test_encode_rows_with_single_row() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 4, &mut transcript);
        let pp = MultilinearZipParams::new(2, 1, code); // Single row edge case

        let poly = DenseMultilinearExtension::from_evaluations_vec(2, vec![Int::from(5); 4]);
        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly,
        );

        assert_eq!(encoded.len(), pp.linear_code.codeword_len());
    }

    #[test]
    fn test_semantic_correctness_encoded_rows_match_linear_code() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);
        let pp = MultilinearZipParams::new(3, 4, code);

        let evaluations = (1..=8).map(Int::from).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, evaluations);
        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly,
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
                "Row {i} encoding mismatch"
            );
        }
    }

    #[test]
    fn test_merkle_root_integrity() {
        let (result, pp) = setup_zip_params(3, 4, 8, 1..9, 3);
        let (data, commitment) = result.unwrap();
        let codeword_len = pp.linear_code.codeword_len();
        for (i, tree) in data.rows_merkle_trees.iter().enumerate() {
            let start = i * codeword_len;
            let end = start + codeword_len;
            let row_data = &data.rows[start..end];
            let independent_tree = MerkleTree::new(tree.depth, row_data);

            assert_eq!(
                tree.root, independent_tree.root,
                "Merkle root mismatch for row {i}"
            );
            assert_eq!(
                commitment.roots[i], independent_tree.root,
                "Commitment root mismatch for row {i}"
            );
        }
    }

    #[test]
    fn test_matrix_dimension_invariant() {
        let test_cases = vec![
            (2, 2), // 2^2 = 4 evaluations, sqrt(4) = 2 rows
            (4, 4), // 2^4 = 16 evaluations, sqrt(16) = 4 rows
            (6, 8), // 2^6 = 64 evaluations, sqrt(64) = 8 rows
        ];

        for (num_vars, expected_rows) in test_cases {
            assert_eq!(
                expected_rows,
                1 << (num_vars / 2),
                "Expected rows should equal sqrt(2^num_vars) for {num_vars} vars"
            );
            let mut transcript = MockTranscript::default();
            let total_evals = 1 << num_vars;
            let code =
                ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, total_evals, &mut transcript);
            let pp = MultilinearZipParams::new(num_vars, expected_rows, code);
            let evaluations = vec![Int::from(1); total_evals];
            let poly = DenseMultilinearExtension::from_evaluations_vec(num_vars, evaluations);
            let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
            assert!(
                result.is_ok(),
                "Commitment should work for {num_vars} vars with {expected_rows} rows"
            );
        }
    }

    #[test]
    fn test_reject_incompatible_dimensions() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 4, &mut transcript);
        let pp = MultilinearZipParams::new(3, 3, code); // 2^3 = 8 evaluations
        let evaluations = vec![Int::from(1); 8];
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, evaluations);
        let result = ark_std::panic::catch_unwind(|| {
            MultilinearZip::<ZT, _>::encode_rows(
                &pp,
                pp.linear_code.codeword_len(),
                pp.linear_code.row_len(),
                &poly,
            )
        });

        assert!(
            result.is_err() || pp.num_rows * pp.linear_code.row_len() != poly.evaluations.len()
        );
    }

    #[test]
    fn test_corrupted_encoding_detection() {
        let (result, pp) = setup_zip_params(3, 4, 8, 1..9, 3);
        let (mut data, commitment) = result.unwrap();
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
    fn test_encoding_consistency_across_threads() {
        #[cfg(feature = "parallel")]
        use rayon::prelude::*;

        use ark_std::cfg_into_iter;

        let results = cfg_into_iter!(0..16)
            .map(|_| {
                let mut transcript = MockTranscript::default();
                let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 16, &mut transcript);
                let pp = MultilinearZipParams::new(6, 8, code);
                let evaluations = (1..=64).map(Int::from).collect();
                let poly = DenseMultilinearExtension::from_evaluations_vec(6, evaluations);
                MultilinearZip::<ZT, _>::encode_rows(
                    &pp,
                    pp.linear_code.codeword_len(),
                    pp.linear_code.row_len(),
                    &poly,
                )
            })
            .collect::<Vec<_>>();

        let first_encoding = results[0].clone();

        assert!(
            results.iter().all(|encoding| encoding == &first_encoding),
            "All encodings should be consistent across threads"
        );
    }

    #[test]
    fn test_commitment_binding_property() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);
        let pp = MultilinearZipParams::new(3, 4, code);
        let poly1 = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(1); 8]);
        let poly2 = DenseMultilinearExtension::from_evaluations_vec(3, vec![Int::from(2); 8]);

        let (_, commitment1) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly1).unwrap();
        let (_, commitment2) = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly2).unwrap();
        assert_ne!(
            commitment1.roots, commitment2.roots,
            "Different polynomials must have different commitments"
        );
    }

    #[test]
    fn test_row_wise_linear_combination_property() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 4, &mut transcript);
        let pp = MultilinearZipParams::new(4, 4, code);

        let evaluations = (1..=16).map(Int::from).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(4, evaluations);

        let encoded = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly,
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

        assert_eq!(
            combined_encoded, expected_combined,
            "Linear code must preserve linearity"
        );
    }

    #[test]
    #[should_panic]
    fn commit_panics_if_poly_evals_not_multiple_of_row_len() {
        let mut transcript = MockTranscript::default();
        let poly_size = 16; // num_vars = 4
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, poly_size, &mut transcript);
        let pp = MultilinearZipParams::new(4, 4, code);
        let evaluations = vec![Int::from(1); 16];
        let mut poly = DenseMultilinearExtension::from_evaluations_vec(4, evaluations);
        poly.evaluations.truncate(15);
        assert_eq!(poly.evaluations.len(), 15);
        let _ = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
    }

    #[test]
    fn commit_with_all_evaluations_identical_and_nonzero() {
        let mut transcript = MockTranscript::default();
        let poly_size = 8; // num_vars = 3
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, poly_size, &mut transcript);
        let pp = MultilinearZipParams::new(3, 4, code);
        let evaluations1 = vec![Int::from(42); 8];
        let poly1 = DenseMultilinearExtension::from_evaluations_vec(3, evaluations1);
        let commit_result1 = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly1);
        assert!(commit_result1.is_ok());
        let (_, commitment1) = commit_result1.unwrap();
        let evaluations2 = vec![Int::from(43); 8];
        let poly2 = DenseMultilinearExtension::from_evaluations_vec(3, evaluations2);
        let commit_result2 = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly2);
        assert!(commit_result2.is_ok());
        let (_, commitment2) = commit_result2.unwrap();
        assert_ne!(
            commitment1.roots, commitment2.roots,
            "Commitments to different constant polynomials should not be the same"
        );
    }

    #[test]
    fn commit_with_many_variables() {
        let mut transcript = MockTranscript::default();
        let num_vars = 16;
        let poly_size = 1 << num_vars; // 65,536 evaluations

        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, poly_size, &mut transcript);
        let pp = MultilinearZip::<ZT, _>::setup(poly_size, code);
        assert_eq!(pp.num_vars, num_vars);
        let evaluations: Vec<_> = (0..poly_size).map(|v| Int::from(v as i32)).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(num_vars, evaluations);
        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_ok());
        let (_, commitment) = result.unwrap();
        assert_eq!(commitment.roots.len(), pp.num_rows);
    }

    #[test]
    fn commit_with_smallest_matrix_arrangement() {
        let mut transcript = MockTranscript::default();
        let num_vars = 2;
        let poly_size = 1 << num_vars; // 4 evaluations

        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, poly_size, &mut transcript);
        let pp = MultilinearZip::<ZT, _>::setup(poly_size, code);
        assert_eq!(pp.num_rows, 2, "Expected a 2x2 matrix arrangement");
        assert_eq!(pp.linear_code.row_len(), 2);
        let evaluations = vec![Int::from(1), Int::from(2), Int::from(3), Int::from(4)];
        let poly = DenseMultilinearExtension::from_evaluations_vec(num_vars, evaluations);
        let result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(result.is_ok());
        let (data, commitment) = result.unwrap();
        assert_eq!(
            commitment.roots.len(),
            2,
            "Commitment for a 2x2 matrix should have two Merkle roots"
        );
        assert_eq!(
            data.rows_merkle_trees.len(),
            2,
            "Data for a 2x2 matrix should have two Merkle trees"
        );
    }

    #[test]
    fn batch_commit_with_single_polynomial() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);
        let pp = MultilinearZipParams::new(3, 4, code);
        let poly =
            DenseMultilinearExtension::from_evaluations_vec(3, (1..=8).map(Int::from).collect());
        let batch_result = MultilinearZip::<ZT, _>::batch_commit::<F>(&pp, &[poly.clone()]);
        assert!(batch_result.is_ok());
        let mut batch_outputs = batch_result.unwrap();
        assert_eq!(batch_outputs.len(), 1);
        let (batch_data, batch_commitment) = batch_outputs.remove(0);
        let single_result = MultilinearZip::<ZT, _>::commit::<F>(&pp, &poly);
        assert!(single_result.is_ok());
        let (single_data, single_commitment) = single_result.unwrap();
        assert_eq!(
            batch_commitment.roots, single_commitment.roots,
            "Commitment roots should be identical"
        );
        assert_eq!(
            batch_data.rows, single_data.rows,
            "Encoded rows data should be identical"
        );
    }

    #[test]
    fn encode_rows_with_large_integer_values() {
        let mut transcript = MockTranscript::default();
        let code = ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, 8, &mut transcript);
        let pp = MultilinearZipParams::new(3, 4, code);
        let max_val = Int::<INT_LIMBS>::from(i64::MAX);
        let evaluations = vec![max_val; 8];
        let poly = DenseMultilinearExtension::from_evaluations_vec(3, evaluations);
        let encoded_rows = MultilinearZip::<ZT, _>::encode_rows(
            &pp,
            pp.linear_code.codeword_len(),
            pp.linear_code.row_len(),
            &poly,
        );
        let expected_len = pp.num_rows * pp.linear_code.codeword_len();
        assert_eq!(
            encoded_rows.len(),
            expected_len,
            "Encoded rows vector has an incorrect length"
        );
    }

    #[test]
    #[should_panic(expected = "leaves.len().is_power_of_two()")]
    fn merkle_tree_new_panics_on_non_power_of_two_leaves() {
        let leaves_data: Vec<Int<INT_LIMBS>> = (0..7).map(Int::from).collect();
        let merkle_depth = 3; // Depth for 8 leaves, but we only provide 7.
        let _ = MerkleTree::new(merkle_depth, &leaves_data);
    }

    #[test]
    fn test_verifier_rejects_commitment_with_bad_proximity() {
        fn evaluate_in_field<'cfg>(
            evaluations: &[Int<INT_LIMBS>],
            point: &[RandomField<'cfg, FIELD_LIMBS>],
            config: ConfigRef<'cfg, FIELD_LIMBS>,
        ) -> F<'cfg> {
            let num_vars = point.len();
            assert_eq!(evaluations.len(), 1 << num_vars);
            let mut current_evals: Vec<F> = evaluations.map_to_field(config);
            for p in point.iter().take(num_vars) {
                let one_minus_p_i = FieldMap::<F>::map_to_field(&1i32, config) - p; // F::one is correctly mapped
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
        let param = TestZip::<LC>::setup(poly_size, linear_code);
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
        assert!(
            verification_result.is_err(),
            "Verifier should reject a proof based on a corrupted codeword"
        );
    }

    #[test]
    fn test_proof_size_is_correct_for_parameters() {
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
            let evaluation_phase_size = pp.linear_code.row_len() * size_of_f_b; // Use the corrected size here
            (proximity_phase_size + column_opening_phase_size + evaluation_phase_size) * 8
        }
        type F<'cfg> = RandomField<'cfg, FIELD_LIMBS>;
        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        let num_vars = 4;
        let poly_size = 1 << num_vars;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code =
            ZipLinearCode::<ZT>::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = MultilinearZip::<ZT, _>::setup(poly_size, linear_code);

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

        println!(
            "Actual proof size:   {} bits ({} bytes)",
            actual_proof_size_bits,
            proof.len()
        );
        println!(
            "Expected proof size: {} bits ({} bytes)",
            expected_proof_size_bits,
            expected_proof_size_bits / 8
        );
        assert_eq!(
            actual_proof_size_bits, expected_proof_size_bits,
            "Actual proof size does not match the expected size calculated from parameters"
        );
    }
}
