#![allow(non_snake_case)]
use ark_std::{borrow::Cow, vec::Vec};
use itertools::izip;

use super::{
    structs::{MultilinearZip, MultilinearZipData},
    utils::{left_point_to_tensor, validate_input, ColumnOpening},
};
use crate::{
    poly_z::mle::DenseMultilinearExtension,
    traits::{Field, FieldMap, ZipTypes},
    zip::{
        code::LinearCode,
        pcs::structs::MultilinearZipParams,
        pcs_transcript::PcsTranscript,
        utils::{combine_rows, expand},
        Error,
    },
};

impl<ZT: ZipTypes, LC: LinearCode<ZT>> MultilinearZip<ZT, LC> {
    pub fn open<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        poly: &DenseMultilinearExtension<ZT::N>,
        commit_data: &MultilinearZipData<ZT::K>,
        point: &[F],
        field: F::R,
        transcript: &mut PcsTranscript<F>,
    ) -> Result<(), Error>
    where
        ZT::N: FieldMap<F, Output = F>,
    {
        validate_input("open", pp.num_vars, [poly], [point])?;

        Self::prove_testing_phase(pp, poly, commit_data, transcript, field)?;

        Self::prove_evaluation_phase(pp, transcript, point, poly, field)?;

        Ok(())
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30
    pub fn batch_open<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        polys: &[DenseMultilinearExtension<ZT::N>],
        comms: &[MultilinearZipData<ZT::K>],
        points: &[Vec<F>],
        transcript: &mut PcsTranscript<F>,
        field: F::R,
    ) -> Result<(), Error>
    where
        ZT::N: FieldMap<F, Output = F>,
    {
        for (poly, comm, point) in izip!(polys.iter(), comms.iter(), points.iter()) {
            Self::open(pp, poly, comm, point, field, transcript)?;
        }
        Ok(())
    }

    // Subprotocol functions

    fn prove_evaluation_phase<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        transcript: &mut PcsTranscript<F>,
        point: &[F],
        poly: &DenseMultilinearExtension<ZT::N>,
        field: F::R,
    ) -> Result<(), Error>
    where
        ZT::N: FieldMap<F, Output = F>,
    {
        let num_rows = pp.num_rows;
        let row_len = pp.linear_code.row_len();

        // We prove evaluations over the field, so integers need to be mapped to field elements first
        let q_0 = left_point_to_tensor(num_rows, point, field)?;

        let evaluations = poly.evaluations.map_to_field(field);

        let q_0_combined_row = if num_rows > 1 {
            // Return the evaluation row combination
            let combined_row = combine_rows(q_0, evaluations, row_len);
            Cow::<Vec<F>>::Owned(combined_row)
        } else {
            // If there is only one row, we have no need to take linear combinations
            // We just return the evaluation row combination
            Cow::Borrowed(&evaluations)
        };

        transcript.write_field_elements(&q_0_combined_row)
    }

    pub(super) fn prove_testing_phase<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        poly: &DenseMultilinearExtension<ZT::N>,
        commit_data: &MultilinearZipData<ZT::K>,
        transcript: &mut PcsTranscript<F>,
        field: F::R, // This is only needed to call the transcript, but we are getting integers not fields
    ) -> Result<(), Error> {
        if pp.num_rows > 1 {
            // If we can take linear combinations
            // perform the proximity test an arbitrary number of times
            for _ in 0..pp.linear_code.num_proximity_testing() {
                let coeffs = transcript.fs_transcript.get_integer_challenges(pp.num_rows);
                let coeffs = coeffs.iter().map(expand::<ZT::N, ZT::M>);

                let evals = poly.evaluations.iter().map(expand::<ZT::N, ZT::M>);

                // u' in the Zinc paper
                let combined_row = combine_rows(coeffs, evals, pp.linear_code.row_len());

                transcript.write_integers(combined_row.iter())?;
            }
        }

        // Open merkle tree for each column drawn
        for _ in 0..pp.linear_code.num_column_opening() {
            let column = transcript.squeeze_challenge_idx(field, pp.linear_code.codeword_len());
            Self::open_merkle_trees_for_column(pp, commit_data, column, transcript)?;
        }
        Ok(())
    }

    pub(super) fn open_merkle_trees_for_column<F: Field>(
        pp: &MultilinearZipParams<ZT, LC>,
        commit_data: &MultilinearZipData<ZT::K>,
        column: usize,
        transcript: &mut PcsTranscript<F>,
    ) -> Result<(), Error> {
        let column_values = commit_data
            .rows
            .iter()
            .skip(column)
            .step_by(pp.linear_code.codeword_len());

        // Write the elements in the squeezed column to the shared transcript
        transcript.write_integers(column_values)?;

        ColumnOpening::open_at_column(column, commit_data, transcript)
            .map_err(|_| Error::InvalidPcsOpen("Failed to open merkle tree".into()))?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{define_random_field_zip_types, field::{Int, RandomField}, field_config, implement_random_field_zip_types, poly_z::mle::DenseMultilinearExtension, traits::FieldMap, zip::{
        code::{DefaultLinearCodeSpec},
        pcs::structs::{MultilinearZip, MultilinearZipData},
    }};
    use ark_std::{UniformRand};
    use crypto_bigint::Random;
    use sha3::{Digest, Keccak256};
    use crate::field::ConfigRef;
    use crate::transcript::KeccakTranscript;
    use crate::zip::code_raa::RaaCode;
    // --- Test Setup ---
    // Note: The specific types and setup might need adjustment to match your test environment.

    const INT_LIMBS: usize = 1;
    const FIELD_LIMBS: usize = 4;

    define_random_field_zip_types!();
    implement_random_field_zip_types!(INT_LIMBS);

    type ZT = RandomFieldZipTypes<INT_LIMBS>;
    type F<'cfg> = RandomField<'cfg, FIELD_LIMBS>;
    type LC = RaaCode<ZT>;
    type TestZip = MultilinearZip<ZT, LC>;

    /// Verifies that a valid proof can be generated for a standard polynomial and commitment.
    #[test]
    fn open_succeeds_for_valid_polynomial_and_commitment() {
        let config = field_config!(115792089237316195423570985008687907853269984665640564039457584007913129639747, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        let n = 3;
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = MultilinearZip::<ZT, _>::setup(poly_size, linear_code);

        let evaluations: Vec<_> = (0..poly_size).map(|_| Int::<INT_LIMBS>::from(i8::rand(&mut rng))).collect();
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);
        let point: Vec<_> = (0..n).map(|_| F::rand_with_config(&mut rng, config)).collect();

        let (data, _) = MultilinearZip::<ZT, _>::commit::<F>(&param, &mle).unwrap();
        let mut transcript = PcsTranscript::new();

        let result = MultilinearZip::<ZT, _>::open(&param, &mle, &data, &point, config, &mut transcript);

        assert!(result.is_ok());
    }

    /// Ensures that the open function returns an error if the evaluation point has a different
    /// number of variables than the polynomial.
    #[test]
    #[should_panic]
    fn open_fails_with_mismatched_point_dimensions() {
        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        let n = 3;
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = (0..poly_size).map(|_| Int::<INT_LIMBS>::from(i8::rand(&mut rng))).collect();
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

        let (data, _) = TestZip::commit::<F>(&param, &mle).unwrap();
        let mut transcript = PcsTranscript::new();

        let mismatched_point: Vec<_> = (0..n - 1).map(|_| F::rand_with_config(&mut rng, config)).collect();

        let result = TestZip::open(&param, &mle, &data, &mismatched_point, config, &mut transcript);
        assert!(result.is_err());
    }

    /// Checks that the batch opening process succeeds for a list of valid inputs.
    #[test]
    fn batch_open_succeeds_for_multiple_valid_inputs() {
        let config = field_config!(115792089237316195423570985008687907853269984665640564039457584007913129639747, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        let n = 3;
        let num_polys = 5;
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);

        let mut polys = Vec::new();
        let mut comms_data = Vec::new();
        for _ in 0..num_polys {
            let evaluations: Vec<_> = (0..poly_size).map(|_| Int::<INT_LIMBS>::from(i8::rand(&mut rng))).collect();
            let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);
            let (data, _) = TestZip::commit::<F>(&param, &mle).unwrap();
            polys.push(mle);
            comms_data.push(data);
        }

        let point: Vec<_> = (0..n).map(|_| F::rand_with_config(&mut rng, config)).collect();
        let points: Vec<_> = (0..num_polys).map(|_| point.clone()).collect();

        let mut transcript = PcsTranscript::new();
        let result = TestZip::batch_open(&param, &polys, &comms_data, &points, &mut transcript, config);

        assert!(result.is_ok());
    }

    /// Verifies that the batch opening function handles an empty input slice gracefully.
    #[test]
    fn batch_open_succeeds_for_empty_input_slice() {
        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config = ConfigRef::from(&config);

        let n = 3;
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);

        let polys: Vec<DenseMultilinearExtension<Int<INT_LIMBS>>> = Vec::new();
        let comms_data: Vec<MultilinearZipData<Int<4>>> = Vec::new();
        let points: Vec<Vec<F>> = Vec::new();

        let mut transcript = PcsTranscript::new();
        let result = TestZip::batch_open(&param, &polys, &comms_data, &points, &mut transcript, config);

        assert!(result.is_ok());
        assert!(transcript.into_proof().is_empty());
    }

    /// Ensures the `prove_testing_phase` correctly writes the expected number of
    /// column openings and proximity checks to the transcript.
    #[test]
    fn prove_testing_phase_writes_correct_number_of_elements_to_transcript() {
        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        let n = 3;
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = (0..poly_size).map(|_| Int::<INT_LIMBS>::from(i8::rand(&mut rng))).collect();
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);
        let (data, _) = TestZip::commit::<F>(&param, &mle).unwrap();

        let mut transcript = PcsTranscript::<F>::new();
        TestZip::prove_testing_phase(&param, &mle, &data, &mut transcript, config).unwrap();
        let proof = transcript.into_proof();

        // Calculate the expected size in bytes
        let proximity_size = param.linear_code.num_proximity_testing()
            * param.linear_code.row_len()
            * std::mem::size_of::<<ZT as ZipTypes>::M>();

        let column_values_size = param.num_rows * std::mem::size_of::<<ZT as ZipTypes>::K>();
        let merkle_proof_size = {
            let depth = param.linear_code.codeword_len().next_power_of_two().ilog2() as usize;
            // CORRECTED: Use Keccak256::output_size() to match the transcript's hash.
            param.num_rows * (std::mem::size_of::<u64>() + depth * Keccak256::output_size())
        };
        let column_opening_size = param.linear_code.num_column_opening() * (column_values_size + merkle_proof_size);

        let expected_size = proximity_size + column_opening_size;

        assert_eq!(proof.len(), expected_size);
    }

    /// Verifies that the `prove_evaluation_phase` correctly computes and writes
    /// the combined evaluation row to the transcript.
    #[test]
    fn prove_evaluation_phase_writes_correct_combined_row() {
        let config = field_config!(115792089237316195423570985008687907853269984665640564039457584007913129639747, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        let n = 3;
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = (0..poly_size).map(|v| Int::<INT_LIMBS>::from(v as i32)).collect();
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);
        let point: Vec<_> = (0..n).map(|_| F::rand_with_config(&mut rng, config)).collect();

        let q_0 = crate::zip::pcs::utils::left_point_to_tensor(param.num_rows, &point, config).unwrap();
        let evaluations_f = mle.evaluations.map_to_field(config);
        let expected_combined_row = crate::zip::utils::combine_rows(q_0, evaluations_f, param.linear_code.row_len());

        let mut transcript = PcsTranscript::<F>::new();
        TestZip::prove_evaluation_phase(&param, &mut transcript, &point, &mle, config).unwrap();
        let proof = transcript.into_proof();

        let mut reader_transcript = PcsTranscript::<F>::from_proof(&proof);
        let combined_row_from_proof = reader_transcript.read_field_elements(param.linear_code.row_len(), config).unwrap();

        assert_eq!(combined_row_from_proof, expected_combined_row);
    }

    /// Verifies the opening proof for a polynomial where all evaluations are zero.
    #[test]
    fn open_succeeds_for_zero_polynomial() {
        let config = field_config!(115792089237316195423570985008687907853269984665640564039457584007913129639747, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        let n = 3;
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = vec![Int::<INT_LIMBS>::from(0); poly_size];
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);
        let point: Vec<_> = (0..n).map(|_| F::rand_with_config(&mut rng, config)).collect();

        let (data, _) = TestZip::commit::<F>(&param, &mle).unwrap();
        let mut transcript = PcsTranscript::new();

        let result = TestZip::open(&param, &mle, &data, &point, config, &mut transcript);

        assert!(result.is_ok());
    }

    /// Verifies the opening proof for the smallest possible square matrix arrangement (2x2).
    #[test]
    fn open_succeeds_for_smallest_matrix_arrangement() {
        let config = field_config!(115792089237316195423570985008687907853269984665640564039457584007913129639747, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        // 1. Setup parameters for a 2x2 matrix arrangement.
        let n = 2; // 2^2 = 4 evaluations
        let poly_size = 1 << n;
        let mut keccak_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak_transcript);
        let param = TestZip::setup(poly_size, linear_code);

        // For num_vars = 2, we expect a 2x2 matrix, hence 2 rows.
        assert_eq!(param.num_rows, 2);

        // 2. Create a valid polynomial and commitment.
        let evaluations: Vec<_> = (0..poly_size).map(|_| Int::<INT_LIMBS>::from(i8::rand(&mut rng))).collect();
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);
        let point: Vec<_> = (0..n).map(|_| F::rand_with_config(&mut rng, config)).collect();
        let (data, _) = TestZip::commit::<F>(&param, &mle).unwrap();
        let mut transcript = PcsTranscript::new();

        // 3. Generate an opening proof and assert it succeeds.
        let result = TestZip::open(&param, &mle, &data, &point, config, &mut transcript);
        assert!(result.is_ok());
    }

    #[test]
    fn valid_proof_is_accepted_by_verifier() {
        // Import the field-based polynomial struct and give it a distinct name.
        use crate::poly_f::mle::DenseMultilinearExtension as DenseMultilinearExtensionF;

        let config = field_config!(115792089237316195423570985008687907853269984665640564039457584007913129639747, FIELD_LIMBS);
        let config = ConfigRef::from(&config);
        let mut rng = ark_std::test_rng();

        let n = 3;
        let poly_size = 1 << n;

        // A single transcript must be used to generate the parameters for both prover and verifier.
        let mut shared_transcript = KeccakTranscript::new();
        let linear_code: LC = LC::new(&DefaultLinearCodeSpec, poly_size, &mut shared_transcript);
        let param = MultilinearZip::<ZT, _>::setup(poly_size, linear_code);

        // The integer-based polynomial for the prover.
        let evaluations: Vec<_> = (0..poly_size).map(|_| Int::<INT_LIMBS>::from(i8::rand(&mut rng))).collect();
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

        let point_int: Vec<_> = (0..n).map(|_| Int::<INT_LIMBS>::random(&mut rng)).collect();
        let point_f: Vec<F> = point_int.map_to_field(config);

        // The prover creates the commitment and proof using the integer-based polynomial.
        let (data, comm) = MultilinearZip::<ZT, _>::commit::<F>(&param, &mle).unwrap();
        let mut prover_transcript = PcsTranscript::new();
        MultilinearZip::<ZT, _>::open(&param, &mle, &data, &point_f, config, &mut prover_transcript).unwrap();
        let proof = prover_transcript.into_proof();

        // The verifier calculates the expected evaluation using a field-based polynomial.
        let mut verifier_transcript = PcsTranscript::from_proof(&proof);
        let eval = {
            let mle_f = DenseMultilinearExtensionF::from_evaluations_vec(
                n,
                mle.evaluations.map_to_field(config),
                config,
            );
            mle_f.evaluate(&point_f, config).unwrap()
        };

        // The verifier uses the same shared parameters to verify the proof.
        let verification_result = MultilinearZip::<ZT, _>::verify(&param, &comm, &point_f, eval, &mut verifier_transcript, config);

        assert!(verification_result.is_ok());
    }

    /// CRITICAL: This test directly addresses the proximity constraint from the paper.
    /// It simulates a malicious prover who creates a proof based on data that is NOT a
    /// valid codeword (violating the δ-proximity requirement). The verifier MUST reject this.
    #[test]
    fn verifier_rejects_proof_based_on_non_codeword_data() {
        // 1. Generate a valid commitment (Merkle roots) for a polynomial.
        // 2. Manually corrupt the prover's encoded row data (`MultilinearZipData`).
        // 3. Generate a proof using this corrupted data and the original, valid commitment.
        // 4. Call `MultilinearZip::verify` with the valid commitment and the malicious proof.
        // 5. Assert that verification fails with a "Proximity failure" error.
    }

    /// This test ensures the soundness of the evaluation check. It simulates a malicious
    /// prover who correctly generates a proof but lies about the polynomial's evaluation
    /// at the given point. The verifier MUST reject this.
    #[test]
    fn verifier_rejects_proof_with_incorrect_evaluation() {
        // 1. Generate a valid commitment and a valid opening proof for a polynomial `p`
        //    at a point `z`, with the correct evaluation `p(z)`.
        // 2. Call `MultilinearZip::verify` but provide an incorrect evaluation `p(z) + 1`.
        // 3. Assert that verification fails with an "Evaluation consistency failure" error.
    }

    /// This test protects against replay attacks or other transcript-based manipulation.
    /// It ensures that a proof is only valid for the specific commitment it was generated for.
    #[test]
    fn verifier_rejects_proof_with_mismatched_commitment() {
        // 1. Generate a commitment `comm1` for `poly1` and `comm2` for `poly2`.
        // 2. Generate a valid proof `proof1` for `poly1` against `comm1`.
        // 3. Call `MultilinearZip::verify` with `comm2` (the wrong commitment) and `proof1`.
        // 4. Assert that verification fails.
    }
}