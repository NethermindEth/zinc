#![allow(non_snake_case)]
use std::borrow::Cow;

use ark_std::iterable::Iterable;

use itertools::izip;

use crate::{
    field::{conversion::FieldMap, RandomField as F},
    field_config::FieldConfig,
    poly_z::mle::DenseMultilinearExtension,
    zip::{
        code::{LinearCodes, Zip, ZipSpec},
        pcs_transcript::PcsTranscript,
        utils::{combine_rows, expand},
        ZipError,
    },
};

use super::{
    structs::{MultilinearZip, MultilinearZipData, ZipTranscript},
    utils::{left_point_to_tensor, validate_input, ColumnOpening},
};

impl<const N: usize, const L: usize, const K: usize, const M: usize, S, T>
    MultilinearZip<N, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
{
    /// Proves an evaluation of the polynomial at a given point.
    ///
    /// This function generates a proof that a polynomial has a specific value
    /// at a given evaluation point. The proof can later be used to verify
    /// the evaluation without revealing the entire polynomial.
    ///
    /// # Arguments
    ///
    /// * `pp` - A reference to the prover parameters. e.g. the number of variables expected
    /// * `poly` - A reference to the polynomial for which the evaluation proof
    ///            is being generated.
    /// * `commit_data` - The data generated during the commitment phase.
    ///                   i.e. The EA code matrix and the merkle tree(s) with the rows of th
    ///                   matrix as leaves.  
    /// * `point` - A slice representing the point at which the polynomial is evaluated.
    ///            
    /// * `field` - A raw pointer to the random field in which we will prove the final evaluation.
    /// * `transcript` - A mutable reference to the transcript used for recording the proof.
    ///
    /// # Returns
    ///
    /// Returns a `Result<(), Error>`. If the proof is generated successfully,
    /// the function returns `Ok(())`. Otherwise, it returns an `Error`
    /// describing the failure reason.
    ///
    pub fn open(
        pp: &Self::Param,
        poly: &Self::Polynomial,
        commit_data: &Self::Data,
        point: &[F<N>],
        field: *const FieldConfig<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), ZipError> {
        validate_input("open", pp.num_vars(), [poly], [point])?;

        Self::prove_testing_phase(pp, poly, commit_data, transcript, field)?;

        Self::prove_evaluation_phase(pp, transcript, point, poly, field)?;

        Ok(())
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30ynomal
    /// Open multiple polynomials at (possibly) different point
    pub fn batch_open<'a>(
        pp: &Self::Param,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension<N>>,
        comms: impl Iterable<Item = &'a MultilinearZipData<N, K>>,
        points: &[Vec<F<N>>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), ZipError> {
        let mut proofs = vec![];
        for (poly, comm, point) in izip!(polys.iter(), comms.iter(), points.iter()) {
            proofs.push(Self::open(pp, poly, comm, point, field, transcript)?);
        }
        Ok(())
    }

    // Subprotocol functions
    fn prove_evaluation_phase(
        pp: &Self::Param,
        transcript: &mut PcsTranscript<N>,
        point: &[F<N>],
        poly: &Self::Polynomial,
        field: *const FieldConfig<N>,
    ) -> Result<(), ZipError> {
        let num_rows = pp.num_rows();
        let row_len = <Zip<N, L> as LinearCodes<N, L>>::row_len(pp.zip());

        // We prove evaluations over the field,so integers need to be mapped to field elements first
        let q_0 = left_point_to_tensor(num_rows, point, field).unwrap();

        let evaluations = poly.evaluations.map_to_field(field);

        let q_0_combined_row = if num_rows > 1 {
            // Return the evaluation row combination
            let combined_row = combine_rows(q_0, evaluations, row_len);
            Cow::<Vec<F<N>>>::Owned(combined_row)
        } else {
            // If there is only one row, we have no need to take linear combinations
            // We just return the evaluation row combination
            Cow::Borrowed(&evaluations)
        };

        transcript.write_field_elements(&q_0_combined_row)
    }

    pub(super) fn prove_testing_phase(
        pp: &Self::Param,
        poly: &Self::Polynomial,
        commit_data: &Self::Data,
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>, // This is only needed to called the transcript but we are getting integers not fields
    ) -> Result<(), ZipError> {
        if pp.num_rows() > 1 {
            // If we can take linear combinations
            // perform the proximity test an arbitrary number of times
            for _ in 0..<Zip<N, L> as LinearCodes<N, L>>::num_proximity_testing(pp.zip()) {
                let coeffs = transcript
                    .fs_transcript
                    .get_integer_challenges::<N>(pp.num_rows());
                let coeffs = coeffs.iter().map(expand::<N, M>);

                let evals = poly.evaluations.iter().map(expand::<N, M>);
                let combined_row = combine_rows(
                    coeffs,
                    evals,
                    <Zip<N, L> as LinearCodes<N, L>>::row_len(pp.zip()),
                );

                transcript.write_integers(&combined_row)?;
            }
        }

        // Open merkle tree for each column drawn
        for _ in 0..<Zip<N, L> as LinearCodes<N, L>>::num_column_opening(pp.zip()) {
            let column = transcript.squeeze_challenge_idx(
                field,
                <Zip<N, L> as LinearCodes<N, L>>::codeword_len(pp.zip()),
            );
            Self::open_merkle_trees_for_column(pp, commit_data, column, transcript)?;
        }
        Ok(())
    }

    pub(super) fn open_merkle_trees_for_column(
        pp: &Self::Param,
        commit_data: &MultilinearZipData<N, K>,
        column: usize,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), ZipError> {
        //Write the elements in the squeezed column to the shared transcript
        transcript.write_integers(
            &commit_data
                .rows()
                .iter()
                .copied()
                .skip(column)
                .step_by(<Zip<N, L> as LinearCodes<N, L>>::codeword_len(pp.zip()))
                .collect::<Vec<_>>(),
        )?;
        ColumnOpening::open_at_column(column, commit_data, transcript)
            .map_err(|_| ZipError::InvalidPcsOpen("Failed to open merkle tree".to_string()))?;
        Ok(())
    }
}
