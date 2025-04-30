#![allow(non_snake_case)]
use std::borrow::Cow;

use ark_std::iterable::Iterable;
use crypto_bigint::Int;
use i256::I256;
use itertools::izip;

use crate::{
    field::{conversion::FieldMap, RandomField as F},
    field_config::FieldConfig,
    poly_z::mle::DenseMultilinearExtension,
    zip::{
        code::{LinearCodes, ZipSpec},
        pcs_transcript::PcsTranscript,
        utils::combine_rows,
        Error,
    },
};

use super::{
    structs::{MultilinearZip, MultilinearZipData, ZipTranscript},
    utils::{left_point_to_tensor, validate_input, ColumnOpening},
};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub fn open(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        commit_data: &Self::Data,
        point: &[F<N>],
        field: *const FieldConfig<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        validate_input("open", pp.num_vars(), [poly], [point])?;

        Self::prove_testing_phase(pp, poly, commit_data, transcript, field)?;

        Self::prove_evaluation_phase(pp, transcript, point, poly, field)?;

        Ok(())
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30
    pub fn batch_open<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension<N>>,
        comms: impl Iterable<Item = &'a MultilinearZipData<N>>,
        points: &[Vec<F<N>>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let mut proofs = vec![];
        for (poly, comm, point) in izip!(polys.iter(), comms.iter(), points.iter()) {
            proofs.push(Self::open(pp, poly, comm, point, field, transcript)?);
        }
        Ok(())
    }

    // Subprotocol functions
    fn prove_evaluation_phase(
        pp: &Self::ProverParam,
        transcript: &mut PcsTranscript<N>,
        point: &[F<N>],
        poly: &Self::Polynomial,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let num_rows = pp.num_rows();
        let row_len = pp.zip().row_len();

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
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        commit_data: &Self::Data,
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>, // This is only needed to called the transcript but we are getting integers not fields
    ) -> Result<(), Error> {
        if pp.num_rows() > 1 {
            // If we can take linear combinations
            // perform the proximity test an arbitrary number of times
            for _ in 0..pp.zip().num_proximity_testing() {
                let coeffs = transcript
                    .fs_transcript
                    .get_integer_challenges(pp.num_rows());
                let coeffs = coeffs.iter().map(|x| I256::from(*x));
                let evals = poly.evaluations.iter().map(|x| Int::<N>::from(*x));
                let combined_row = combine_rows(coeffs, evals, pp.zip().row_len());

                transcript.write_I256_vec(&combined_row)?;
            }
        }

        // Open merkle tree for each column drawn
        for _ in 0..pp.zip().num_column_opening() {
            let column = transcript.squeeze_challenge_idx(field, pp.zip().codeword_len());
            Self::open_merkle_trees_for_column(pp, commit_data, column, transcript)?;
        }
        Ok(())
    }

    pub(super) fn open_merkle_trees_for_column(
        pp: &Self::ProverParam,
        commit_data: &MultilinearZipData<N>,
        column: usize,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        //Write the elements in the squeezed column to the shared transcript
        transcript.write_I512_vec(
            &commit_data
                .rows()
                .iter()
                .copied()
                .skip(column)
                .step_by(pp.zip().codeword_len())
                .collect::<Vec<_>>(),
        )?;
        ColumnOpening::open_at_column(column, commit_data, transcript)
            .map_err(|_| Error::InvalidPcsOpen("Failed to open merkle tree".to_string()))?;
        Ok(())
    }
}
