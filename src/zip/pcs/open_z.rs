#![allow(non_snake_case)]
use ark_std::{borrow::Cow, iterable::Iterable, vec::Vec};
use itertools::izip;

use super::{
    structs::{MultilinearZip, MultilinearZipData},
    utils::{left_point_to_tensor, validate_input, ColumnOpening},
};
use crate::{
    poly_z::mle::DenseMultilinearExtension,
    traits::{Field, FieldMap, Integer},
    zip::{
        code::{LinearCodes, Zip},
        pcs::structs::MultilinearZipParams,
        pcs_transcript::PcsTranscript,
        utils::{combine_rows, expand},
        Error,
    },
};

impl<I: Integer, L: Integer, K: Integer, M: Integer> MultilinearZip<I, L, K, M>
where
    L: for<'a> From<&'a I> + for<'a> From<&'a L>,
    M: for<'a> From<&'a I>,
    Zip<I, L>: LinearCodes<I, M>,
{
    pub fn open<F: Field>(
        pp: &MultilinearZipParams<I, L>,
        poly: &DenseMultilinearExtension<I>,
        commit_data: &MultilinearZipData<K>,
        point: &[F],
        field: F::R,
        transcript: &mut PcsTranscript<F>,
    ) -> Result<(), Error>
    where
        I: FieldMap<F, Output = F>,
    {
        validate_input("open", pp.num_vars, [poly], [point])?;

        Self::prove_testing_phase(pp, poly, commit_data, transcript, field)?;

        Self::prove_evaluation_phase(pp, transcript, point, poly, field)?;

        Ok(())
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30
    pub fn batch_open<'a, F: Field>(
        pp: &MultilinearZipParams<I, L>,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension<I>>,
        comms: impl Iterable<Item = &'a MultilinearZipData<K>>,
        points: &[Vec<F>],
        transcript: &mut PcsTranscript<F>,
        field: F::R,
    ) -> Result<(), Error>
    where
        I: FieldMap<F, Output = F> + 'a,
        K: 'a,
    {
        for (poly, comm, point) in izip!(polys.iter(), comms.iter(), points.iter()) {
            Self::open(pp, poly, comm, point, field, transcript)?;
        }
        Ok(())
    }

    // Subprotocol functions
    fn prove_evaluation_phase<F: Field>(
        pp: &MultilinearZipParams<I, L>,
        transcript: &mut PcsTranscript<F>,
        point: &[F],
        poly: &DenseMultilinearExtension<I>,
        field: F::R,
    ) -> Result<(), Error>
    where
        I: FieldMap<F, Output = F>,
    {
        let num_rows = pp.num_rows;
        let row_len = pp.zip.row_len();

        // We prove evaluations over the field,so integers need to be mapped to field elements first
        let q_0 = left_point_to_tensor(num_rows, point, field).unwrap();

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
        pp: &MultilinearZipParams<I, L>,
        poly: &DenseMultilinearExtension<I>,
        commit_data: &MultilinearZipData<K>,
        transcript: &mut PcsTranscript<F>,
        field: F::R, // This is only needed to call the transcript, but we are getting integers not fields
    ) -> Result<(), Error> {
        if pp.num_rows > 1 {
            // If we can take linear combinations
            // perform the proximity test an arbitrary number of times
            for _ in 0..pp.zip.num_proximity_testing() {
                let coeffs = transcript.fs_transcript.get_integer_challenges(pp.num_rows);
                let coeffs = coeffs.iter().map(expand::<I, M>);

                let evals = poly.evaluations.iter().map(expand::<I, M>);
                let combined_row = combine_rows(coeffs, evals, pp.zip.row_len());

                transcript.write_integers(combined_row.iter())?;
            }
        }

        // Open merkle tree for each column drawn
        for _ in 0..pp.zip.num_column_opening() {
            let column = transcript.squeeze_challenge_idx(field, pp.zip.codeword_len());
            Self::open_merkle_trees_for_column(pp, commit_data, column, transcript)?;
        }
        Ok(())
    }

    pub(super) fn open_merkle_trees_for_column<F: Field>(
        pp: &MultilinearZipParams<I, L>,
        commit_data: &MultilinearZipData<K>,
        column: usize,
        transcript: &mut PcsTranscript<F>,
    ) -> Result<(), Error> {
        // Write the elements in the squeezed column to the shared transcript
        transcript.write_integers(
            commit_data
                .rows
                .iter()
                .skip(column)
                .step_by(pp.zip.codeword_len()),
        )?;
        ColumnOpening::open_at_column(column, commit_data, transcript)
            .map_err(|_| Error::InvalidPcsOpen("Failed to open merkle tree".into()))?;
        Ok(())
    }
}
