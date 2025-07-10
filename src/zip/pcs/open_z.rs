#![allow(non_snake_case)]
use crate::field_config::ConfigRef;
use ark_std::borrow::Cow;
use ark_std::vec::Vec;

use ark_std::iterable::Iterable;

use itertools::izip;

use crate::{
    field::{conversion::FieldMap, RandomField},
    poly_z::mle::DenseMultilinearExtension,
    zip::{
        code::{LinearCodes, Zip, ZipSpec},
        pcs_transcript::PcsTranscript,
        utils::{combine_rows, expand},
        Error,
    },
};

use super::{
    structs::{MultilinearZip, MultilinearZipData, ZipTranscript},
    utils::{left_point_to_tensor, validate_input, ColumnOpening},
};

impl<const I: usize, const L: usize, const K: usize, const M: usize, S, T>
    MultilinearZip<I, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
{
    pub fn open<'cfg, const N: usize>(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        commit_data: &Self::Data,
        point: &[RandomField<'cfg, N>],
        field: ConfigRef<'cfg, N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        validate_input("open", pp.num_vars(), [poly], [point])?;

        Self::prove_testing_phase(pp, poly, commit_data, transcript, field)?;

        Self::prove_evaluation_phase(pp, transcript, point, poly, field)?;

        Ok(())
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30
    pub fn batch_open<'cfg, 'a, const N: usize>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension<I>>,
        comms: impl Iterable<Item = &'a MultilinearZipData<I, K>>,
        points: &[Vec<RandomField<'cfg, N>>],
        transcript: &mut PcsTranscript<N>,
        field: ConfigRef<'cfg, N>,
    ) -> Result<(), Error> {
        for (poly, comm, point) in izip!(polys.iter(), comms.iter(), points.iter()) {
            Self::open(pp, poly, comm, point, field, transcript)?;
        }
        Ok(())
    }

    // Subprotocol functions
    fn prove_evaluation_phase<'cfg, const N: usize>(
        pp: &Self::ProverParam,
        transcript: &mut PcsTranscript<N>,
        point: &[RandomField<'cfg, N>],
        poly: &Self::Polynomial,
        field: ConfigRef<'cfg, N>,
    ) -> Result<(), Error> {
        let num_rows = pp.num_rows();
        let row_len = <Zip<I, L> as LinearCodes<I, L>>::row_len(pp.zip());

        // We prove evaluations over the field,so integers need to be mapped to field elements first
        let q_0 = left_point_to_tensor(num_rows, point, field).unwrap();

        let evaluations = poly.evaluations.map_to_field(field);

        let q_0_combined_row = if num_rows > 1 {
            // Return the evaluation row combination
            let combined_row = combine_rows(q_0, evaluations, row_len);
            Cow::<Vec<RandomField<N>>>::Owned(combined_row)
        } else {
            // If there is only one row, we have no need to take linear combinations
            // We just return the evaluation row combination
            Cow::Borrowed(&evaluations)
        };

        transcript.write_field_elements(&q_0_combined_row)
    }

    pub(super) fn prove_testing_phase<const N: usize>(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        commit_data: &Self::Data,
        transcript: &mut PcsTranscript<N>,
        field: ConfigRef<N>, // This is only needed to called the transcript but we are getting integers not fields
    ) -> Result<(), Error> {
        if pp.num_rows() > 1 {
            // If we can take linear combinations
            // perform the proximity test an arbitrary number of times
            for _ in 0..<Zip<I, L> as LinearCodes<I, L>>::num_proximity_testing(pp.zip()) {
                let coeffs = transcript
                    .fs_transcript
                    .get_integer_challenges::<I>(pp.num_rows());
                let coeffs = coeffs.iter().map(expand::<I, M>);

                let evals = poly.evaluations.iter().map(expand::<I, M>);
                let combined_row = combine_rows(
                    coeffs,
                    evals,
                    <Zip<I, L> as LinearCodes<I, L>>::row_len(pp.zip()),
                );

                transcript.write_integers(&combined_row)?;
            }
        }

        // Open merkle tree for each column drawn
        for _ in 0..<Zip<I, L> as LinearCodes<I, L>>::num_column_opening(pp.zip()) {
            let column = transcript.squeeze_challenge_idx(
                field,
                <Zip<I, L> as LinearCodes<I, L>>::codeword_len(pp.zip()),
            );
            Self::open_merkle_trees_for_column(pp, commit_data, column, transcript)?;
        }
        Ok(())
    }

    pub(super) fn open_merkle_trees_for_column<const N: usize>(
        pp: &Self::ProverParam,
        commit_data: &MultilinearZipData<I, K>,
        column: usize,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        //Write the elements in the squeezed column to the shared transcript
        transcript.write_integers(
            &commit_data
                .rows()
                .iter()
                .copied()
                .skip(column)
                .step_by(<Zip<I, L> as LinearCodes<I, L>>::codeword_len(pp.zip()))
                .collect::<Vec<_>>(),
        )?;
        ColumnOpening::open_at_column(column, commit_data, transcript)
            .map_err(|_| Error::InvalidPcsOpen("Failed to open merkle tree".into()))?;
        Ok(())
    }
}
