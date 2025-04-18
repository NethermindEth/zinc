#![allow(non_snake_case)]

use std::borrow::Cow;

use ark_std::iterable::Iterable;
use i256::I256;
use itertools::izip;

use crate::{
    field::{conversion::FieldMap, RandomField},
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
    utils::{point_to_tensor_f, ColumnOpening, MerkleProof},
};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub fn open_f(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        comm: &Self::Data,
        point: &[RandomField<N>],
        field: *const FieldConfig<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        // validate_input("open", pp.num_vars(), [poly], [point])?;

        Self::prove_test(pp, poly, comm, transcript, field)?;

        let row_len = pp.zip().row_len();
        Self::prove_evaluation_f(pp.num_rows(), row_len, transcript, point, poly, field)?;

        Ok(())
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30
    pub fn batch_open_f<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension>,
        comms: impl Iterable<Item = &'a MultilinearZipData<N>>,
        points: &[Vec<RandomField<N>>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        //	use std::env;
        //	let key = "RAYON_NUM_THREADS";
        //	env::set_var(key, "8");
        for (poly, comm, point) in izip!(polys.iter(), comms.iter(), points.iter()) {
            Self::open_f(pp, poly, comm, point, field, transcript)?;
        }
        Ok(())
    }

    pub(super) fn prove_test(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        commitment_data: &Self::Data,
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>, // This is only needed to called the trasncript but we are getting integers not fields
    ) -> Result<(), Error> {
        if pp.num_rows() > 1 {
            // If we can take linear combinations
            // perform the proximity test an arbitrary number of times
            for _ in 0..pp.zip().num_proximity_testing() {
                let coeffs = transcript
                    .fs_transcript
                    .get_integer_challenges(pp.num_rows());
                let coeffs = coeffs.iter().map(|x| I256::from(*x));
                let evals = poly.evaluations.iter().map(|x| I256::from(*x));
                let combined_row = combine_rows(coeffs, evals, pp.zip().row_len());

                transcript.write_I256_vec(&combined_row)?;
            }
        }

        // Open merkle tree for each column drawn
        for _ in 0..pp.zip().num_column_opening() {
            let column = transcript.squeeze_challenge_idx(field, pp.zip().codeword_len());
            Self::open_merkle_trees_for_column(pp, commitment_data, column, transcript)?;
        }
        Ok(())
    }

    pub(super) fn open_merkle_trees_for_column(
        pp: &Self::ProverParam,
        comm: &MultilinearZipData<N>,
        column: usize,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        //Write the elements in the squeezed column to the shared transcript
        transcript.write_I512_vec(
            &comm
                .rows()
                .iter()
                .copied()
                .skip(column)
                .step_by(pp.zip().codeword_len())
                .collect::<Vec<_>>(),
        )?;

        let merkle_depth = pp.zip().row_len().next_power_of_two().ilog2() as usize;
        ColumnOpening::open_at_column(merkle_depth, column, comm, transcript);
        Ok(())
    }
    // Subprotocol functions
    fn prove_evaluation_f(
        num_rows: usize,
        row_len: usize,
        transcript: &mut PcsTranscript<N>,

        point: &[RandomField<N>],
        poly: &Self::Polynomial,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let (t_0_f, _) = point_to_tensor_f(num_rows, point, field).unwrap();

        let evaluations: Vec<RandomField<N>> = poly
            .evaluations
            .iter()
            .map(|i| i.map_to_field(field))
            .collect();

        let t_0_combined_row = if num_rows > 1 {
            // Return the evaluation row combination
            let combined_row = combine_rows(t_0_f, evaluations, row_len);
            Cow::<Vec<RandomField<N>>>::Owned(combined_row)
        } else {
            // If there is only one row, we have no need to take linear combinations
            // We just return the evaluation row combination
            Cow::Borrowed(&evaluations)
        };

        transcript.write_field_elements(&t_0_combined_row)
    }
}
