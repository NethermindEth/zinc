use ark_std::iterable::Iterable;

use i256::{I256, I512};
use itertools::Itertools;
use sha3::{digest::Output, Digest, Keccak256};

use crate::{
    field::{conversion::FieldMap, RandomField},
    field_config::FieldConfig,
    zip::{
        code::{I256_to_I512, LinearCodes, ZipSpec},
        pcs_transcript::PcsTranscript,
        utils::inner_product,
        Error,
    },
};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment, ZipTranscript},
    utils::{point_to_tensor_z, validate_input, ColumnOpening},
};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub fn verify_z(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &Vec<i64>,
        eval: &i64,
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        validate_input("verify", vp.num_vars(), [], [point])?;

        let columns_opened = Self::verify_testing(vp, comm.roots(), transcript, field)?;

        Self::verify_evaluation_z(vp, point, eval, &columns_opened, transcript, field)?;

        Ok(())
    }

    pub fn batch_verify_z<'a>(
        vp: &Self::VerifierParam,
        comms: impl Iterable<Item = &'a MultilinearZipCommitment<N>>,
        points: &[Vec<i64>],
        evals: &[i64],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        for (i, (eval, comm)) in evals.iter().zip(comms.iter()).enumerate() {
            Self::verify_z(vp, comm, &points[i], eval, transcript, field)?;
        }
        Ok(())
    }

    pub(super) fn verify_testing(
        vp: &Self::VerifierParam,
        roots: &[Output<Keccak256>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<Vec<(usize, Vec<I512>)>, Error> {
        // Gather the coeffs and encoded combined rows per proximity test
        let mut encoded_combined_rows = Vec::with_capacity(vp.zip().num_proximity_testing());
        if vp.num_rows() > 1 {
            for _ in 0..vp.zip().num_proximity_testing() {
                let coeffs: Vec<_> = transcript
                    .fs_transcript
                    .get_integer_challenges(vp.num_rows())
                    .iter()
                    .map(|i| I256::from(*i))
                    .collect();

                let combined_row = transcript.read_I256_vec(vp.zip().row_len())?;

                let encoded_combined_row = vp.zip().encode(&combined_row);
                encoded_combined_rows.push((coeffs, encoded_combined_row));
            }
        }

        let mut columns_opened: Vec<(usize, Vec<I512>)> =
            Vec::with_capacity(vp.zip().num_column_opening());
        for _ in 0..vp.zip().num_column_opening() {
            let column_idx = transcript.squeeze_challenge_idx(field, vp.zip().codeword_len());
            let column_values = transcript.read_I512_vec(vp.num_rows())?;

            for (coeffs, encoded_combined_row) in encoded_combined_rows.iter() {
                Self::verify_column_testing(
                    coeffs,
                    encoded_combined_row,
                    &column_values,
                    column_idx,
                    vp.num_rows(),
                )?;
            }

            let _ = ColumnOpening::verify_column(roots, &column_values, column_idx, transcript);
            // TODO: Verify column opening is taking a long time.
            columns_opened.push((column_idx, column_values));
        }
        Ok(columns_opened)
    }

    pub(super) fn verify_column_testing(
        coeffs: &[I256],
        encoded_combined_row: &[I512],
        column_entries: &[I512],
        column: usize,
        num_rows: usize,
    ) -> Result<(), Error> {
        let column_entries_comb = if num_rows > 1 {
            let coeff: Vec<_> = coeffs.iter().map(|i| I256_to_I512(*i)).collect();
            inner_product(&coeff, column_entries)
        } else {
            column_entries[0]
        };

        if column_entries_comb != encoded_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
        }
        Ok(())
    }

    fn verify_evaluation_z(
        vp: &Self::VerifierParam,
        point: &Vec<i64>,
        eval: &i64,
        columns_opened: &[(usize, Vec<I512>)],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let t_0_combined_row = transcript.read_field_elements(vp.zip().row_len(), field)?;
        let encoded_combined_row = vp.zip().encode_f(&t_0_combined_row, field);
        let (t_0, t_1) = point_to_tensor_z(vp.num_rows(), point)?;
        let t_0_f = t_0.map_to_field(field);
        let t_1_f = t_1.map_to_field(field);
        if inner_product(&t_0_combined_row, &t_1_f) != eval.map_to_field(field) {
            return Err(Error::InvalidPcsOpen(
                "Evaluation consistency failure".to_string(),
            ));
        }

        for (column_idx, column_values) in columns_opened.iter() {
            Self::verify_proximity_t_0(
                &t_0_f,
                &encoded_combined_row,
                column_values,
                *column_idx,
                vp.num_rows(),
                field,
            )?;
        }

        Ok(())
    }

    fn verify_proximity_t_0(
        t_0_f: &Vec<RandomField<N>>,
        t_0_combined_row: &[RandomField<N>],
        column_entries: &[I512],
        column: usize,
        num_rows: usize,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let column_entries_comb = if num_rows > 1 {
            let column_entries = column_entries.map_to_field(field);
            inner_product(t_0_f, &column_entries)
            // TODO: this inner product is taking a long time.
        } else {
            column_entries.first().unwrap().map_to_field(field)
        };
        if column_entries_comb != t_0_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
        }

        Ok(())
    }
}
