use ark_std::iterable::Iterable;

use i256::{I256, I512};

use crate::{
    field::RandomField as F,
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
    utils::point_to_tensor_f,
};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub fn verify_f(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &[F<N>],
        eval_f: &F<N>,
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        // validate_input("verify", vp.num_vars(), [], [point])?;

        let row_len = vp.zip().row_len();
        let codeword_len = vp.zip().codeword_len();

        // Retrieve the row combinations from the transcript
        // Pair the combinations with the coefficients that generated them
        let mut encoded_combined_rows = Vec::with_capacity(vp.zip().num_proximity_testing());

        if vp.num_rows() > 1 {
            for _ in 0..vp.zip().num_proximity_testing() {
                let coeffs: Vec<_> = transcript
                    .fs_transcript
                    .get_integer_challenges(vp.num_rows())
                    .iter()
                    .map(|i| I256::from(*i))
                    .collect();

                let combined_row = transcript.read_I256_vec(row_len)?;

                let code = vp.zip().encode(&combined_row);

                encoded_combined_rows.push((coeffs, code));
            }
        }

        let depth = codeword_len.next_power_of_two().ilog2() as usize;
        let t_0_combined_row = transcript.read_field_elements(row_len, field)?;
        let (t_0_f, t_1_f) = point_to_tensor_f(vp.num_rows(), point, field)?;

        // Ensure that the test combinations are valid codewords
        for _ in 0..vp.zip().num_column_opening() {
            let column = transcript.squeeze_challenge_idx(field, codeword_len);

            let items = transcript.read_I512_vec(vp.num_rows())?;

            let merkle_path = transcript.read_commitments(depth)?;

            Self::verify_proximity_f(&encoded_combined_rows, &items, column, vp.num_rows())?;
            Self::verify_proximity_t_0_f(
                &t_0_f,
                &vp.zip().encode_f(&t_0_combined_row, field),
                &items,
                column,
                vp.num_rows(),
                field,
            )?;
            Self::verify_merkle_path(&items, &merkle_path, column, comm)?;
        }

        // verify consistency

        if inner_product(&t_0_combined_row, &t_1_f) != *eval_f {
            return Err(Error::InvalidPcsOpen("Consistency failure".to_string()));
        }

        Ok(())
    }

    pub fn batch_verify_f<'a>(
        vp: &Self::VerifierParam,
        comms: impl Iterable<Item = &'a MultilinearZipCommitment<N>>,
        points: &[Vec<F<N>>],
        evals: &[F<N>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        for (i, (eval, comm)) in evals.iter().zip(comms.iter()).enumerate() {
            Self::verify_f(vp, comm, &points[i], eval, transcript, field)?;
        }
        Ok(())
    }

    pub(super) fn verify_proximity_f(
        combined_rows: &[(Vec<I256>, Vec<I512>)],
        column_entries: &[I512],
        column: usize,
        num_rows: usize,
    ) -> Result<(), Error> {
        for (coeff, encoded) in combined_rows.iter() {
            let column_entries_comb = if num_rows > 1 {
                let coeff: Vec<_> = coeff.iter().map(|i| I256_to_I512(*i)).collect();

                inner_product(&coeff, column_entries)
            } else {
                column_entries[0]
            };

            if column_entries_comb != encoded[column] {
                return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
            }
        }
        Ok(())
    }

    fn verify_proximity_t_0_f(
        t_0_f: &Vec<F<N>>,
        t_0_combined_row: &[F<N>],
        column_entries: &[I512],
        column: usize,
        num_rows: usize,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let column_entries_comb = if num_rows > 1 {
            let column_entries = column_entries
                .iter()
                .map(|i| F::from_I512(*i, field))
                .collect::<Vec<_>>();
            inner_product(t_0_f, &column_entries)
        } else {
            F::from_I512(column_entries[0], field)
        };
        if column_entries_comb != t_0_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
        }

        Ok(())
    }
}
