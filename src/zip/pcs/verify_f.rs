use ark_ff::Zero;
use ark_std::iterable::Iterable;

use crate::{
    field::RandomField as F,
    field_config::FieldConfig,
    zip::{
        code::{LinearCodes, ZipSpec},
        pcs_transcript::PcsTranscript,
        utils::inner_product,
        Error,
    },
};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment},
    utils::point_to_tensor_f,
};

impl<const N: usize, S> MultilinearZip<N, S>
where
    S: ZipSpec,
{
    pub fn verify_f(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &[F<N>],
        eval: &F<N>,
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        // TODO bring this back
        // validate_input("verify", vp.num_vars(), [], [point])?;

        let row_len = vp.zip().row_len();
        let codeword_len = vp.zip().codeword_len();

        // Retrieve the row combinations from the transcript
        // Pair the combinations with the coefficients that generated them
        let mut combined_rows = Vec::with_capacity(vp.zip().num_proximity_testing() + 1);

        if vp.num_rows() > 1 {
            for _ in 0..vp.zip().num_proximity_testing() {
                let coeffs = transcript
                    .fs_transcript
                    .get_challenges(vp.num_rows(), field);
                let mut combined_row = transcript.read_integers(row_len)?;

                combined_row.resize(codeword_len, 0i64);
                vp.zip().encode(&combined_row);
                let combined_row_f = combined_row
                    .iter()
                    .map(|i| F::from_i64(*i, field).unwrap())
                    .collect::<Vec<_>>();
                combined_rows.push((coeffs, combined_row_f));
            }
        }

        let (t_0, t_1) = point_to_tensor_f(vp.num_rows(), point, field)?;

        combined_rows.push({
            let mut t_0_combined_row = transcript.read_field_elements(row_len, field)?;

            t_0_combined_row.resize(codeword_len, F::zero());

            (t_0, t_0_combined_row)
        });

        let depth = codeword_len.next_power_of_two().ilog2() as usize;

        // Ensure that the test combinations are valid codewords
        for _ in 0..vp.zip().num_column_opening() {
            let column = transcript.squeeze_challenge_idx(field, codeword_len);

            let items = transcript.read_I256_vec(vp.num_rows())?;

            let merkle_path = transcript.read_commitments(depth)?;

            Self::verify_proximity(&combined_rows, &items, column, vp.num_rows(), field)?;

            Self::verify_merkle_path(&items, &merkle_path, column, comm)?;
        }

        // verify consistency
        let t_0_combined_row = combined_rows
            .last()
            .map(|(_, combined_row)| &combined_row[..row_len])
            .unwrap();

        if inner_product(t_0_combined_row, &t_1) != *eval {
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
}
