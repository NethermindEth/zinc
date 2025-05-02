use ark_std::iterable::Iterable;

use crypto_bigint::Int;

use sha3::{digest::Output, Keccak256};

use crate::{
    field::{conversion::FieldMap, RandomField as F},
    field_config::FieldConfig,
    zip::{
        code::{LinearCodes, Zip, ZipSpec},
        pcs_transcript::PcsTranscript,
        utils::{expand, inner_product},
        Error,
    },
};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment, ZipTranscript},
    utils::{point_to_tensor, validate_input, ColumnOpening},
};

impl<const N: usize, const L: usize, const K: usize, const M: usize, S, T>
    MultilinearZip<N, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
{
    pub fn verify(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &[F<N>],
        eval: F<N>,
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        validate_input::<N>("verify", vp.num_vars(), [], [point])?;

        let columns_opened = Self::verify_testing(vp, comm.roots(), transcript, field)?;

        Self::verify_evaluation_z(vp, point, eval, &columns_opened, transcript, field)?;

        Ok(())
    }

    pub fn batch_verify_z<'a>(
        vp: &Self::VerifierParam,
        comms: impl Iterable<Item = &'a MultilinearZipCommitment<N>>,
        points: &[Vec<F<N>>],
        evals: &[F<N>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        for (i, (eval, comm)) in evals.iter().zip(comms.iter()).enumerate() {
            Self::verify(vp, comm, &points[i], *eval, transcript, field)?;
        }
        Ok(())
    }

    pub(super) fn verify_testing(
        vp: &Self::VerifierParam,
        roots: &[Output<Keccak256>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<Vec<(usize, Vec<Int<K>>)>, Error> {
        // Gather the coeffs and encoded combined rows per proximity test
        let mut encoded_combined_rows = Vec::with_capacity(
            <Zip<N, L> as LinearCodes<N, L>>::num_proximity_testing(vp.zip()),
        );
        if vp.num_rows() > 1 {
            for _ in 0..<Zip<N, L> as LinearCodes<N, L>>::num_proximity_testing(vp.zip()) {
                let coeffs = transcript
                    .fs_transcript
                    .get_integer_challenges::<N>(vp.num_rows());

                let combined_row: Vec<Int<M>> = transcript
                    .read_integers(<Zip<N, L> as LinearCodes<N, L>>::row_len(vp.zip()))?;

                let encoded_combined_row: Vec<Int<M>> = vp.zip().encode_wide(&combined_row);
                encoded_combined_rows.push((coeffs, encoded_combined_row));
            }
        }

        let mut columns_opened: Vec<(usize, Vec<Int<K>>)> = Vec::with_capacity(
            <Zip<N, L> as LinearCodes<N, L>>::num_column_opening(vp.zip()),
        );
        for _ in 0..<Zip<N, L> as LinearCodes<N, L>>::num_column_opening(vp.zip()) {
            let column_idx = transcript.squeeze_challenge_idx(
                field,
                <Zip<N, L> as LinearCodes<N, L>>::codeword_len(vp.zip()),
            );
            let column_values = transcript.read_integers(vp.num_rows())?;

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
        coeffs: &[Int<N>],
        encoded_combined_row: &[Int<M>],
        column_entries: &[Int<K>],
        column: usize,
        num_rows: usize,
    ) -> Result<(), Error> {
        let column_entries_comb = if num_rows > 1 {
            let coeffs: Vec<_> = coeffs.iter().map(expand::<N, M>).collect();
            let column_entries: Vec<_> = column_entries.iter().map(expand::<K, M>).collect();
            inner_product(coeffs.iter(), column_entries.iter())
        } else {
            expand::<K, M>(&column_entries[0])
        };

        if column_entries_comb != encoded_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
        }
        Ok(())
    }

    fn verify_evaluation_z(
        vp: &Self::VerifierParam,
        point: &[F<N>],
        eval: F<N>,
        columns_opened: &[(usize, Vec<Int<K>>)],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let q_0_combined_row = transcript
            .read_field_elements(<Zip<N, L> as LinearCodes<N, L>>::row_len(vp.zip()), field)?;
        let encoded_combined_row = vp.zip().encode_f(&q_0_combined_row, field);

        let (q_0, q_1) = point_to_tensor(vp.num_rows(), point, field)?;

        if inner_product(&q_0_combined_row, &q_1) != eval {
            return Err(Error::InvalidPcsOpen(
                "Evaluation consistency failure".to_string(),
            ));
        }

        for (column_idx, column_values) in columns_opened.iter() {
            Self::verify_proximity_q_0(
                &q_0,
                &encoded_combined_row,
                column_values,
                *column_idx,
                vp.num_rows(),
                field,
            )?;
        }

        Ok(())
    }

    fn verify_proximity_q_0(
        q_0: &Vec<F<N>>,
        q_0_combined_row: &[F<N>],
        column_entries: &[Int<K>],
        column: usize,
        num_rows: usize,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let column_entries_comb = if num_rows > 1 {
            let column_entries = column_entries.map_to_field(field);
            inner_product(q_0, &column_entries)
            // TODO: this inner product is taking a long time.
        } else {
            column_entries.first().unwrap().map_to_field(field)
        };
        if column_entries_comb != q_0_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
        }

        Ok(())
    }
}
