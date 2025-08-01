use ark_std::{iterable::Iterable, vec::Vec};
use sha3::{digest::Output, Keccak256};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment, ZipTranscript},
    utils::{point_to_tensor, validate_input, ColumnOpening},
};
use crate::{
    traits::{CryptoInt, Field, FieldMap},
    zip::{
        code::{LinearCodes, Zip, ZipSpec},
        pcs::utils::ToBytes,
        pcs_transcript::PcsTranscript,
        utils::{expand, inner_product},
        Error,
    },
};

impl<I: CryptoInt, L: CryptoInt, K: CryptoInt + ToBytes, M: CryptoInt, S, T>
    MultilinearZip<I, L, K, M, S, T>
where
    S: ZipSpec,
    T: ZipTranscript<L>,
    L: for<'a> From<&'a I> + for<'a> From<&'a L>,
    M: for<'a> From<&'a I> + for<'a> From<&'a K> + for<'a> From<&'a L> + for<'a> From<&'a M>,
    Zip<I, L>: LinearCodes<I, M>,
{
    pub fn verify<F: Field>(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &[F],
        eval: F,
        transcript: &mut PcsTranscript<F>,
        field: F::Cr,
    ) -> Result<(), Error>
    where
        L: FieldMap<F, Output = F>,
        K: FieldMap<F, Output = F>,
    {
        validate_input::<I, F>("verify", vp.num_vars(), [], [point])?;

        let columns_opened = Self::verify_testing(vp, comm.roots(), transcript, field)?;

        Self::verify_evaluation_z(vp, point, eval, &columns_opened, transcript, field)?;

        Ok(())
    }

    pub fn batch_verify_z<'a, F: Field>(
        vp: &Self::VerifierParam,
        comms: impl Iterable<Item = &'a MultilinearZipCommitment<I>>,
        points: &[Vec<F>],
        evals: &[F],
        transcript: &mut PcsTranscript<F>,
        field: F::Cr,
    ) -> Result<(), Error>
    where
        L: FieldMap<F, Output = F>,
        K: FieldMap<F, Output = F>,
        I: 'a,
    {
        for (i, (eval, comm)) in evals.iter().zip(comms.iter()).enumerate() {
            Self::verify(vp, comm, &points[i], eval.clone(), transcript, field)?;
        }
        Ok(())
    }

    pub(super) fn verify_testing<F: Field>(
        vp: &Self::VerifierParam,
        roots: &[Output<Keccak256>],
        transcript: &mut PcsTranscript<F>,
        field: F::Cr,
    ) -> Result<Vec<(usize, Vec<K>)>, Error> {
        // Gather the coeffs and encoded combined rows per proximity test
        let mut encoded_combined_rows = Vec::with_capacity(
            <Zip<I, L> as LinearCodes<I, L>>::num_proximity_testing(vp.zip()),
        );
        if vp.num_rows() > 1 {
            for _ in 0..<Zip<I, L> as LinearCodes<I, L>>::num_proximity_testing(vp.zip()) {
                let coeffs = transcript
                    .fs_transcript
                    .get_integer_challenges(vp.num_rows());

                let combined_row: Vec<M> = transcript
                    .read_integers(<Zip<I, L> as LinearCodes<I, L>>::row_len(vp.zip()))?;

                let encoded_combined_row: Vec<M> = vp.zip().encode_wide(&combined_row);
                encoded_combined_rows.push((coeffs, encoded_combined_row));
            }
        }

        let mut columns_opened: Vec<(usize, Vec<K>)> = Vec::with_capacity(
            <Zip<I, L> as LinearCodes<I, L>>::num_column_opening(vp.zip()),
        );
        for _ in 0..<Zip<I, L> as LinearCodes<I, L>>::num_column_opening(vp.zip()) {
            let column_idx = transcript.squeeze_challenge_idx(
                field,
                <Zip<I, L> as LinearCodes<I, L>>::codeword_len(vp.zip()),
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
        coeffs: &[I],
        encoded_combined_row: &[M],
        column_entries: &[K],
        column: usize,
        num_rows: usize,
    ) -> Result<(), Error> {
        let column_entries_comb = if num_rows > 1 {
            let coeffs: Vec<_> = coeffs.iter().map(expand::<I, M>).collect();
            let column_entries: Vec<_> = column_entries.iter().map(expand::<K, M>).collect();
            inner_product(coeffs.iter(), column_entries.iter())
        } else {
            expand::<K, M>(&column_entries[0])
        };

        if column_entries_comb != encoded_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".into()));
        }
        Ok(())
    }

    fn verify_evaluation_z<F: Field>(
        vp: &Self::VerifierParam,
        point: &[F],
        eval: F,
        columns_opened: &[(usize, Vec<K>)],
        transcript: &mut PcsTranscript<F>,
        field: F::Cr,
    ) -> Result<(), Error>
    where
        L: FieldMap<F, Output = F>,
        K: FieldMap<F, Output = F>,
    {
        let q_0_combined_row = transcript
            .read_field_elements(<Zip<I, L> as LinearCodes<I, L>>::row_len(vp.zip()), field)?;
        let encoded_combined_row = vp.zip().encode_f(&q_0_combined_row, field);

        let (q_0, q_1) = point_to_tensor(vp.num_rows(), point, field)?;

        if inner_product(&q_0_combined_row, &q_1) != eval {
            return Err(Error::InvalidPcsOpen(
                "Evaluation consistency failure".into(),
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

    fn verify_proximity_q_0<F: Field>(
        q_0: &Vec<F>,
        encoded_q_0_combined_row: &[F],
        column_entries: &[K],
        column: usize,
        num_rows: usize,
        field: F::Cr,
    ) -> Result<(), Error>
    where
        K: FieldMap<F, Output = F>,
    {
        let column_entries_comb = if num_rows > 1 {
            let column_entries = column_entries.map_to_field(field);
            inner_product(q_0, &column_entries)
            // TODO: this inner product is taking a long time.
        } else {
            column_entries.first().unwrap().map_to_field(field)
        };
        if column_entries_comb != encoded_q_0_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".into()));
        }

        Ok(())
    }
}
