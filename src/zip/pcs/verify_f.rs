use ark_std::iterable::Iterable;

use i256::{I256, I512};

use crate::{
    field::{conversion::FieldMap, RandomField as F},
    field_config::FieldConfig,
    zip::{
        code::{I256_to_I512, ZipSpec},
        pcs_transcript::PcsTranscript,
        utils::inner_product,
        Error,
    },
};

use super::structs::{MultilinearZip, MultilinearZipCommitment, ZipTranscript};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
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
                .map(|i| i.map_to_field(field))
                .collect::<Vec<_>>();
            inner_product(t_0_f, &column_entries)
        } else {
            column_entries[0].map_to_field(field)
        };
        if column_entries_comb != t_0_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
        }

        Ok(())
    }
}
