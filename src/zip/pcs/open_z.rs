#![allow(non_snake_case)]
use std::borrow::Cow;

use ark_std::iterable::Iterable;
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
    utils::{point_to_tensor_z, validate_input, ColumnOpening},
};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub fn open_z(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        comm: &Self::Data,
        point: &Vec<i64>,
        field: *const FieldConfig<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<ColumnOpening>, Error> {
        validate_input("open", pp.num_vars(), [poly], [point])?;

        let proof = Self::prove_test(pp, poly, comm, transcript, field)?;

        let row_len = pp.zip().row_len();
        Self::prove_evaluation_z(pp.num_rows(), row_len, transcript, point, poly, field)?;

        Ok(proof)
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30
    pub fn batch_open<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension>,
        comms: impl Iterable<Item = &'a MultilinearZipData<N>>,
        points: &[Vec<i64>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<Vec<Vec<ColumnOpening>>, Error> {
        let mut proofs = vec![];
        for (poly, comm, point) in izip!(polys.iter(), comms.iter(), points.iter()) {
            proofs.push(Self::open_z(pp, poly, comm, point, field, transcript)?);
        }
        Ok(proofs)
    }

    // Subprotocol functions
    fn prove_evaluation_z(
        num_rows: usize,
        row_len: usize,
        transcript: &mut PcsTranscript<N>,
        point: &[i64],
        poly: &Self::Polynomial,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        // We prove evaluations over the field,so integers need to be mapped to field elements first
        let (t_0, _) = point_to_tensor_z(num_rows, point).unwrap();

        let t_O_f = t_0.map_to_field(field);

        let evaluations = poly.evaluations.map_to_field(field);

        let t_0_combined_row = if num_rows > 1 {
            // Return the evaluation row combination
            let combined_row = combine_rows(t_O_f, evaluations, row_len);
            Cow::<Vec<RandomField<N>>>::Owned(combined_row)
        } else {
            // If there is only one row, we have no need to take linear combinations
            // We just return the evaluation row combination
            Cow::Borrowed(&evaluations)
        };

        transcript.write_field_elements(&t_0_combined_row)
    }
}
