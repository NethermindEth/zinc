#![allow(non_snake_case)]
use std::borrow::Cow;

use ark_std::iterable::Iterable;
use itertools::izip;
use sha3::{digest::Output, Keccak256};

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
    structs::{MultilinearZip, MultilinearZipCommitment, ZipTranscript},
    utils::{point_to_tensor_z, validate_input},
};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub fn open_z(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        comm: &Self::Commitment,
        point: &Vec<i64>,
        field: *const FieldConfig<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<Output<Keccak256>>, Error> {
        validate_input("open", pp.num_vars(), [poly], [point])?;

        let row_len = pp.zip().row_len();

        let codeword_len = pp.zip().codeword_len();
        Self::prove_test(
            pp.num_rows(),
            row_len,
            pp.zip().num_proximity_testing(),
            poly,
            transcript,
        )?;
        Self::prove_evaluation_z(pp.num_rows(), row_len, transcript, point, poly, field)?;

        let merkle_depth = codeword_len.next_power_of_two().ilog2() as usize;
        let mut proof: Vec<Output<Keccak256>> = vec![];

        Self::open_merkle_tree(
            merkle_depth,
            &mut proof,
            pp.zip().num_column_opening(),
            transcript,
            codeword_len,
            comm,
            field,
        )?;

        Ok(proof)
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30
    pub fn batch_open<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension>,
        comms: impl Iterable<Item = &'a MultilinearZipCommitment<N>>,
        points: &[Vec<i64>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<Vec<Vec<Output<Keccak256>>>, Error> {
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

        let t_O_f: Vec<RandomField<N>> = t_0.iter().map(|i| i.map_to_field(field)).collect();

        let evaluations: Vec<RandomField<N>> = poly
            .evaluations
            .iter()
            .map(|i| i.map_to_field(field))
            .collect();

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

    pub(crate) fn open_merkle_tree(
        merkle_depth: usize,
        proof: &mut Vec<Output<Keccak256>>,
        num_col_opening: usize,
        transcript: &mut PcsTranscript<N>,

        codeword_len: usize,
        comm: &Self::Commitment,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        for _ in 0..num_col_opening {
            let column = transcript.squeeze_challenge_idx(field, codeword_len);

            //Write the elements in the squeezed column to the shared transcript
            transcript.write_I512_vec(
                &comm
                    .rows()
                    .iter()
                    .copied()
                    .skip(column)
                    .step_by(codeword_len)
                    .collect::<Vec<_>>(),
            )?;

            //Write the neighbour hash path to the shared transcript
            let mut offset = 0;
            for (idx, width) in (1..=merkle_depth).rev().map(|depth| 1 << depth).enumerate() {
                let neighbor_idx = (column >> idx) ^ 1;
                transcript.write_commitment(&comm.intermediate_hashes()[offset + neighbor_idx])?;

                proof.push(comm.intermediate_hashes()[offset + neighbor_idx]);
                offset += width;
            }
        }
        Ok(())
    }
}
