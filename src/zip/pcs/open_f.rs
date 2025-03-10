#![allow(non_snake_case)]
use std::borrow::Cow;

use ark_ff::Zero;
use ark_std::iterable::Iterable;
use itertools::izip;
use sha3::{digest::Output, Keccak256};

use crate::{
    field::RandomField as F,
    field_config::FieldConfig,
    poly_z::mle::DenseMultilinearExtension,
    zip::{
        code::{LinearCodes, ZipSpec},
        pcs_transcript::PcsTranscript,
        utils::parallelize,
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
    pub fn open_f(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        comm: &Self::Commitment,
        point: &[F<N>],
        field: *const FieldConfig<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<Output<Keccak256>>, Error> {
        // TODO put this back as when we have a function
        // validate_input("open", pp.num_vars(), [poly], [point])?;

        let row_len = pp.zip().row_len();

        let codeword_len = pp.zip().codeword_len();

        Self::prove_proximity_f(
            pp.num_rows(),
            row_len,
            transcript,
            pp.zip().num_proximity_testing(),
            point,
            poly,
            field,
        )?;

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
    pub fn batch_open_f<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension>,
        comms: impl Iterable<Item = &'a MultilinearZipCommitment<N>>,
        points: &[Vec<F<N>>],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<Vec<Vec<Output<Keccak256>>>, Error> {
        //	use std::env;
        //	let key = "RAYON_NUM_THREADS";
        //	env::set_var(key, "8");

        let mut proofs = vec![];
        for (poly, comm, point) in izip!(polys.iter(), comms.iter(), points.iter()) {
            proofs.push(Self::open_f(pp, poly, comm, point, field, transcript)?);
        }
        Ok(proofs)
    }

    // Subprotocol functions
    fn prove_proximity_f(
        num_rows: usize,
        row_len: usize,
        transcript: &mut PcsTranscript<N>,
        num_proximity_testing: usize,
        point: &[F<N>],

        poly: &Self::Polynomial,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        let (t_0, _) = point_to_tensor_f(num_rows, point, field).unwrap();

        let evaluations: Vec<F<N>> = poly
            .evaluations
            .iter()
            .map(|i| F::from_i64(*i, field).unwrap())
            .collect();

        if num_rows > 1 {
            // If we can take linear combinations evaluation.config_ptr(
            // perform the proximity test an arbitrary number of times
            for _ in 0..num_proximity_testing {
                let coeffs = transcript.fs_transcript.get_challenges(num_rows, field);
                let combined_row = combine_rows(&coeffs, &evaluations, row_len);
                transcript.write_field_elements(&combined_row)?;
            }
        }

        let t_0_combined_row = if num_rows > 1 {
            // Return the evalauation row combination
            let combined_row = combine_rows(&t_0, &evaluations, row_len);
            Cow::<Vec<F<N>>>::Owned(combined_row)
        } else {
            // If there is only one row, we have no need to take linear combinations
            // We just return the evaluation row combination
            Cow::Borrowed(&evaluations)
        };

        transcript.write_field_elements(&t_0_combined_row)
    }
}

// Define function that performs a row operation on the evaluation matrix
// [t_0]^T * M]
fn combine_rows<const N: usize>(
    coeffs: &[F<N>],
    evaluations: &[F<N>],
    row_len: usize,
) -> Vec<F<N>> {
    let mut combined_row = Vec::with_capacity(row_len);
    parallelize(&mut combined_row, |(combined_row, offset)| {
        combined_row
            .iter_mut()
            .zip(offset..)
            .for_each(|(combined, column)| {
                *combined = F::zero();
                coeffs
                    .iter()
                    .zip(evaluations.iter().skip(column).step_by(row_len))
                    .for_each(|(coeff, eval)| {
                        *combined += &(*coeff * eval);
                    });
            })
    });

    combined_row
}
