use std::borrow::Cow;

use ark_ff::Zero;
use ark_std::iterable::Iterable;
use itertools::izip;
use sha3::{digest::Output, Keccak256};

use crate::{
    brakedown::{
        code::{BrakedownSpec, LinearCodes},
        pcs_transcript::PcsTranscript,
        utils::parallelize,
        Error,
    },
    field::RandomField as F,
    poly::mle::DenseMultilinearExtension,
};

use super::{
    structs::{MultilinearBrakedown, MultilinearBrakedownCommitment},
    utils::{point_to_tensor, validate_input},
};

impl<const N: usize, S> MultilinearBrakedown<N, S>
where
    S: BrakedownSpec,
{
    pub fn open(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        comm: &Self::Commitment,
        point: &Vec<F<N>>,
        eval: &F<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<Output<Keccak256>>, Error> {
        validate_input("open", pp.num_vars(), [poly], [point])?;

        let row_len = pp.brakedown().row_len();

        let codeword_len = pp.brakedown().codeword_len();

        Self::prove_proximity(
            pp.num_rows(),
            row_len,
            transcript,
            pp.brakedown().num_proximity_testing(),
            point,
            eval,
            poly,
        )?;

        let merkle_depth = codeword_len.next_power_of_two().ilog2() as usize;
        let mut proof: Vec<Output<Keccak256>> = vec![];

        Self::open_merkle_tree(
            merkle_depth,
            &mut proof,
            pp.brakedown().num_column_opening(),
            transcript,
            eval,
            codeword_len,
            comm,
        )?;

        Ok(proof)
    }

    // TODO Apply 2022/1355 https://eprint.iacr.org/2022/1355.pdf#page=30
    pub fn batch_open<'a>(
        pp: &Self::ProverParam,
        polys: impl Iterable<Item = &'a DenseMultilinearExtension<N>>,
        comms: impl Iterable<Item = &'a MultilinearBrakedownCommitment<N>>,
        points: &[Vec<F<N>>],
        evals: &[F<N>],
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<Vec<Output<Keccak256>>>, Error> {
        //	use std::env;
        //	let key = "RAYON_NUM_THREADS";
        //	env::set_var(key, "8");

        let mut proofs = vec![];
        for (eval, poly, comm, point) in
            izip!(evals.iter(), polys.iter(), comms.iter(), points.iter())
        {
            proofs.push(Self::open(pp, poly, comm, point, eval, transcript)?);
        }
        Ok(proofs)
    }
    fn prove_proximity(
        num_rows: usize,
        row_len: usize,
        transcript: &mut PcsTranscript<N>,
        num_proximity_testing: usize,
        point: &[F<N>],
        evaluation: &F<N>,
        poly: &Self::Polynomial,
    ) -> Result<(), Error> {
        let (t_0, _) = point_to_tensor(num_rows, point, evaluation.config_ptr()).unwrap();

        if num_rows > 1 {
            // If we can take linear combinations
            // perform the proximity test an arbitrary number of times
            for _ in 0..num_proximity_testing {
                let coeffs = transcript
                    .fs_transcript
                    .get_challenges(evaluation.config_ptr(), num_rows);
                let combined_row = combine_rows(&coeffs, &poly.evaluations, row_len);
                transcript.write_field_elements(&combined_row)?;
            }
        }

        let t_0_combined_row = if num_rows > 1 {
            // Return the evalauation row combination
            let combined_row = combine_rows(&t_0, &poly.evaluations, row_len);
            Cow::<Vec<F<N>>>::Owned(combined_row)
        } else {
            // If there is only one row, we have no need to take linear combinations
            // We just return the evaluation row combination
            Cow::Borrowed(&poly.evaluations)
        };

        transcript.write_field_elements(&t_0_combined_row)
    }

    fn open_merkle_tree(
        merkle_depth: usize,
        proof: &mut Vec<Output<Keccak256>>,
        num_col_opening: usize,
        transcript: &mut PcsTranscript<N>,
        eval: &F<N>,
        codeword_len: usize,
        comm: &Self::Commitment,
    ) -> Result<(), Error> {
        for _ in 0..num_col_opening {
            let column = transcript.squeeze_challenge_idx(eval.config_ptr(), codeword_len);

            transcript.write_field_elements(
                &comm
                    .rows()
                    .iter()
                    .copied()
                    .skip(column)
                    .step_by(codeword_len)
                    .collect::<Vec<_>>(),
            )?;

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
