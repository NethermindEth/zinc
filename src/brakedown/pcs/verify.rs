use ark_ff::Zero;
use ark_std::iterable::Iterable;
use itertools::Itertools;
use sha3::{digest::Output, Digest, Keccak256};

use crate::{
    brakedown::{
        code::{BrakedownSpec, LinearCodes},
        pcs_transcript::PcsTranscript,
        utils::inner_product,
        Error,
    },
    field::RandomField as F,
};

use super::{
    structs::{MultilinearBrakedown, MultilinearBrakedownCommitment},
    utils::{point_to_tensor, validate_input},
};

impl<const N: usize, S> MultilinearBrakedown<N, S>
where
    S: BrakedownSpec,
{
    pub fn read_commitments(
        _: &Self::VerifierParam,
        num_polys: usize,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<Self::Commitment>, Error> {
        transcript.read_commitments(num_polys).map(|roots| {
            roots
                .into_iter()
                .map(MultilinearBrakedownCommitment::from_root)
                .collect_vec()
        })
    }

    pub fn verify(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &Vec<F<N>>,
        eval: &F<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        validate_input("verify", vp.num_vars(), [], [point])?;

        let row_len = vp.brakedown().row_len();
        let codeword_len = vp.brakedown().codeword_len();

        // Retrieve the row combinations from the transcript
        // Pair the combinations with the coefficients that generated them
        let mut combined_rows = Vec::with_capacity(vp.brakedown().num_proximity_testing() + 1);

        if vp.num_rows() > 1 {
            for _ in 0..vp.brakedown().num_proximity_testing() {
                let coeffs = transcript
                    .fs_transcript
                    .get_challenges(eval.config_ptr(), vp.num_rows());
                let mut combined_row =
                    transcript.read_field_elements(row_len, eval.config_ptr())?;

                combined_row.resize(codeword_len, F::zero());
                vp.brakedown().encode(&mut combined_row);
                combined_rows.push((coeffs, combined_row));
            }
        }

        let (t_0, t_1) = point_to_tensor(vp.num_rows(), point, eval.config_ptr())?;
        combined_rows.push({
            let mut t_0_combined_row =
                transcript.read_field_elements(row_len, eval.config_ptr())?;

            t_0_combined_row.resize(codeword_len, F::zero());
            vp.brakedown().encode(&mut t_0_combined_row);
            (t_0, t_0_combined_row)
        });

        let depth = codeword_len.next_power_of_two().ilog2() as usize;

        // Ensure that the test combinations are valid codewords
        for _ in 0..vp.brakedown().num_column_opening() {
            let column = transcript.squeeze_challenge_idx(eval.config_ptr(), codeword_len);

            let items = transcript.read_field_elements(vp.num_rows(), eval.config_ptr())?;
            let merkle_path = transcript.read_commitments(depth)?;

            Self::verify_proximity(&combined_rows, &items, column, vp.num_rows())?;
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

    pub fn batch_verify<'a>(
        vp: &Self::VerifierParam,
        comms: impl Iterable<Item = &'a MultilinearBrakedownCommitment<N>>,
        points: &[Vec<F<N>>],
        evals: &[F<N>],
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        for (i, (eval, comm)) in evals.iter().zip(comms.iter()).enumerate() {
            Self::verify(vp, comm, &points[i], eval, transcript)?;
        }
        Ok(())
    }

    fn verify_merkle_path(
        items: &[F<N>],
        path: &[Output<Keccak256>],
        column: usize,
        comm: &Self::Commitment,
    ) -> Result<(), Error> {
        let mut hasher = Keccak256::default();
        let mut output = {
            for item in items.iter() {
                <Keccak256 as sha3::digest::Update>::update(
                    &mut hasher,
                    &item.value().to_bytes_be(),
                );
            }

            hasher.clone().finalize()
        };
        for (idx, neighbor) in path.iter().enumerate() {
            if (column >> idx) & 1 == 0 {
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, &output);
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, neighbor);
            } else {
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, neighbor);
                <Keccak256 as sha3::digest::Update>::update(&mut hasher, &output);
            }
            output = hasher.clone().finalize();
        }

        if &output != comm.root() {
            return Err(Error::InvalidPcsOpen(
                "Invalid merkle tree opening".to_string(),
            ));
        }
        Ok(())
    }

    fn verify_proximity(
        combined_rows: &[(Vec<F<N>>, Vec<F<N>>)],
        items: &[F<N>],
        column: usize,
        num_rows: usize,
    ) -> Result<(), Error> {
        for (coeff, encoded) in combined_rows.iter() {
            let item = if num_rows > 1 {
                inner_product(coeff, items)
            } else {
                items[0]
            };
            if item != encoded[column] {
                return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
            }
        }
        Ok(())
    }
}
