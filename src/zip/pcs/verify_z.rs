use ark_std::iterable::Iterable;

use i256::{I256, I512};
use itertools::Itertools;
use sha3::{digest::Output, Digest, Keccak256};

use crate::{
    field::{conversion::FieldMap, RandomField},
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
    utils::{point_to_tensor_z, validate_input, ColumnOpening},
};

impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    // pub fn read_commitments(
    //     _: &Self::VerifierParam,
    //     num_polys: usize,
    //     transcript: &mut PcsTranscript<N>,
    // ) -> Result<Vec<Self::Commitment>, Error> {
    //     transcript.read_commitments(num_polys).map(|roots| {
    //         roots
    //             .into_iter()
    //             .map(MultilinearZipCommitment::new)
    //             .collect_vec()
    //     })
    // }

    pub fn verify_z(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &Vec<i64>,
        eval: &i64,
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        validate_input("verify", vp.num_vars(), [], [point])?;

        let row_len = vp.zip().row_len();
        let codeword_len = vp.zip().codeword_len();

        let columns_opened = Self::verify_testing(vp, comm.roots(), transcript, row_len, field)?;

        // // Retrieve the row combinations from the transcript
        // // Pair the combinations with the coefficients that generated them
        // let mut encoded_combined_rows = Vec::with_capacity(vp.zip().num_proximity_testing());

        // if vp.num_rows() > 1 {
        //     for _ in 0..vp.zip().num_proximity_testing() {
        //         let coeffs: Vec<_> = transcript
        //             .fs_transcript
        //             .get_integer_challenges(vp.num_rows())
        //             .iter()
        //             .map(|i| I256::from(*i))
        //             .collect();

        //         let combined_row = transcript.read_I256_vec(row_len)?;

        //         let code = vp.zip().encode(&combined_row);

        //         encoded_combined_rows.push((coeffs, code));
        //     }
        // }

        // let depth = codeword_len.next_power_of_two().ilog2() as usize;
        // let t_0_combined_row = transcript.read_field_elements(row_len, field)?;
        // let (t_0, t_1) = point_to_tensor_z(vp.num_rows(), point)?;
        // let t_0_f = t_0.map_to_field(field);
        // let mut columns_to_open = Vec::with_capacity(vp.zip().num_column_opening());
        // // Ensure that the test combinations are valid codewords
        // for _ in 0..vp.zip().num_column_opening() {
        //     let column = transcript.squeeze_challenge_idx(field, codeword_len);
        //     columns_to_open.push(column);

        //     let items = transcript.read_I512_vec(vp.num_rows())?;

        //     let merkle_path = transcript.read_commitments(depth)?;

        //     Self::verify_proximity_z(&encoded_combined_rows, &items, column, vp.num_rows())?;
        //     Self::verify_proximity_t_0(
        //         &t_0_f,
        //         &vp.zip().encode_f(&t_0_combined_row, field),
        //         &items,
        //         column,
        //         vp.num_rows(),
        //         field,
        //     )?;
        //     Self::verify_merkle_path(&items, &merkle_path, column, comm)?;
        // }

        // // verify consistency

        // let t_1_f = t_1
        //     .iter()
        //     .map(|i| i.map_to_field(field))
        //     .collect::<Vec<_>>();

        // let eval_f = eval.map_to_field(field);

        // if inner_product(&t_0_combined_row, &t_1_f) != eval_f {
        //     return Err(Error::InvalidPcsOpen("Consistency failure".to_string()));
        // }

        // Ok(())
        todo!()
    }

    pub fn batch_verify_z<'a>(
        vp: &Self::VerifierParam,
        comms: impl Iterable<Item = &'a MultilinearZipCommitment<N>>,
        points: &[Vec<i64>],
        evals: &[i64],
        transcript: &mut PcsTranscript<N>,
        field: *const FieldConfig<N>,
    ) -> Result<(), Error> {
        for (i, (eval, comm)) in evals.iter().zip(comms.iter()).enumerate() {
            Self::verify_z(vp, comm, &points[i], eval, transcript, field)?;
        }
        Ok(())
    }

    // pub(super) fn verify_merkle_path(
    //     items: &[I512],
    //     path: &[Output<Keccak256>],
    //     column: usize,
    //     comm: &Self::Commitment,
    // ) -> Result<(), Error> {
    //     let mut hasher = Keccak256::default();

    //     let mut output = {
    //         for item in items.iter() {
    //             <Keccak256 as sha3::digest::Update>::update(&mut hasher, &item.to_be_bytes());
    //         }

    //         hasher.clone().finalize()
    //     };

    //     hasher.reset();
    //     for (idx, neighbor) in path.iter().enumerate() {
    //         if (column >> idx) & 1 == 0 {
    //             <Keccak256 as sha3::digest::Update>::update(&mut hasher, &output);
    //             <Keccak256 as sha3::digest::Update>::update(&mut hasher, neighbor);
    //         } else {
    //             <Keccak256 as sha3::digest::Update>::update(&mut hasher, neighbor);
    //             <Keccak256 as sha3::digest::Update>::update(&mut hasher, &output);
    //         }
    //         output = hasher.clone().finalize();
    //         hasher.reset();
    //     }
    //     if output != comm.root() {
    //         return Err(Error::InvalidPcsOpen(
    //             "Invalid merkle tree opening".to_string(),
    //         ));
    //     }
    //     Ok(())
    // }

    pub(super) fn verify_testing(
        vp: &Self::VerifierParam,
        roots: &[Output<Keccak256>],
        transcript: &mut PcsTranscript<N>,
        row_len: usize,
        field: *const FieldConfig<N>,
    ) -> Result<Vec<usize>, Error> {
        // Gather the coeffs and encoded combined rows per proximity test
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

                let encoded_combined_row = vp.zip().encode(&combined_row);
                encoded_combined_rows.push((coeffs, encoded_combined_row));
            }
        }

        let mut columns_opened = vec![0; vp.zip().num_column_opening()];
        for i in 0..vp.zip().num_column_opening() {
            let column_idx = transcript.squeeze_challenge_idx(field, vp.zip().codeword_len());
            columns_opened[i] = column_idx;

            let column_values = transcript.read_I512_vec(vp.num_rows())?;

            for (coeffs, encoded_combined_row) in encoded_combined_rows.iter() {
                Self::verify_proximity_z(
                    coeffs,
                    encoded_combined_row,
                    &column_values,
                    column_idx,
                    vp.num_rows(),
                )?;
            }

            let _ =ColumnOpening::verify_column(roots, &column_values, transcript);
        }
        Ok(columns_opened)
    }

    pub(super) fn verify_proximity_z(
        coeffs: &[I256],
        encoded_combined_row: &[I512],
        column_entries: &[I512],
        column: usize,
        num_rows: usize,
    ) -> Result<(), Error> {
        let column_entries_comb = if num_rows > 1 {
            let coeff: Vec<_> = coeffs.iter().map(|i| I256_to_I512(*i)).collect();

            inner_product(&coeff, column_entries)
        } else {
            column_entries[0]
        };

        if column_entries_comb != encoded_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
        }
        Ok(())
    }

    fn verify_proximity_t_0(
        t_0_f: &Vec<RandomField<N>>,
        t_0_combined_row: &[RandomField<N>],
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
