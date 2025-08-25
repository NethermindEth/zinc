use ark_std::{format, iterable::Iterable, vec::Vec};

use super::{
    structs::{MultilinearZip, MultilinearZipCommitment},
    utils::{ColumnOpening, point_to_tensor, validate_input},
};
use crate::{
    traits::{Field, FieldMap, ZipTypes},
    zip::{
        Error,
        code::LinearCode,
        pcs::{structs::MultilinearZipParams, utils::MtHash},
        pcs_transcript::PcsTranscript,
        utils::{expand, inner_product},
    },
};

impl<ZT: ZipTypes, LC: LinearCode<ZT>> MultilinearZip<ZT, LC> {
    pub fn verify<F: Field>(
        vp: &MultilinearZipParams<ZT, LC>,
        comm: &MultilinearZipCommitment,
        point: &[F],
        eval: F,
        transcript: &mut PcsTranscript<F>,
        field: F::R,
    ) -> Result<(), Error>
    where
        ZT::L: FieldMap<F, Output = F>,
        ZT::K: FieldMap<F, Output = F>,
    {
        validate_input::<ZT::N, F>("verify", vp.num_vars, [], [point])?;

        let columns_opened = Self::verify_testing(vp, &comm.root, transcript, field)?;

        Self::verify_evaluation_z(vp, point, eval, &columns_opened, transcript, field)?;

        Ok(())
    }

    pub fn batch_verify_z<'a, F: Field>(
        vp: &MultilinearZipParams<ZT, LC>,
        comms: impl Iterable<Item = &'a MultilinearZipCommitment>,
        points: &[Vec<F>],
        evals: &[F],
        transcript: &mut PcsTranscript<F>,
        field: F::R,
    ) -> Result<(), Error>
    where
        ZT::L: FieldMap<F, Output = F>,
        ZT::K: FieldMap<F, Output = F>,
        ZT::N: 'a,
    {
        for (i, (eval, comm)) in evals.iter().zip(comms.iter()).enumerate() {
            Self::verify(vp, comm, &points[i], eval.clone(), transcript, field)?;
        }
        Ok(())
    }

    #[allow(clippy::type_complexity)]
    pub(super) fn verify_testing<F: Field>(
        vp: &MultilinearZipParams<ZT, LC>,
        root: &MtHash,
        transcript: &mut PcsTranscript<F>,
        field: F::R,
    ) -> Result<Vec<(usize, Vec<ZT::K>)>, Error> {
        // Gather the coeffs and encoded combined rows per proximity test
        let mut encoded_combined_rows: Vec<(Vec<ZT::N>, Vec<ZT::M>)> =
            Vec::with_capacity(vp.linear_code.num_proximity_testing());

        if vp.num_rows > 1 {
            for _ in 0..vp.linear_code.num_proximity_testing() {
                let coeffs = transcript.fs_transcript.get_integer_challenges(vp.num_rows);

                let combined_row: Vec<ZT::M> =
                    transcript.read_integers(vp.linear_code.row_len())?;

                let encoded_combined_row: Vec<ZT::M> = vp.linear_code.encode_wide(&combined_row);
                encoded_combined_rows.push((coeffs, encoded_combined_row));
            }
        }

        let mut columns_opened: Vec<(usize, Vec<ZT::K>)> =
            Vec::with_capacity(vp.linear_code.num_column_opening());

        for _ in 0..vp.linear_code.num_column_opening() {
            let column_idx = transcript.squeeze_challenge_idx(field, vp.linear_code.codeword_len());
            let column_values = transcript.read_integers(vp.num_rows)?;

            for (coeffs, encoded_combined_row) in encoded_combined_rows.iter() {
                Self::verify_column_testing(
                    coeffs,
                    encoded_combined_row,
                    &column_values,
                    column_idx,
                    vp.num_rows,
                )?;
            }

            ColumnOpening::verify_column(root, &column_values, column_idx, transcript).map_err(
                |e| Error::InvalidPcsOpen(format!("Column opening verification failed: {e}")),
            )?;
            // TODO: Verify column opening is taking a long time.
            columns_opened.push((column_idx, column_values));
        }

        Ok(columns_opened)
    }

    pub(super) fn verify_column_testing(
        coeffs: &[ZT::N],
        encoded_combined_row: &[ZT::M],
        column_entries: &[ZT::K],
        column: usize,
        num_rows: usize,
    ) -> Result<(), Error> {
        let column_entries_comb: ZT::M = if num_rows > 1 {
            let coeffs: Vec<ZT::M> = coeffs.iter().map(expand::<ZT::N, ZT::M>).collect();
            let column_entries: Vec<ZT::M> =
                column_entries.iter().map(expand::<ZT::K, ZT::M>).collect();
            inner_product(coeffs.iter(), column_entries.iter())
        } else {
            expand(&column_entries[0])
        };

        if column_entries_comb != encoded_combined_row[column] {
            return Err(Error::InvalidPcsOpen("Proximity failure".into()));
        }
        Ok(())
    }

    fn verify_evaluation_z<F: Field>(
        vp: &MultilinearZipParams<ZT, LC>,
        point: &[F],
        eval: F,
        columns_opened: &[(usize, Vec<ZT::K>)],
        transcript: &mut PcsTranscript<F>,
        field: F::R,
    ) -> Result<(), Error>
    where
        ZT::L: FieldMap<F, Output = F>,
        ZT::K: FieldMap<F, Output = F>,
    {
        let q_0_combined_row = transcript.read_field_elements(vp.linear_code.row_len(), field)?;
        let encoded_combined_row = vp.linear_code.encode_f(&q_0_combined_row, field);

        let (q_0, q_1) = point_to_tensor(vp.num_rows, point, field)?;

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
                vp.num_rows,
                field,
            )?;
        }

        Ok(())
    }

    fn verify_proximity_q_0<F: Field>(
        q_0: &Vec<F>,
        encoded_q_0_combined_row: &[F],
        column_entries: &[ZT::K],
        column: usize,
        num_rows: usize,
        field: F::R,
    ) -> Result<(), Error>
    where
        ZT::K: FieldMap<F, Output = F>,
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

#[cfg(test)]
mod tests {
    use ark_std::{boxed::Box, vec, vec::Vec};

    use super::*;
    use crate::{
        field::{ConfigRef, Int, RandomField},
        field_config,
        poly_z::mle::DenseMultilinearExtension,
        traits::{ConfigReference, Integer, Words, ZipTypes},
        zip::{
            code::DefaultLinearCodeSpec,
            code_raa::RaaCode,
            pcs,
            pcs::{
                structs::{MultilinearZip, MultilinearZipParams},
                tests::MockTranscript,
            },
        },
    };

    const INT_LIMBS: usize = 1;
    const FIELD_LIMBS: usize = 4;

    type ZT = pcs::tests::RandomFieldZipTypes<1>;
    type F<'cfg> = RandomField<'cfg, FIELD_LIMBS>;
    type LC = RaaCode<ZT>;
    type TestZip = MultilinearZip<ZT, LC>;

    #[allow(clippy::type_complexity)]
    fn setup_full_protocol(
        num_vars: usize,
    ) -> (
        MultilinearZipParams<ZT, LC>,
        MultilinearZipCommitment,
        Vec<F<'static>>,
        F<'static>,
        Vec<u8>,
        ConfigRef<'static, FIELD_LIMBS>,
    ) {
        let poly_size = 1 << num_vars;
        let evaluations: Vec<_> = (0..poly_size as i32).map(Int::<INT_LIMBS>::from).collect();
        let poly = DenseMultilinearExtension::from_evaluations_vec(num_vars, evaluations);

        let mut keccak = MockTranscript::default();
        let linear_code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak);
        let pp = TestZip::setup(poly_size, linear_code);

        let config: &'static crate::field::FieldConfig<FIELD_LIMBS> = Box::leak(Box::new(
            field_config!(57316695564490278656402085503, FIELD_LIMBS),
        ));

        let config_ref = ConfigRef::from(config);

        let (data, comm) = TestZip::commit::<F>(&pp, &poly).unwrap();

        let point_int: Vec<Int<INT_LIMBS>> =
            (0..num_vars).map(|i| Int::from(i as i32 + 2)).collect();
        let point_f: Vec<F> = point_int.map_to_field(config_ref);

        let mut prover_transcript = PcsTranscript::new();
        TestZip::open(
            &pp,
            &poly,
            &data,
            &point_f,
            config_ref,
            &mut prover_transcript,
        )
        .unwrap();
        let proof = prover_transcript.into_proof();

        let eval: F = match poly.evaluate(&point_int) {
            None => panic!("failed to evaluate polynomial"),
            Some(p) => p,
        }
        .map_to_field(config_ref);

        (pp, comm, point_f, eval, proof, config_ref)
    }

    #[test]
    fn successful_verification_of_valid_proof() {
        let num_vars = 4;
        let (pp, comm, point_f, eval, proof, config) = setup_full_protocol(num_vars);

        let mut verifier_transcript = PcsTranscript::from_proof(&proof);
        let result = TestZip::verify(&pp, &comm, &point_f, eval, &mut verifier_transcript, config);

        assert!(result.is_ok());
    }

    #[test]
    fn verification_fails_with_incorrect_evaluation() {
        let num_vars = 4;
        let (pp, comm, point_f, eval, proof, config) = setup_full_protocol(num_vars);

        let one: F = 1i32.map_to_field(config);

        let incorrect_eval = eval + one;
        let mut verifier_transcript = PcsTranscript::from_proof(&proof);
        let result = TestZip::verify(
            &pp,
            &comm,
            &point_f,
            incorrect_eval,
            &mut verifier_transcript,
            config,
        );

        assert!(result.is_err());
    }

    #[test]
    fn verification_fails_with_tampered_proof() {
        let num_vars = 4;
        let (pp, comm, point_f, eval, proof, config) = setup_full_protocol(num_vars);

        let mut tampered = proof.clone();
        tampered[0] ^= 0x01;

        let mut verifier_transcript = PcsTranscript::from_proof(&tampered);
        let result = TestZip::verify(&pp, &comm, &point_f, eval, &mut verifier_transcript, config);

        assert!(result.is_err());
    }

    #[test]
    fn verification_fails_with_wrong_commitment() {
        let num_vars = 4;
        let (pp, _comm_poly1, point_f, eval, proof_poly1, config) = setup_full_protocol(num_vars);

        let different_evals: Vec<_> = (20..(20 + (1 << num_vars))).map(Int::from).collect();
        let poly2 = DenseMultilinearExtension::from_evaluations_vec(num_vars, different_evals);
        let (_, comm_poly2) = TestZip::commit::<F>(&pp, &poly2).unwrap();

        let mut transcript = PcsTranscript::from_proof(&proof_poly1);
        let result = TestZip::verify(&pp, &comm_poly2, &point_f, eval, &mut transcript, config);

        assert!(result.is_err());
    }

    #[test]
    fn verification_fails_with_invalid_point_size() {
        let num_vars = 4;
        let (pp, comm, _point_f, eval, proof, config) = setup_full_protocol(num_vars);
        let mut invalid_point = vec![];
        for i in 0..=num_vars {
            invalid_point.push((100 + i as i32).map_to_field(config));
        }

        let mut transcript = PcsTranscript::from_proof(&proof);
        let result = TestZip::verify(&pp, &comm, &invalid_point, eval, &mut transcript, config);

        assert!(matches!(result, Err(..)));
    }

    #[test]
    fn verification_fails_if_proximity_check_is_invalid() {
        let mut keccak = MockTranscript::default();
        let poly_size = 8; // row_len=4, num_rows=2 -> proximity checks are active
        let linear_code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak);
        let pp = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = (0..poly_size as i32).map(Int::<INT_LIMBS>::from).collect();
        let n = 3;
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

        let (data, comm) = TestZip::commit::<F>(&pp, &mle).expect("commit should succeed");

        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config_ref = ConfigRef::from(&config);
        let point_int = [0i64, 0i64, 0i64]
            .into_iter()
            .map(Int::from)
            .collect::<Vec<_>>();
        let point: Vec<F> = point_int.map_to_field(config_ref);
        let eval = mle.evaluate(&point_int).unwrap().map_to_field(config_ref);

        let mut prover_tr = PcsTranscript::<F>::new();
        TestZip::open(&pp, &mle, &data, &point, config_ref, &mut prover_tr)
            .expect("open should succeed");
        let mut proof = prover_tr.into_proof();

        let row_len = pp.linear_code.row_len();
        let bytes_per_int = <<ZT as ZipTypes>::M as Integer>::W::num_words() * 8;
        let first_combined_row_bytes = row_len * bytes_per_int;
        assert!(
            first_combined_row_bytes <= proof.len(),
            "proof too small to tamper"
        );

        let flip_at = bytes_per_int * (row_len / 2);
        proof[flip_at] ^= 0x01;

        let mut ver_tr = PcsTranscript::<F>::from_proof(&proof);
        config_ref.reference().expect("Field config cannot be none");
        let res = TestZip::verify(&pp, &comm, &point, eval, &mut ver_tr, config_ref);

        match res {
            Err(crate::zip::Error::InvalidPcsOpen(msg)) => {
                assert_eq!(msg, "Proximity failure");
            }
            Ok(()) => panic!("verification unexpectedly succeeded"),
            Err(e) => panic!("unexpected error: {e:?}"),
        }
    }

    #[test]
    fn verification_fails_if_evaluation_consistency_check_is_invalid() {
        let mut keccak = MockTranscript::default();
        let poly_size = 8;
        let linear_code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak);
        let pp = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = (0..poly_size as i32).map(Int::<INT_LIMBS>::from).collect();
        let n = 3;
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

        let (data, comm) = TestZip::commit::<F>(&pp, &mle).expect("commit should succeed");

        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config_ref = ConfigRef::from(&config);
        let point_int = [0i64, 0i64, 0i64]
            .into_iter()
            .map(Int::from)
            .collect::<Vec<_>>();
        let point: Vec<F> = point_int.map_to_field(config_ref);
        let eval = mle.evaluate(&point_int).unwrap().map_to_field(config_ref);

        let mut prover_tr = PcsTranscript::<F>::new();
        TestZip::open(&pp, &mle, &data, &point, config_ref, &mut prover_tr)
            .expect("open should succeed");
        let mut proof = prover_tr.into_proof();

        let row_len = pp.linear_code.row_len();
        let bytes_per_field = <F<'static> as Field>::W::num_words() * 8;
        let q0_bytes = row_len * bytes_per_field;
        assert!(
            proof.len() >= q0_bytes,
            "proof too small to contain q_0_combined_row"
        );

        let tail_start = proof.len() - q0_bytes;
        let flip_at = tail_start + (bytes_per_field / 2);
        proof[flip_at] ^= 0x01;

        let mut ver_tr = PcsTranscript::<F>::from_proof(&proof);
        let res = TestZip::verify(&pp, &comm, &point, eval, &mut ver_tr, config_ref);

        match res {
            Err(crate::zip::Error::InvalidPcsOpen(msg)) => {
                assert_eq!(msg, "Evaluation consistency failure");
            }
            Ok(()) => panic!("verification unexpectedly succeeded"),
            Err(e) => panic!("unexpected error: {e:?}"),
        }
    }

    #[test]
    fn verification_succeeds_for_zero_polynomial() {
        let mut keccak = MockTranscript::default();
        let poly_size = 8;
        let linear_code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak);
        let pp = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = vec![Int::<INT_LIMBS>::from(0); poly_size];
        let n = 3;
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

        let (data, comm) = TestZip::commit::<F>(&pp, &mle).expect("commit should succeed");

        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config_ref = ConfigRef::from(&config);
        let point_int = [0i64, 0i64, 0i64]
            .into_iter()
            .map(Int::from)
            .collect::<Vec<_>>();
        let point: Vec<F> = point_int.map_to_field(config_ref);
        let eval = mle.evaluate(&point_int).unwrap().map_to_field(config_ref);
        let mut prover_tr = PcsTranscript::<F>::new();
        TestZip::open(&pp, &mle, &data, &point, config_ref, &mut prover_tr)
            .expect("open should succeed");
        let proof = prover_tr.into_proof();

        let mut ver_tr = PcsTranscript::<F>::from_proof(&proof);
        let res = TestZip::verify(&pp, &comm, &point, eval, &mut ver_tr, config_ref);
        assert!(res.is_ok());
    }

    #[test]
    fn verification_succeeds_at_zero_point() {
        let mut keccak = MockTranscript::default();
        let poly_size = 8;
        let linear_code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak);
        let pp = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = (1..=poly_size as i32).map(Int::<INT_LIMBS>::from).collect();
        let n = 3;
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

        let (data, comm) = TestZip::commit::<F>(&pp, &mle).expect("commit should succeed");

        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config_ref = ConfigRef::from(&config);
        let point_int = vec![Int::from(0i64); n];
        let point: Vec<F> = point_int.map_to_field(config_ref);

        let eval = mle.evaluate(&point_int).unwrap().map_to_field(config_ref);

        let mut prover_tr = PcsTranscript::<F>::new();
        TestZip::open(&pp, &mle, &data, &point, config_ref, &mut prover_tr)
            .expect("open should succeed");
        let proof = prover_tr.into_proof();

        let mut ver_tr = PcsTranscript::<F>::from_proof(&proof);
        let res = TestZip::verify(&pp, &comm, &point, eval, &mut ver_tr, config_ref);
        assert!(res.is_ok());
    }

    #[test]
    fn verification_fails_if_proximity_values_are_too_large() {
        let mut keccak = MockTranscript::default();
        let poly_size = 8;
        let linear_code = LC::new(&DefaultLinearCodeSpec, poly_size, &mut keccak);
        let pp = TestZip::setup(poly_size, linear_code);

        let evaluations: Vec<_> = (1..=poly_size as i32).map(Int::<INT_LIMBS>::from).collect();
        let n = 3;
        let mle = DenseMultilinearExtension::from_evaluations_slice(n, &evaluations);

        let (data, comm) = TestZip::commit::<F>(&pp, &mle).expect("commit should succeed");

        let config = field_config!(57316695564490278656402085503, FIELD_LIMBS);
        let config_ref = ConfigRef::from(&config);
        let point_int = [0i64, 0i64, 0i64]
            .into_iter()
            .map(Int::from)
            .collect::<Vec<_>>();
        let point: Vec<F> = point_int.map_to_field(config_ref);
        let eval = mle.evaluate(&point_int).unwrap().map_to_field(config_ref);

        let mut prover_tr = PcsTranscript::<F>::new();
        TestZip::open(&pp, &mle, &data, &point, config_ref, &mut prover_tr)
            .expect("open should succeed");
        let mut proof = prover_tr.into_proof();

        let row_len = pp.linear_code.row_len();
        let bytes_per_int = <<ZT as ZipTypes>::M as Integer>::W::num_words() * 8;
        let first_section_bytes = row_len * bytes_per_int;
        assert!(
            first_section_bytes <= proof.len(),
            "proof too small to tamper u'"
        );

        for b in &mut proof[0..bytes_per_int] {
            *b = 0xFF;
        }

        let mut ver_tr = PcsTranscript::<F>::from_proof(&proof);
        let res = TestZip::verify(&pp, &comm, &point, eval, &mut ver_tr, config_ref);
        assert!(res.is_err());
    }
}
