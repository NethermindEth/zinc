use ark_ff::Zero;

use crate::{
    brakedown::{
        code::BrakedownSpec,
        pcs::structs::{MultilinearBrakedown, MultilinearBrakedownCommitment},
        pcs_transcript::PcsTranscript,
    },
    ccs::{
        ccs_f::{Instance_F, Statement_F, Witness_F, CCS_F},
        ccs_z::{Statement_Z, Witness_Z, CCS_Z},
    },
    field::{conversion::FieldMap, RandomField},
    field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension,
    sparse_matrix::SparseMatrix,
    sumcheck::{utils::build_eq_x_r, MLSumcheck, Proof},
    transcript::KeccakTranscript,
};

use super::{
    errors::{MleEvaluationError, SpartanError, ZincError},
    structs::{LookupProof, SpartanProof, ZincProof, ZincProver},
    utils::{
        calculate_Mz_mles, draw_random_field, prepare_lin_sumcheck_polynomial,
        sumcheck_polynomial_comb_fn_1, SqueezeBeta, SqueezeGamma,
    },
};

pub trait Prover<const N: usize> {
    fn prove(
        &self,
        statement: &Statement_Z,
        wit: &Witness_Z,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z,
    ) -> Result<ZincProof<N>, ZincError<N>>;
}

impl<const N: usize, S: BrakedownSpec> Prover<N> for ZincProver<N, S> {
    fn prove(
        &self,
        statement: &Statement_Z,
        wit: &Witness_Z,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z,
    ) -> Result<ZincProof<N>, ZincError<N>> {
        let field_config = draw_random_field::<N>(&statement.public_input, transcript);
        // TODO: Write functionality to let the verifier know that there are no denominators that can be divided by q(As an honest prover)
        let ccs_F = ccs.map_to_field(field_config);
        let wit_F = wit.map_to_field(field_config);
        let statement_F = statement.map_to_field(field_config);
        let spartan_proof = SpartanProver::<N>::prove(
            self,
            &statement_F,
            &wit_F,
            transcript,
            &ccs_F,
            field_config,
        )?;
        let lookup_proof = LookupProver::<N>::prove(self, wit)?;
        Ok(ZincProof {
            spartan_proof,
            lookup_proof,
        })
    }
}
/// Prover for the Spartan protocol
pub trait SpartanProver<const N: usize> {
    /// Generates a proof for the spartan protocol
    ///
    /// # Arguments
    ///
    /// * `cm_i` - A reference to a committed CCS statement to be linearized, i.e. a CCCS<C, NTT>.
    /// * `wit` - A reference to a CCS witness for the statement cm_i.
    /// * `transcript` - A mutable reference to a sponge for generating NI challenges.
    /// * `ccs` - A reference to a Customizable Constraint System circuit representation.
    ///
    /// # Returns
    ///
    /// On success, returns a tuple `(LCCCS<C, NTT>, LinearizationProof<NTT>)` where:
    ///   * `LCCCS<C, NTT>` is a linearized version of the CCS witness commitment.
    ///   * `LinearizationProof<NTT>` is a proof that the linearization subprotocol was executed correctly.
    ///
    /// # Errors
    ///
    /// Returns an error if asked to evaluate MLEs with incorrect number of variables
    ///
    fn prove(
        &self,
        statement: &Statement_F<N>,
        wit: &Witness_F<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
        config: *const FieldConfig<N>,
    ) -> Result<SpartanProof<N>, SpartanError<N>>;
}

impl<const N: usize, S: BrakedownSpec> SpartanProver<N> for ZincProver<N, S> {
    fn prove(
        &self,
        statement: &Statement_F<N>,
        wit: &Witness_F<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
        config: *const FieldConfig<N>,
    ) -> Result<SpartanProof<N>, SpartanError<N>> {
        // z_ccs vector, i.e. concatenation x || 1 || w.
        let (z_ccs, z_mle) = Self::get_z_ccs_and_z_mle(&statement, &wit, config, &ccs);

        // Do first Sumcheck
        let (sumcheck_proof_1, r_a, mz_mles) =
            Self::sumcheck_1(&z_ccs, transcript, &statement, &ccs, config)?;

        // Do second sumcheck
        let (sumcheck_proof_2, r_y) =
            Self::sumcheck_2(&r_a, &ccs, &statement, config, &z_mle, transcript)?;

        // Commit to z_mle and prove its evaluation at v
        let (z_comm, v, pcs_proof) =
            Self::commit_z_mle_and_prove_evaluation(&z_mle, &ccs, config, &r_y)?;

        // Calculate V_s
        let V_s = Self::calculate_V_s(&mz_mles, &r_a, config)?;

        // TODO: Add lookup argument for enforcing integers

        // Return proof
        Ok(SpartanProof {
            linearization_sumcheck: sumcheck_proof_1,
            second_sumcheck: sumcheck_proof_2,
            V_s,
            v,
            z_comm,
            pcs_proof,
        })
    }
}

pub trait LookupProver<const N: usize> {
    fn prove(&self, wit: &Witness_Z) -> Result<LookupProof<N>, SpartanError<N>>;
}

impl<const N: usize, S: BrakedownSpec> LookupProver<N> for ZincProver<N, S> {
    fn prove(&self, _wit: &Witness_Z) -> Result<LookupProof<N>, SpartanError<N>> {
        todo!()
    }
}

impl<const N: usize, S: BrakedownSpec> ZincProver<N, S> {
    /// Step 2 of Fig 5: Construct polynomial $g$ and generate $\beta$ challenges.
    fn construct_polynomial_g(
        z_ccs: &[RandomField<N>],
        transcript: &mut KeccakTranscript,
        constraints: &[SparseMatrix<RandomField<N>>],
        ccs: &CCS_F<N>,
        config: *const FieldConfig<N>,
    ) -> Result<
        (
            Vec<DenseMultilinearExtension<N>>,
            usize,
            Vec<DenseMultilinearExtension<N>>,
        ),
        SpartanError<N>,
    > {
        // Generate beta challenges from Step 1
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, config);

        // Prepare MLEs
        let Mz_mles = calculate_Mz_mles::<SpartanError<N>, N>(constraints, ccs.s, z_ccs, config)?;

        // Construct the sumcheck polynomial g
        let (g_mles, g_degree) =
            prepare_lin_sumcheck_polynomial(&ccs.c, &ccs.d, &Mz_mles, &ccs.S, &beta_s, config)?;

        Ok((g_mles, g_degree, Mz_mles))
    }

    fn get_z_ccs_and_z_mle(
        statement: &Statement_F<N>,
        wit: &Witness_F<N>,
        config: *const FieldConfig<N>,
        ccs: &CCS_F<N>,
    ) -> (Vec<RandomField<N>>, DenseMultilinearExtension<N>) {
        let mut z_ccs = statement.get_z_vector(&wit.w_ccs, config);

        if z_ccs.len() <= ccs.m {
            z_ccs.resize(
                ccs.m,
                RandomField::new_unchecked(unsafe { *ccs.config.as_ptr() }, 0u32.into()),
            )
        }
        let z_mle = DenseMultilinearExtension::from_evaluations_slice(ccs.s_prime, &z_ccs, config);

        (z_ccs, z_mle)
    }

    fn sumcheck_1(
        z_ccs: &[RandomField<N>],
        transcript: &mut KeccakTranscript,
        statement: &Statement_F<N>,
        ccs: &CCS_F<N>,
        config: *const FieldConfig<N>,
    ) -> Result<
        (
            Proof<N>,
            Vec<RandomField<N>>,
            Vec<DenseMultilinearExtension<N>>,
        ),
        SpartanError<N>,
    > {
        let (g_mles, g_degree, mz_mles) =
            Self::construct_polynomial_g(z_ccs, transcript, &statement.constraints, ccs, config)?;

        let comb_fn = |vals: &[RandomField<N>]| -> RandomField<N> {
            sumcheck_polynomial_comb_fn_1(vals, ccs)
        };

        // Run sumcheck protocol.
        let (sumcheck_proof_1, r_a) =
            Self::generate_sumcheck_proof(transcript, g_mles, ccs.s, g_degree, comb_fn, config)?;
        Ok((sumcheck_proof_1, r_a, mz_mles))
    }

    fn sumcheck_2(
        r_a: &[RandomField<N>],
        ccs: &CCS_F<N>,
        statement: &Statement_F<N>,
        config: *const FieldConfig<N>,
        z_mle: &DenseMultilinearExtension<N>,
        transcript: &mut KeccakTranscript,
    ) -> Result<(Proof<N>, Vec<RandomField<N>>), SpartanError<N>> {
        let gamma = transcript.squeeze_gamma_challenge(config);
        let mut sumcheck_2_mles = Vec::with_capacity(2);

        let eq_r_a = build_eq_x_r(r_a, config)?;
        let evals = {
            // compute the initial evaluation table for R(r_a, x)

            let evals_vec =
                statement.compute_eval_table_sparse(ccs.n, ccs.m, ccs, &eq_r_a.evaluations);

            (0..evals_vec[0].len())
                .map(|i| {
                    evals_vec
                        .iter()
                        .rev()
                        .fold(RandomField::zero(), |mut lin_comb, eval_vec| {
                            lin_comb *= gamma;
                            lin_comb += &eval_vec[i];
                            lin_comb
                        })
                })
                .collect::<Vec<RandomField<N>>>()
        };

        let evals_mle =
            DenseMultilinearExtension::from_evaluations_vec(ccs.s_prime, evals, unsafe {
                *(ccs.config.as_ptr())
            });

        sumcheck_2_mles.push(evals_mle);
        sumcheck_2_mles.push(z_mle.clone());
        let comb_fn_2 = |vals: &[RandomField<N>]| -> RandomField<N> { vals[0] * vals[1] };

        Self::generate_sumcheck_proof(transcript, sumcheck_2_mles, ccs.s, 2, comb_fn_2, config)
    }

    fn commit_z_mle_and_prove_evaluation(
        z_mle: &DenseMultilinearExtension<N>,
        ccs: &CCS_F<N>,
        config: *const FieldConfig<N>,
        r_y: &Vec<RandomField<N>>,
    ) -> Result<(MultilinearBrakedownCommitment<N>, RandomField<N>, Vec<u8>), SpartanError<N>> {
        let rng = ark_std::test_rng();
        let param = MultilinearBrakedown::<N, S>::setup(ccs.m, rng, config);
        let z_comm = MultilinearBrakedown::<N, S>::commit(&param, z_mle)?;
        let mut pcs_transcript = PcsTranscript::new();
        let v = z_mle
            .evaluate(r_y, config)
            .ok_or(MleEvaluationError::IncorrectLength(
                r_y.len(),
                z_mle.num_vars,
            ))?;
        MultilinearBrakedown::<N, S>::open(&param, z_mle, &z_comm, r_y, &v, &mut pcs_transcript)?;

        let pcs_proof = pcs_transcript.into_proof();
        Ok((z_comm, v, pcs_proof))
    }

    fn calculate_V_s(
        mz_mles: &[DenseMultilinearExtension<N>],
        r_a: &[RandomField<N>],
        config: *const FieldConfig<N>,
    ) -> Result<Vec<RandomField<N>>, SpartanError<N>> {
        let V_s: Result<Vec<RandomField<N>>, MleEvaluationError> = mz_mles
            .iter()
            .map(
                |mle: &DenseMultilinearExtension<N>| -> Result<RandomField<N>, MleEvaluationError> {
                    mle.evaluate(r_a, config)
                        .ok_or(MleEvaluationError::IncorrectLength(r_a.len(), mle.num_vars))
                },
            )
            .collect();

        let V_s = V_s?;
        Ok(V_s)
    }
    /// Step 2: Run linearization sum-check protocol.
    fn generate_sumcheck_proof(
        transcript: &mut KeccakTranscript,
        mles: Vec<DenseMultilinearExtension<N>>,
        nvars: usize,
        degree: usize,
        comb_fn: impl Fn(&[RandomField<N>]) -> RandomField<N> + Send + Sync,
        config: *const FieldConfig<N>,
    ) -> Result<(Proof<N>, Vec<RandomField<N>>), SpartanError<N>> {
        let (sum_check_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles, nvars, degree, comb_fn, config);

        Ok((sum_check_proof, prover_state.randomness))
    }
}
