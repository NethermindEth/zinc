use ark_std::{sync::atomic::Ordering, vec::Vec};
use num_traits::Zero;

use super::{
    errors::{MleEvaluationError, SpartanError, ZincError},
    structs::{SpartanProof, ZincProof, ZincProver, ZipProof},
    utils::{
        SqueezeBeta, SqueezeGamma, calculate_Mz_mles, prepare_lin_sumcheck_polynomial,
        sumcheck_polynomial_comb_fn_1,
    },
};
use crate::{
    ccs::{
        ccs_f::{CCS_F, Statement_F},
        ccs_z::{CCS_Z, Instance_Z, Statement_Z, Witness_Z},
    },
    field::RandomField,
    poly_f::mle::DenseMultilinearExtension,
    poly_z::mle::DenseMultilinearExtension as DenseMultilinearExtensionZ,
    sparse_matrix::SparseMatrix,
    sumcheck::{MLSumcheck, SumcheckProof, utils::build_eq_x_r},
    traits::{ConfigReference, FieldMap, FromRef, Integer, MapsToField, ZipTypes},
    transcript::KeccakTranscript,
    zip::{
        code::LinearCodeSpec, code_raa::RaaCode, pcs::structs::MultilinearZip,
        pcs_transcript::PcsTranscript,
    },
};

pub type SpartanResult<T> = Result<T, SpartanError>;
pub type ZincResult<T> = Result<T, ZincError>;

pub trait Prover<I: Integer, C: ConfigReference> {
    fn prove(
        &self,
        statement: &Statement_Z<I>,
        wit: &Witness_Z<I>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<I>,
        config: C,
    ) -> Result<ZincProof<C>, ZincError>;
}

impl<ZT: ZipTypes, C: ConfigReference, S: LinearCodeSpec> Prover<ZT::N, C> for ZincProver<ZT, C, S>
where
    ZT::N: FromRef<C::I>,
    C::I: FromRef<<ZT::N as Integer>::I>, // TODO
    C::I: FromRef<ZT::N>,
    ZT::N: MapsToField<C>,
{
    fn prove(
        &self,
        statement: &Statement_Z<ZT::N>,
        wit: &Witness_Z<ZT::N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<ZT::N>,
        config: C,
    ) -> ZincResult<ZincProof<C>> {
        // TODO: Write functionality to let the verifier know that there are no denominators that can be divided by q(As an honest prover)
        let (z_ccs, z_mle, ccs_f, statement_f) =
            Self::prepare_for_random_field_piop(statement, wit, ccs, config)?;

        // Prove Spartan protocol over random field
        let (spartan_proof, r_y) = SpartanProver::<ZT::N, C>::prove(
            self,
            &statement_f,
            &z_ccs,
            &z_mle,
            &ccs_f,
            transcript,
            config,
        )?;

        // Commit to z_mle and prove its evaluation at v
        let zip_proof = Self::commit_z_mle_and_prove_evaluation(
            &self.lc_spec,
            &z_mle,
            &ccs_f,
            &r_y,
            transcript,
            config,
        )?;

        // Return proof
        Ok(ZincProof {
            spartan_proof,
            zip_proof,
        })
    }
}

/// Prover for the Spartan protocol
pub trait SpartanProver<I: Integer, C: ConfigReference> {
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
        statement_f: &Statement_F<C>,
        z_ccs: &[RandomField<C>],
        z_mle: &DenseMultilinearExtensionZ<I>,
        ccs_f: &CCS_F<C>,
        transcript: &mut KeccakTranscript,
        config: C,
    ) -> SpartanResult<(SpartanProof<C>, Vec<RandomField<C>>)>;
}

impl<ZT: ZipTypes, C: ConfigReference, S: LinearCodeSpec> SpartanProver<ZT::N, C>
    for ZincProver<ZT, C, S>
where
    ZT::N: FromRef<C::I>,
    C::I: FromRef<<ZT::N as Integer>::I>, // TODO
    C::I: FromRef<ZT::N>,
    ZT::N: MapsToField<C>,
{
    fn prove(
        &self,
        statement_f: &Statement_F<C>,
        z_ccs: &[RandomField<C>],
        z_mle: &DenseMultilinearExtensionZ<ZT::N>,
        ccs_f: &CCS_F<C>,
        transcript: &mut KeccakTranscript,
        config: C,
    ) -> SpartanResult<(SpartanProof<C>, Vec<RandomField<C>>)> {
        // Do first Sumcheck
        let (sumcheck_proof_1, r_x, mz_mles) =
            Self::sumcheck_1(z_ccs, transcript, statement_f, ccs_f, config)?;

        // Do second sumcheck
        let (sumcheck_proof_2, r_y) = Self::sumcheck_2(
            &r_x,
            ccs_f,
            statement_f,
            config,
            &z_mle.map_to_field(config),
            transcript,
        )?;

        let V_s = Self::calculate_V_s(&mz_mles, &r_x, config)?;

        let proof = SpartanProof {
            linearization_sumcheck: sumcheck_proof_1,
            second_sumcheck: sumcheck_proof_2,
            V_s,
        };
        Ok((proof, r_y))
    }
}

impl<ZT: ZipTypes, C: ConfigReference, S: LinearCodeSpec> ZincProver<ZT, C, S>
where
    ZT::N: FromRef<C::I>,
    C::I: FromRef<<ZT::N as Integer>::I>, // TODO
    C::I: FromRef<ZT::N>,
    ZT::N: MapsToField<C>,
{
    #[allow(clippy::type_complexity)] // TODO refactor this out
    pub fn prepare_for_random_field_piop(
        statement: &Statement_Z<ZT::N>,
        wit: &Witness_Z<ZT::N>,
        ccs: &CCS_Z<ZT::N>,
        config: C,
    ) -> SpartanResult<(
        Vec<RandomField<C>>,
        DenseMultilinearExtensionZ<ZT::N>,
        CCS_F<C>,
        Statement_F<C>,
    )> {
        // z_ccs vector, i.e. concatenation x || 1 || w.
        let (z_ccs, z_mle) = Self::get_z_ccs_and_z_mle(statement, wit, ccs, config);
        let ccs_f = ccs.map_to_field(config);
        let statement_f = statement.map_to_field(config);
        Ok((z_ccs, z_mle, ccs_f, statement_f))
    }

    /// Step 2 of Fig 5: Construct polynomial $g$ and generate $\beta$ challenges.
    #[allow(clippy::type_complexity)] // TODO refactor this out
    fn construct_polynomial_g(
        z_ccs: &[RandomField<C>],
        transcript: &mut KeccakTranscript,
        constraints: &[SparseMatrix<RandomField<C>>],
        ccs: &CCS_F<C>,
        config: C,
    ) -> SpartanResult<(
        Vec<DenseMultilinearExtension<C>>,
        usize,
        Vec<DenseMultilinearExtension<C>>,
    )> {
        // Generate beta challenges from Step 1
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, config);

        // Prepare MLEs
        let Mz_mles = calculate_Mz_mles::<SpartanError, C>(constraints, ccs.s, z_ccs, config)?;

        // Construct the sumcheck polynomial g
        let (g_mles, g_degree) =
            prepare_lin_sumcheck_polynomial(&ccs.c, &ccs.d, &Mz_mles, &ccs.S, &beta_s, config)?;

        Ok((g_mles, g_degree, Mz_mles))
    }

    fn get_z_ccs_and_z_mle(
        statement: &Statement_Z<ZT::N>,
        wit: &Witness_Z<ZT::N>,
        ccs: &CCS_Z<ZT::N>,
        config: C,
    ) -> (Vec<RandomField<C>>, DenseMultilinearExtensionZ<ZT::N>) {
        let mut z_ccs = statement.get_z_vector(&wit.w_ccs);

        if z_ccs.len() <= ccs.m {
            z_ccs.resize(ccs.m, ZT::N::zero());
        }
        let z_mle = DenseMultilinearExtensionZ::from_evaluations_slice(ccs.s_prime, &z_ccs);

        (
            z_ccs.into_iter().map(|x| x.map_to_field(config)).collect(),
            z_mle,
        )
    }

    #[allow(clippy::type_complexity)] // TODO refactor this out
    fn sumcheck_1(
        z_ccs: &[RandomField<C>],
        transcript: &mut KeccakTranscript,
        statement: &Statement_F<C>,
        ccs: &CCS_F<C>,
        config: C,
    ) -> SpartanResult<(
        SumcheckProof<C>,
        Vec<RandomField<C>>,
        Vec<DenseMultilinearExtension<C>>,
    )> {
        let (g_mles, g_degree, mz_mles) = {
            Self::construct_polynomial_g(z_ccs, transcript, &statement.constraints, ccs, config)?
        };

        let comb_fn = {
            move |vals: &[RandomField<C>]| -> RandomField<C> {
                sumcheck_polynomial_comb_fn_1(vals, ccs)
            }
        };

        let (sumcheck_proof_1, r_x) =
            Self::generate_sumcheck_proof(transcript, g_mles, ccs.s, g_degree, comb_fn, config)?;

        Ok((sumcheck_proof_1, r_x, mz_mles))
    }

    fn sumcheck_2(
        r_a: &[RandomField<C>],
        ccs: &CCS_F<C>,
        statement: &Statement_F<C>,
        config: C,
        z_mle: &DenseMultilinearExtension<C>,
        transcript: &mut KeccakTranscript,
    ) -> SpartanResult<(SumcheckProof<C>, Vec<RandomField<C>>)> {
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
                            lin_comb *= &gamma;
                            lin_comb += &eval_vec[i];
                            lin_comb
                        })
                })
                .collect::<Vec<_>>()
        };

        let evals_mle =
            DenseMultilinearExtension::from_evaluations_vec(ccs.s_prime, evals, unsafe {
                C::new(ccs.config.load(Ordering::Acquire))
            });

        sumcheck_2_mles.push(evals_mle);
        sumcheck_2_mles.push(z_mle.clone());
        let comb_fn_2 = |vals: &[RandomField<C>]| -> RandomField<C> { vals[0].clone() * &vals[1] };

        Self::generate_sumcheck_proof(transcript, sumcheck_2_mles, ccs.s, 2, comb_fn_2, config)
    }

    fn commit_z_mle_and_prove_evaluation(
        lc_spec: &S,
        z_mle: &DenseMultilinearExtensionZ<ZT::N>,
        ccs: &CCS_F<C>,
        r_y: &[RandomField<C>],
        transcript: &mut KeccakTranscript,
        config: C,
    ) -> SpartanResult<ZipProof<C>> {
        let linear_code = RaaCode::<ZT>::new(lc_spec, ccs.m, transcript);
        let param = MultilinearZip::<ZT, _>::setup(ccs.m, linear_code);
        let (z_data, z_comm) = MultilinearZip::<ZT, _>::commit::<C>(&param, z_mle)?;
        let mut pcs_transcript = PcsTranscript::new();
        let v = z_mle.map_to_field(config).evaluate(r_y, config).ok_or(
            MleEvaluationError::IncorrectLength(r_y.len(), z_mle.num_vars),
        )?;
        MultilinearZip::<ZT, _>::open(&param, z_mle, &z_data, r_y, config, &mut pcs_transcript)?;

        let pcs_proof = pcs_transcript.into_proof();
        Ok(ZipProof {
            z_comm,
            v,
            pcs_proof,
        })
    }

    fn calculate_V_s(
        mz_mles: &[DenseMultilinearExtension<C>],
        r_x: &[RandomField<C>],
        config: C,
    ) -> SpartanResult<Vec<RandomField<C>>> {
        let V_s: Result<Vec<RandomField<C>>, MleEvaluationError> = mz_mles
            .iter()
            .map(
                |mle: &DenseMultilinearExtension<C>| -> Result<RandomField<C>, MleEvaluationError> {
                    mle.evaluate(r_x, config)
                        .ok_or(MleEvaluationError::IncorrectLength(r_x.len(), mle.num_vars))
                },
            )
            .collect();

        let V_s = V_s?;
        Ok(V_s)
    }
    /// Step 2: Run linearization sum-check protocol.
    fn generate_sumcheck_proof(
        transcript: &mut KeccakTranscript,
        mles: Vec<DenseMultilinearExtension<C>>,
        nvars: usize,
        degree: usize,
        comb_fn: impl Fn(&[RandomField<C>]) -> RandomField<C> + Send + Sync,
        config: C,
    ) -> SpartanResult<(SumcheckProof<C>, Vec<RandomField<C>>)> {
        let (sum_check_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles, nvars, degree, &comb_fn, config);

        Ok((sum_check_proof, prover_state.randomness))
    }
}
