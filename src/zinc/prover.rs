use crate::field_config::ConfigRef;
use ark_ff::Zero;
use crypto_bigint::Int;
use std::sync::atomic::Ordering;

use crate::{
    ccs::{
        ccs_f::{Statement_F, CCS_F},
        ccs_z::{Instance_Z, Statement_Z, Witness_Z, CCS_Z},
    },
    field::{conversion::FieldMap, RandomField},
    field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension,
    poly_z::mle::DenseMultilinearExtension as DenseMultilinearExtensionZ,
    sparse_matrix::SparseMatrix,
    sumcheck::{utils::build_eq_x_r, MLSumcheck, SumcheckProof},
    transcript::KeccakTranscript,
    zip::{code::ZipSpec, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};

use super::{
    errors::{MleEvaluationError, SpartanError, ZincError},
    structs::{SpartanProof, ZincProof, ZincProver, ZipProof},
    utils::{
        calculate_Mz_mles, prepare_lin_sumcheck_polynomial, sumcheck_polynomial_comb_fn_1,
        SqueezeBeta, SqueezeGamma,
    },
};

pub trait Prover<'cfg, const I: usize, const N: usize> {
    fn prove(
        &self,
        statement: &Statement_Z<I>,
        wit: &Witness_Z<I>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<I>,
        config: &'cfg FieldConfig<N>,
    ) -> Result<ZincProof<'cfg, I, N>, ZincError<N>>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:;
}

impl<'cfg, const I: usize, const N: usize, S: ZipSpec> Prover<'cfg, I, N> for ZincProver<I, N, S> {
    fn prove(
        &self,
        statement: &Statement_Z<I>,
        wit: &Witness_Z<I>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<I>,
        config: &'cfg FieldConfig<N>,
    ) -> Result<ZincProof<'cfg, I, N>, ZincError<N>>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:,
    {
        // TODO: Write functionality to let the verifier know that there are no denominators that can be divided by q(As an honest prover)
        let (z_ccs, z_mle, ccs_f, statement_f) =
            Self::prepare_for_random_field_piop(statement, wit, ccs, ConfigRef::from(config))?;

        // Prove Spartan protocol over random field
        let (spartan_proof, r_y) = SpartanProver::<I, N>::prove(
            self,
            &statement_f,
            &z_ccs,
            &z_mle,
            &ccs_f,
            transcript,
            ConfigRef::from(config),
        )?;

        // Commit to z_mle and prove its evaluation at v
        let zip_proof = Self::commit_z_mle_and_prove_evaluation(
            &z_mle,
            &ccs_f,
            &r_y,
            transcript,
            ConfigRef::from(config),
        )?;

        // Return proof
        Ok(ZincProof {
            spartan_proof,
            zip_proof,
        })
    }
}
/// Prover for the Spartan protocol
pub trait SpartanProver<'cfg, const I: usize, const N: usize> {
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
        statement_f: &Statement_F<'cfg, N>,
        z_ccs: &[RandomField<'cfg, N>],
        z_mle: &DenseMultilinearExtensionZ<I>,
        ccs_f: &CCS_F<'cfg, N>,
        transcript: &mut KeccakTranscript,
        config: ConfigRef<'cfg, N>,
    ) -> Result<(SpartanProof<'cfg, N>, Vec<RandomField<'cfg, N>>), SpartanError<N>>;
}

impl<'cfg, const I: usize, const N: usize, S: ZipSpec> SpartanProver<'cfg, I, N>
    for ZincProver<I, N, S>
{
    fn prove(
        &self,
        statement_f: &Statement_F<'cfg, N>,
        z_ccs: &[RandomField<'cfg, N>],
        z_mle: &DenseMultilinearExtensionZ<I>,
        ccs_f: &CCS_F<'cfg, N>,
        transcript: &mut KeccakTranscript,
        config: ConfigRef<'cfg, N>,
    ) -> Result<(SpartanProof<'cfg, N>, Vec<RandomField<'cfg, N>>), SpartanError<N>> {
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

impl<const I: usize, const N: usize, S: ZipSpec> ZincProver<I, N, S> {
    pub fn prepare_for_random_field_piop<'cfg>(
        statement: &Statement_Z<I>,
        wit: &Witness_Z<I>,
        ccs: &CCS_Z<I>,
        config: ConfigRef<'cfg, N>,
    ) -> Result<
        (
            Vec<RandomField<'cfg, N>>,
            DenseMultilinearExtensionZ<I>,
            CCS_F<'cfg, N>,
            Statement_F<'cfg, N>,
        ),
        SpartanError<N>,
    > {
        // z_ccs vector, i.e. concatenation x || 1 || w.
        let (z_ccs, z_mle) = Self::get_z_ccs_and_z_mle(statement, wit, ccs, config);
        let ccs_f = ccs.map_to_field(config);
        let statement_f = statement.map_to_field(config);
        Ok((z_ccs, z_mle, ccs_f, statement_f))
    }
    /// Step 2 of Fig 5: Construct polynomial $g$ and generate $\beta$ challenges.
    fn construct_polynomial_g<'cfg>(
        z_ccs: &[RandomField<'cfg, N>],
        transcript: &mut KeccakTranscript,
        constraints: &[SparseMatrix<RandomField<'cfg, N>>],
        ccs: &CCS_F<'cfg, N>,
        config: ConfigRef<'cfg, N>,
    ) -> Result<
        (
            Vec<DenseMultilinearExtension<'cfg, N>>,
            usize,
            Vec<DenseMultilinearExtension<'cfg, N>>,
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

    fn get_z_ccs_and_z_mle<'cfg>(
        statement: &Statement_Z<I>,
        wit: &Witness_Z<I>,
        ccs: &CCS_Z<I>,
        config: ConfigRef<'cfg, N>,
    ) -> (Vec<RandomField<'cfg, N>>, DenseMultilinearExtensionZ<I>) {
        let mut z_ccs = statement.get_z_vector(&wit.w_ccs);

        if z_ccs.len() <= ccs.m {
            z_ccs.resize(ccs.m, Int::<I>::ZERO);
        }
        let z_mle = DenseMultilinearExtensionZ::from_evaluations_slice(ccs.s_prime, &z_ccs);

        (
            z_ccs.into_iter().map(|x| x.map_to_field(config)).collect(),
            z_mle,
        )
    }

    fn sumcheck_1<'cfg>(
        z_ccs: &[RandomField<'cfg, N>],
        transcript: &mut KeccakTranscript,
        statement: &Statement_F<'cfg, N>,
        ccs: &CCS_F<'cfg, N>,
        config: ConfigRef<'cfg, N>,
    ) -> Result<
        (
            SumcheckProof<'cfg, N>,
            Vec<RandomField<'cfg, N>>,
            Vec<DenseMultilinearExtension<'cfg, N>>,
        ),
        SpartanError<N>,
    > {
        let (g_mles, g_degree, mz_mles) = {
            Self::construct_polynomial_g(z_ccs, transcript, &statement.constraints, ccs, config)?
        };

        let comb_fn = {
            move |vals: &[RandomField<'cfg, N>]| -> RandomField<'cfg, N> {
                sumcheck_polynomial_comb_fn_1(vals, ccs)
            }
        };

        let (sumcheck_proof_1, r_x) =
            Self::generate_sumcheck_proof(transcript, g_mles, ccs.s, g_degree, comb_fn, config)?;

        Ok((sumcheck_proof_1, r_x, mz_mles))
    }

    fn sumcheck_2<'cfg>(
        r_a: &[RandomField<'cfg, N>],
        ccs: &CCS_F<'cfg, N>,
        statement: &Statement_F<'cfg, N>,
        config: ConfigRef<'cfg, N>,
        z_mle: &DenseMultilinearExtension<'cfg, N>,
        transcript: &mut KeccakTranscript,
    ) -> Result<(SumcheckProof<'cfg, N>, Vec<RandomField<'cfg, N>>), SpartanError<N>> {
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
                ConfigRef::new(ccs.config.load(Ordering::Acquire))
            });

        sumcheck_2_mles.push(evals_mle);
        sumcheck_2_mles.push(z_mle.clone());
        let comb_fn_2 =
            |vals: &[RandomField<'cfg, N>]| -> RandomField<'cfg, N> { vals[0] * vals[1] };

        Self::generate_sumcheck_proof(transcript, sumcheck_2_mles, ccs.s, 2, comb_fn_2, config)
    }

    fn commit_z_mle_and_prove_evaluation<'cfg>(
        z_mle: &DenseMultilinearExtensionZ<I>,
        ccs: &CCS_F<'cfg, N>,
        r_y: &[RandomField<'cfg, N>],
        transcript: &mut KeccakTranscript,
        config: ConfigRef<'cfg, N>,
    ) -> Result<ZipProof<'cfg, I, N>, SpartanError<N>>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:,
    {
        let param =
            MultilinearZip::<I, { 2 * I }, { 4 * I }, { 8 * I }, S, _>::setup(ccs.m, transcript);
        let (z_data, z_comm) =
            MultilinearZip::<I, { 2 * I }, { 4 * I }, { 8 * I }, S, KeccakTranscript>::commit::<N>(
                &param, z_mle,
            )?;
        let mut pcs_transcript = PcsTranscript::new();
        let v = z_mle.map_to_field(config).evaluate(r_y, config).ok_or(
            MleEvaluationError::IncorrectLength(r_y.len(), z_mle.num_vars),
        )?;
        MultilinearZip::<I, { 2 * I }, { 4 * I }, { 8 * I }, S, KeccakTranscript>::open(
            &param,
            z_mle,
            &z_data,
            r_y,
            config,
            &mut pcs_transcript,
        )?;

        let pcs_proof = pcs_transcript.into_proof();
        Ok(ZipProof {
            z_comm,
            v,
            pcs_proof,
        })
    }

    fn calculate_V_s<'cfg>(
        mz_mles: &[DenseMultilinearExtension<'cfg, N>],
        r_x: &[RandomField<'cfg, N>],
        config: ConfigRef<'cfg, N>,
    ) -> Result<Vec<RandomField<'cfg, N>>, SpartanError<N>> {
        let V_s: Result<Vec<RandomField<'cfg, N>>, MleEvaluationError> = mz_mles
            .iter()
            .map(
                |mle: &DenseMultilinearExtension<N>| -> Result<RandomField<N>, MleEvaluationError> {
                    mle.evaluate(r_x, config)
                        .ok_or(MleEvaluationError::IncorrectLength(r_x.len(), mle.num_vars))
                },
            )
            .collect();

        let V_s = V_s?;
        Ok(V_s)
    }
    /// Step 2: Run linearization sum-check protocol.
    fn generate_sumcheck_proof<'cfg>(
        transcript: &mut KeccakTranscript,
        mles: Vec<DenseMultilinearExtension<'cfg, N>>,
        nvars: usize,
        degree: usize,
        comb_fn: impl Fn(&[RandomField<'cfg, N>]) -> RandomField<'cfg, N> + Send + Sync,
        config: ConfigRef<'cfg, N>,
    ) -> Result<(SumcheckProof<'cfg, N>, Vec<RandomField<'cfg, N>>), SpartanError<N>> {
        let (sum_check_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles, nvars, degree, &comb_fn, config);

        Ok((sum_check_proof, prover_state.randomness))
    }
}
