use ark_std::{sync::atomic::Ordering, vec::Vec};
use crypto_bigint::Int;

use super::{
    errors::{MleEvaluationError, SpartanError, ZincError},
    structs::{SpartanProof, ZincProof, ZincProver, ZipProof},
    utils::{
        calculate_Mz_mles, prepare_lin_sumcheck_polynomial, sumcheck_polynomial_comb_fn_1,
        SqueezeBeta, SqueezeGamma,
    },
};
use crate::{
    biginteger::BigInt,
    ccs::{
        ccs_f::{Statement_F, CCS_F},
        ccs_z::{Instance_Z, Statement_Z, Witness_Z, CCS_Z},
    },
    poly_f::mle::DenseMultilinearExtension,
    poly_z::mle::DenseMultilinearExtension as DenseMultilinearExtensionZ,
    sparse_matrix::SparseMatrix,
    sumcheck::{utils::build_eq_x_r, MLSumcheck, SumcheckProof},
    traits::{ConfigReference, CryptoInt, Field, FieldMap},
    transcript::KeccakTranscript,
    zip::{code::ZipSpec, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};

pub type SpartanResult<T, F> = Result<T, SpartanError<F>>;
pub type ZincResult<T, F> = Result<T, ZincError<F>>;

pub trait Prover<const I: usize, F: Field> {
    fn prove(
        &self,
        statement: &Statement_Z<Int<I>>,
        wit: &Witness_Z<Int<I>>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<Int<I>>,
        config: F::Cr,
    ) -> Result<ZincProof<Int<I>, F>, ZincError<F>>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:;
}

impl<const I: usize, F: Field, S: ZipSpec> Prover<I, F> for ZincProver<Int<I>, F, S>
where
    for<'a> Int<I>: From<&'a F::I>,
    for<'a> F::CryptoInt: From<&'a BigInt<I>>, // TODO
    F::I: From<Int<I>>,
{
    fn prove(
        &self,
        statement: &Statement_Z<Int<I>>,
        wit: &Witness_Z<Int<I>>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<Int<I>>,
        config: F::Cr,
    ) -> ZincResult<ZincProof<Int<I>, F>, F>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:,
    {
        // TODO: Write functionality to let the verifier know that there are no denominators that can be divided by q(As an honest prover)
        let (z_ccs, z_mle, ccs_f, statement_f) =
            Self::prepare_for_random_field_piop(statement, wit, ccs, config)?;

        // Prove Spartan protocol over random field
        let (spartan_proof, r_y) = SpartanProver::<Int<I>, F>::prove(
            self,
            &statement_f,
            &z_ccs,
            &z_mle,
            &ccs_f,
            transcript,
            config,
        )?;

        // Commit to z_mle and prove its evaluation at v
        let zip_proof =
            Self::commit_z_mle_and_prove_evaluation(&z_mle, &ccs_f, &r_y, transcript, config)?;

        // Return proof
        Ok(ZincProof {
            spartan_proof,
            zip_proof,
        })
    }
}
/// Prover for the Spartan protocol
pub trait SpartanProver<I: CryptoInt, F: Field> {
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
        statement_f: &Statement_F<F>,
        z_ccs: &[F],
        z_mle: &DenseMultilinearExtensionZ<I>,
        ccs_f: &CCS_F<F>,
        transcript: &mut KeccakTranscript,
        config: F::Cr,
    ) -> SpartanResult<(SpartanProof<F>, Vec<F>), F>;
}

impl<const I: usize, F: Field, S: ZipSpec> SpartanProver<Int<I>, F> for ZincProver<Int<I>, F, S>
where
    for<'a> Int<I>: From<&'a F::I>,
    for<'a> F::CryptoInt: From<&'a BigInt<I>>, // TODO
    F::I: From<Int<I>>,
{
    fn prove(
        &self,
        statement_f: &Statement_F<F>,
        z_ccs: &[F],
        z_mle: &DenseMultilinearExtensionZ<Int<I>>,
        ccs_f: &CCS_F<F>,
        transcript: &mut KeccakTranscript,
        config: F::Cr,
    ) -> SpartanResult<(SpartanProof<F>, Vec<F>), F> {
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

impl<const I: usize, F: Field, S: ZipSpec> ZincProver<Int<I>, F, S>
where
    for<'a> Int<I>: From<&'a F::I>,
    for<'a> F::CryptoInt: From<&'a BigInt<I>>, // TODO
    F::I: From<Int<I>>,
{
    pub type DenseMleF = DenseMultilinearExtension<F>;
    pub type DenseMleZ = DenseMultilinearExtensionZ<Int<I>>;
    pub type CcsF = CCS_F<F>;

    pub type StatementF = Statement_F<F>;

    #[allow(clippy::type_complexity)] // TODO refactor this out
    pub fn prepare_for_random_field_piop(
        statement: &Statement_Z<Int<I>>,
        wit: &Witness_Z<Int<I>>,
        ccs: &CCS_Z<Int<I>>,
        config: F::Cr,
    ) -> SpartanResult<(Vec<F>, Self::DenseMleZ, Self::CcsF, Self::StatementF), F> {
        // z_ccs vector, i.e. concatenation x || 1 || w.
        let (z_ccs, z_mle) = Self::get_z_ccs_and_z_mle(statement, wit, ccs, config);
        let ccs_f = ccs.map_to_field(config);
        let statement_f = statement.map_to_field(config);
        Ok((z_ccs, z_mle, ccs_f, statement_f))
    }

    /// Step 2 of Fig 5: Construct polynomial $g$ and generate $\beta$ challenges.
    #[allow(clippy::type_complexity)] // TODO refactor this out
    fn construct_polynomial_g(
        z_ccs: &[F],
        transcript: &mut KeccakTranscript,
        constraints: &[SparseMatrix<F>],
        ccs: &Self::CcsF,
        config: F::Cr,
    ) -> SpartanResult<(Vec<Self::DenseMleF>, usize, Vec<Self::DenseMleF>), F> {
        // Generate beta challenges from Step 1
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, config);

        // Prepare MLEs
        let Mz_mles = calculate_Mz_mles::<SpartanError<F>, F>(constraints, ccs.s, z_ccs, config)?;

        // Construct the sumcheck polynomial g
        let (g_mles, g_degree) =
            prepare_lin_sumcheck_polynomial(&ccs.c, &ccs.d, &Mz_mles, &ccs.S, &beta_s, config)?;

        Ok((g_mles, g_degree, Mz_mles))
    }

    fn get_z_ccs_and_z_mle(
        statement: &Statement_Z<Int<I>>,
        wit: &Witness_Z<Int<I>>,
        ccs: &CCS_Z<Int<I>>,
        config: F::Cr,
    ) -> (Vec<F>, Self::DenseMleZ) {
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

    #[allow(clippy::type_complexity)] // TODO refactor this out
    fn sumcheck_1(
        z_ccs: &[F],
        transcript: &mut KeccakTranscript,
        statement: &Self::StatementF,
        ccs: &Self::CcsF,
        config: F::Cr,
    ) -> SpartanResult<(SumcheckProof<F>, Vec<F>, Vec<Self::DenseMleF>), F> {
        let (g_mles, g_degree, mz_mles) = {
            Self::construct_polynomial_g(z_ccs, transcript, &statement.constraints, ccs, config)?
        };

        let comb_fn = { move |vals: &[F]| -> F { sumcheck_polynomial_comb_fn_1(vals, ccs) } };

        let (sumcheck_proof_1, r_x) =
            Self::generate_sumcheck_proof(transcript, g_mles, ccs.s, g_degree, comb_fn, config)?;

        Ok((sumcheck_proof_1, r_x, mz_mles))
    }

    fn sumcheck_2(
        r_a: &[F],
        ccs: &Self::CcsF,
        statement: &Self::StatementF,
        config: F::Cr,
        z_mle: &Self::DenseMleF,
        transcript: &mut KeccakTranscript,
    ) -> SpartanResult<(SumcheckProof<F>, Vec<F>), F> {
        let gamma: F = transcript.squeeze_gamma_challenge(config);
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
                        .fold(F::zero(), |mut lin_comb, eval_vec| {
                            lin_comb *= &gamma;
                            lin_comb += &eval_vec[i];
                            lin_comb
                        })
                })
                .collect::<Vec<F>>()
        };

        let evals_mle =
            DenseMultilinearExtension::from_evaluations_vec(ccs.s_prime, evals, unsafe {
                F::Cr::new(ccs.config.load(Ordering::Acquire))
            });

        sumcheck_2_mles.push(evals_mle);
        sumcheck_2_mles.push(z_mle.clone());
        let comb_fn_2 = |vals: &[F]| -> F { vals[0].clone() * &vals[1] };

        Self::generate_sumcheck_proof(transcript, sumcheck_2_mles, ccs.s, 2, comb_fn_2, config)
    }

    fn commit_z_mle_and_prove_evaluation(
        z_mle: &Self::DenseMleZ,
        ccs: &Self::CcsF,
        r_y: &[F],
        transcript: &mut KeccakTranscript,
        config: F::Cr,
    ) -> SpartanResult<ZipProof<Int<I>, F>, F>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:,
    {
        let param =
            MultilinearZip::<Int<I>, Int<{ 2 * I }>, Int<{ 4 * I }>, Int<{ 8 * I }>, S, _>::setup(
                ccs.m, transcript,
            );
        let (z_data, z_comm) = MultilinearZip::<
            Int<I>,
            Int<{ 2 * I }>,
            Int<{ 4 * I }>,
            Int<{ 8 * I }>,
            S,
            KeccakTranscript,
        >::commit::<F>(&param, z_mle)?;
        let mut pcs_transcript = PcsTranscript::new();
        let v = z_mle.map_to_field(config).evaluate(r_y, config).ok_or(
            MleEvaluationError::IncorrectLength(r_y.len(), z_mle.num_vars),
        )?;
        MultilinearZip::<
            Int<I>,
            Int<{ 2 * I }>,
            Int<{ 4 * I }>,
            Int<{ 8 * I }>,
            S,
            KeccakTranscript,
        >::open(&param, z_mle, &z_data, r_y, config, &mut pcs_transcript)?;

        let pcs_proof = pcs_transcript.into_proof();
        Ok(ZipProof {
            z_comm,
            v,
            pcs_proof,
        })
    }

    fn calculate_V_s(
        mz_mles: &[Self::DenseMleF],
        r_x: &[F],
        config: F::Cr,
    ) -> SpartanResult<Vec<F>, F> {
        let V_s: Result<Vec<F>, MleEvaluationError> = mz_mles
            .iter()
            .map(
                |mle: &DenseMultilinearExtension<F>| -> Result<F, MleEvaluationError> {
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
        mles: Vec<Self::DenseMleF>,
        nvars: usize,
        degree: usize,
        comb_fn: impl Fn(&[F]) -> F + Send + Sync,
        config: F::Cr,
    ) -> SpartanResult<(SumcheckProof<F>, Vec<F>), F> {
        let (sum_check_proof, prover_state) =
            MLSumcheck::prove_as_subprotocol(transcript, mles, nvars, degree, &comb_fn, config);

        Ok((sum_check_proof, prover_state.randomness))
    }
}
