#![allow(dead_code, non_snake_case)]
use ark_ff::Zero;
use ark_std::cfg_iter;

use errors::{MleEvaluationError, SpartanError};
use structs::{SpartanProof, ZincProver, ZincVerifier};
use utils::{
    sumcheck_polynomial_comb_fn_1, sumcheck_polynomial_comb_fn_2, SqueezeBeta, SqueezeGamma,
};

use crate::{
    brakedown::{
        code::{BrakedownSpec, BrakedownSpec1},
        pcs::MultilinearBrakedown,
    },
    ccs::{
        ccs_f::{Commitment, Instance_F, Statement, Witness, CCS_F},
        error::CSError,
        utils::mat_vec_mul,
    },
    field::RandomField,
    field_config::FieldConfig,
    poly::mle::DenseMultilinearExtension,
    sparse_matrix::SparseMatrix,
    spartan::utils::prepare_lin_sumcheck_polynomial,
    sumcheck::{utils::eq_eval, MLSumcheck, Proof, SumCheckError::SumCheckFailed},
    transcript::KeccakTranscript,
};

mod errors;
mod structs;
mod utils;

/// Prover for the Linearization subprotocol
pub trait SpartanProver<const N: usize> {
    /// Generates a proof for the linearization subprotocol
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
        statement: &Statement<N>,
        wit: &Witness<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<SpartanProof<N>, SpartanError<N>>;
}

/// Verifier for the Linearization subprotocol.
pub trait SpartanVerifier<const N: usize> {
    /// Verifies a proof for the linearization subprotocol.
    ///
    /// # Arguments
    ///
    /// * `cm_i` - A reference to a `CCCS<C, NTT>`, which represents a CCS statement and a commitment to a witness.
    /// * `proof` - A reference to a `LinearizationProof<NTT>` containing the linearization proof.
    /// * `transcript` - A mutable reference to a sponge for generating NI challenges.
    /// * `ccs` - A reference to a Customizable Constraint System instance used in the protocol.
    ///
    /// # Returns
    ///
    /// * `Ok(LCCCS<C, NTT>)` - On success, returns a linearized version of the CCS witness commitment.
    /// * `Err(LinearizationError<NTT>)` - If verification fails, returns a `LinearizationError<NTT>`.
    ///
    fn verify(
        &self,
        cm_i: &Statement<N>,
        w_commitment: &Commitment<N>,
        proof: &SpartanProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(), SpartanError<N>>;
}

impl<const N: usize, S: BrakedownSpec> SpartanProver<N> for ZincProver<N, S> {
    fn prove(
        &self,
        statement: &Statement<N>,
        wit: &Witness<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<SpartanProof<N>, SpartanError<N>> {
        // Step 1: Generate tau challenges (done in construct_polynomial_g because they are not needed
        // elsewhere.

        // Step 2: Sum check protocol.
        // z_ccs vector, i.e. concatenation x || 1 || w.
        let z_ccs = statement.get_z_vector(&wit.w_ccs);
        let w_mle = DenseMultilinearExtension::from_evaluations_slice(
            ccs.m - ccs.l - 1,
            &wit.w_ccs,
            unsafe { *ccs.config.as_ptr() },
        );
        let rng = ark_std::test_rng();
        let param = MultilinearBrakedown::<N, S>::setup(ccs.m - ccs.l - 1, ccs.m, rng);
        let w_comm = MultilinearBrakedown::<N, S>::commit(&param, &w_mle)?;
        let (g_mles, g_degree, mz_mles) = Self::construct_polynomial_g(
            &z_ccs,
            transcript,
            &statement.constraints,
            ccs,
            self.config,
        )?;

        let comb_fn = |vals: &[RandomField<N>]| -> RandomField<N> {
            sumcheck_polynomial_comb_fn_1(vals, ccs)
        };

        // Run sumcheck protocol.
        let (sumcheck_proof_1, r_a) = Self::generate_sumcheck_proof(
            transcript,
            g_mles,
            ccs.s,
            g_degree,
            comb_fn,
            self.config,
        )?;

        let (all_mles, z_mle) = calculate_poly_2_mles(
            &statement.constraints,
            &r_a,
            &z_ccs,
            ccs.s,
            ccs.s_prime,
            self.config,
        )?;
        let gamma = transcript.squeeze_gamma_challenge(self.config);
        let comb_fn_2 = |vals: &[RandomField<N>]| -> RandomField<N> {
            sumcheck_polynomial_comb_fn_2(vals, ccs, &gamma)
        };

        let V_s: Result<Vec<RandomField<N>>, MleEvaluationError> = mz_mles
            .iter()
            .map(
                |mle: &DenseMultilinearExtension<N>| -> Result<RandomField<N>, MleEvaluationError> {
                    mle.evaluate(&r_a, self.config)
                        .ok_or(MleEvaluationError::IncorrectLength(r_a.len(), mle.num_vars))
                },
            )
            .collect();

        let V_s = V_s?;

        let (sumcheck_proof_2, r_y) =
            Self::generate_sumcheck_proof(transcript, all_mles, ccs.s, 2, comb_fn_2, self.config)?;

        let v = z_mle
            .evaluate(&r_y, self.config)
            .ok_or(MleEvaluationError::IncorrectLength(
                r_y.len(),
                z_mle.num_vars,
            ));
        let v = v?;

        Ok(SpartanProof {
            linearization_sumcheck: sumcheck_proof_1,
            second_sumcheck: sumcheck_proof_2,
            V_s,
            v,
            w_comm,
        })
    }
}

impl<const N: usize, S: BrakedownSpec> SpartanVerifier<N> for ZincVerifier<N, S> {
    fn verify(
        &self,
        cm_i: &Statement<N>,
        w_commitment: &Commitment<N>,
        proof: &SpartanProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(), SpartanError<N>> {
        // Step 1: Generate the beta challenges.
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, self.config);

        //Step 2: The sumcheck.
        let (point_r, s) =
            self.verify_linearization_proof(&proof.linearization_sumcheck, transcript, ccs)?;

        // Step 3. Check V_s is congruent to s
        Self::verify_linearization_claim(&beta_s, &point_r, s, proof, ccs)?;

        let gamma = transcript.squeeze_gamma_challenge(self.config);

        let second_sumcheck_claimed_sum = Self::lin_comb_V_s(&gamma, &proof.V_s);

        let (r_y, e_y) = self.verify_second_sumcheck_proof(
            &proof.second_sumcheck,
            transcript,
            ccs,
            second_sumcheck_claimed_sum,
        )?;
        Ok(())
    }
}

impl<const N: usize, S: BrakedownSpec> ZincVerifier<N, S> {
    fn verify_linearization_proof(
        &self,
        proof: &Proof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(Vec<RandomField<N>>, RandomField<N>), SpartanError<N>> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let nvars = ccs.s;
        let degree = ccs.d + 1;

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            nvars,
            degree,
            RandomField::zero(),
            proof,
            self.config,
        )?;

        Ok((subclaim.point, subclaim.expected_evaluation))
    }

    fn verify_linearization_claim(
        beta_s: &[RandomField<N>],
        point_r: &[RandomField<N>],
        s: RandomField<N>,
        proof: &SpartanProof<N>,
        ccs: &CCS_F<N>,
    ) -> Result<(), SpartanError<N>> {
        let e = eq_eval(point_r, beta_s)?;
        let should_equal_s = e * ccs // e * (\sum c_i * \Pi_{j \in S_i} u_j)
            .c
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                c * ccs.S[i]
                    .iter()
                    .map(|&j| &proof.V_s[j])
                    .product::<RandomField<N>>()
            }) // c_i * \Pi_{j \in S_i} u_j
            .sum::<RandomField<N>>(); // \sum c_i * \Pi_{j \in S_i} u_j

        if should_equal_s != s {
            return Err(SpartanError::SumCheckError(SumCheckFailed(
                should_equal_s,
                s,
            )));
        }

        Ok(())
    }

    fn verify_second_sumcheck_proof(
        &self,
        proof: &Proof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
        claimed_sum: RandomField<N>,
    ) -> Result<(Vec<RandomField<N>>, RandomField<N>), SpartanError<N>> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let nvars = ccs.s_prime;
        let degree = 2;

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            nvars,
            degree,
            claimed_sum,
            proof,
            self.config,
        )?;

        Ok((subclaim.point, subclaim.expected_evaluation))
    }

    fn lin_comb_V_s(gamma: &RandomField<N>, V_s: &[RandomField<N>]) -> RandomField<N> {
        let mut res = RandomField::zero();
        for V_i in V_s.iter().rev() {
            res *= gamma;
            res += V_i;
        }
        res
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

    fn construct_polynomial_2(
        constraints: &[SparseMatrix<RandomField<N>>],
        r_a: &[RandomField<N>],
        z: &[RandomField<N>],
        ccs: &CCS_F<N>,
        config: *const FieldConfig<N>,
    ) -> Result<
        (
            Vec<DenseMultilinearExtension<N>>,
            DenseMultilinearExtension<N>,
        ),
        SpartanError<N>,
    > {
        let (mles, z) = calculate_poly_2_mles(constraints, r_a, z, ccs.s, ccs.s_prime, config)?;
        Ok((mles, z))
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

// Prepare MLE's of the form mle[M_i \cdot z_ccs](x), a.k.a. \sum mle[M_i](x, b) * mle[z_ccs](b).
pub fn calculate_Mz_mles<E, const N: usize>(
    constraints: &[SparseMatrix<RandomField<N>>],
    ccs_s: usize,
    z_ccs: &[RandomField<N>],
    config: *const FieldConfig<N>,
) -> Result<Vec<DenseMultilinearExtension<N>>, E>
where
    E: From<MleEvaluationError> + From<CSError> + Sync + Send,
{
    to_mles_err::<N, _, E, CSError>(
        ccs_s,
        cfg_iter!(constraints).map(|M| mat_vec_mul(M, z_ccs)),
        config,
    )
}

pub fn calculate_poly_2_mles<const N: usize>(
    constraints: &[SparseMatrix<RandomField<N>>],
    r_a: &[RandomField<N>],
    z: &[RandomField<N>],
    ccs_s: usize,
    log_n: usize,
    config: *const FieldConfig<N>,
) -> Result<
    (
        Vec<DenseMultilinearExtension<N>>,
        DenseMultilinearExtension<N>,
    ),
    MleEvaluationError,
> {
    let mut mles: Vec<DenseMultilinearExtension<N>> = constraints
        .iter()
        .map(|constraint| {
            let evaluated_vec: Vec<_> = constraint
                .to_dense()
                .iter()
                .map(|constraint_vec| {
                    DenseMultilinearExtension::<N>::from_evaluations_slice(
                        ccs_s,
                        constraint_vec,
                        config,
                    )
                    .evaluate(r_a, config)
                    .unwrap()
                })
                .collect();

            DenseMultilinearExtension::from_evaluations_vec(log_n, evaluated_vec, config)
        })
        .collect();
    let z = DenseMultilinearExtension::from_evaluations_slice(log_n, z, config);
    mles.push(z.clone());
    Ok((mles, z))
}

pub fn to_mles_err<const N: usize, I, E, E1>(
    n_vars: usize,
    mle_s: I,
    config: *const FieldConfig<N>,
) -> Result<Vec<DenseMultilinearExtension<N>>, E>
where
    I: IntoIterator<Item = Result<Vec<RandomField<N>>, E1>>,
    E: From<MleEvaluationError> + From<E1>,
{
    mle_s
        .into_iter()
        .map(|m| {
            let m = m?;
            if 1 << n_vars < m.len() {
                Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
            } else {
                Ok(DenseMultilinearExtension::from_evaluations_vec(
                    n_vars, m, config,
                ))
            }
        })
        .collect::<Result<_, E>>()
}
