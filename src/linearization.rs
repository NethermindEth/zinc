#![allow(dead_code, non_snake_case)]
use ark_std::cfg_iter;

use errors::{LinearizationError, MleEvaluationError};
use structs::{LinearizationProof, ZincLinearizationProver};
use utils::{
    sumcheck_polynomial_comb_fn_1, sumcheck_polynomial_comb_fn_2, SqueezeBeta, SqueezeGamma,
};

use crate::{
    ccs::{
        ccs_f::{Instance_F, LStatement, Statement, Witness, CCS_F},
        error::CSError,
        utils::mat_vec_mul,
    },
    field::RandomField,
    field_config::FieldConfig,
    linearization::utils::prepare_lin_sumcheck_polynomial,
    poly::mle::DenseMultilinearExtension,
    sparse_matrix::SparseMatrix,
    sumcheck::{MLSumcheck, Proof},
    transcript::KeccakTranscript,
};

mod errors;
mod structs;
mod utils;
/// Prover for the Linearization subprotocol
pub trait LinearizationProver<const N: usize> {
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
        statement: &Statement<N>,
        wit: &Witness<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
        config: *const FieldConfig<N>,
    ) -> Result<(Proof<N>, Proof<N>), LinearizationError<N>>;
}

/// Verifier for the Linearization subprotocol.
pub trait LinearizationVerifier<const N: usize> {
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
        cm_i: &Statement<N>,
        proof: &LinearizationProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<LStatement<N>, LinearizationError<N>>;
}

impl<const N: usize> LinearizationProver<N> for ZincLinearizationProver<N> {
    fn prove(
        statement: &Statement<N>,
        wit: &Witness<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
        config: *const FieldConfig<N>,
    ) -> Result<(Proof<N>, Proof<N>), LinearizationError<N>> {
        // Step 1: Generate beta challenges (done in construct_polynomial_g because they are not needed
        // elsewhere.

        // Step 2: Sum check protocol.
        // z_ccs vector, i.e. concatenation x || 1 || w.
        let z_ccs = statement.get_z_vector(&wit.w_ccs);
        let (g_mles, g_degree, _) =
            Self::construct_polynomial_g(&z_ccs, transcript, &statement.constraints, ccs, config)?;

        let comb_fn = |vals: &[RandomField<N>]| -> RandomField<N> {
            sumcheck_polynomial_comb_fn_1(vals, ccs)
        };

        // Run sumcheck protocol.
        let (sumcheck_proof_1, r_a) =
            Self::generate_sumcheck_proof(transcript, g_mles, ccs.s, g_degree, comb_fn, config)?;

        let mles = calculate_poly_2_mles(
            &statement.constraints,
            &r_a,
            &z_ccs,
            ccs.s,
            ccs.s_prime,
            config,
        )?;
        let gamma = transcript.squeeze_gamma_challenge(config);
        let comb_fn_2 = |vals: &[RandomField<N>]| -> RandomField<N> {
            sumcheck_polynomial_comb_fn_2(vals, ccs, &gamma)
        };

        let (sumcheck_proof_2, _) =
            Self::generate_sumcheck_proof(transcript, mles, ccs.s, 2, comb_fn_2, config)?;

        Ok((sumcheck_proof_1, sumcheck_proof_2))
    }
}

impl<const N: usize> ZincLinearizationProver<N> {
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
        LinearizationError<N>,
    > {
        // Generate beta challenges from Step 1
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, config);

        // Prepare MLEs
        let Mz_mles =
            calculate_Mz_mles::<LinearizationError<N>, N>(constraints, ccs.s, z_ccs, config)?;

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
    ) -> Result<Vec<DenseMultilinearExtension<N>>, LinearizationError<N>> {
        let mles = calculate_poly_2_mles(constraints, r_a, z, ccs.s, ccs.s_prime, config)?;
        Ok(mles)
    }

    /// Step 2: Run linearization sum-check protocol.
    fn generate_sumcheck_proof(
        transcript: &mut KeccakTranscript,
        mles: Vec<DenseMultilinearExtension<N>>,
        nvars: usize,
        degree: usize,
        comb_fn: impl Fn(&[RandomField<N>]) -> RandomField<N>,
        config: *const FieldConfig<N>,
    ) -> Result<(Proof<N>, Vec<RandomField<N>>), LinearizationError<N>> {
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
) -> Result<Vec<DenseMultilinearExtension<N>>, MleEvaluationError> {
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

    mles.push(DenseMultilinearExtension::from_evaluations_slice(
        log_n, z, config,
    ));
    Ok(mles)
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
