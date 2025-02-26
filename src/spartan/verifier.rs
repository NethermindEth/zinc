use ark_ff::Zero;

use crate::{
    brakedown::{
        code::BrakedownSpec, pcs::structs::MultilinearBrakedown, pcs_transcript::PcsTranscript,
    },
    ccs::ccs_f::{Statement, CCS_F},
    field::RandomField,
    poly::mle::DenseMultilinearExtension,
    sumcheck::{utils::eq_eval, MLSumcheck, Proof, SumCheckError::SumCheckFailed},
    transcript::KeccakTranscript,
};

use super::{
    errors::{MleEvaluationError, SpartanError},
    structs::{SpartanProof, ZincVerifier},
    utils::{SqueezeBeta, SqueezeGamma},
};

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
        proof: SpartanProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(), SpartanError<N>>;
}

impl<const N: usize, S: BrakedownSpec> SpartanVerifier<N> for ZincVerifier<N, S> {
    fn verify(
        &self,
        cm_i: &Statement<N>,
        proof: SpartanProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(), SpartanError<N>> {
        let rng = ark_std::test_rng();
        let param =
            MultilinearBrakedown::<N, S>::setup(ccs.m, rng, unsafe { *ccs.config.as_ptr() });
        // Step 1: Generate the beta challenges.
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, self.config);

        //Step 2: The sumcheck.
        let (r_x, s) =
            self.verify_linearization_proof(&proof.linearization_sumcheck, transcript, ccs)?;

        // Step 3. Check V_s is congruent to s
        Self::verify_linearization_claim(&beta_s, &r_x, s, &proof, ccs)?;

        let gamma = transcript.squeeze_gamma_challenge(self.config);

        let second_sumcheck_claimed_sum = Self::lin_comb_V_s(&gamma, &proof.V_s);

        let (r_y, e_y) = self.verify_second_sumcheck_proof(
            &proof.second_sumcheck,
            transcript,
            ccs,
            second_sumcheck_claimed_sum,
        )?;

        let mut pcs_transcript = PcsTranscript::from_proof(&proof.pcs_proof);
        MultilinearBrakedown::<N, S>::verify(
            &param,
            &proof.z_comm,
            &r_y,
            &proof.v,
            &mut pcs_transcript,
        )?;

        let mut rx_ry = r_x;
        rx_ry.extend_from_slice(&r_y);

        let V_x: Result<Vec<RandomField<N>>, MleEvaluationError> = cm_i
            .constraints
            .iter()
            .map(|M| -> Result<RandomField<N>, MleEvaluationError> {
                let mle = DenseMultilinearExtension::from_matrix(M, self.config);
                mle.evaluate(&rx_ry, self.config)
                    .ok_or(MleEvaluationError::IncorrectLength(
                        rx_ry.len(),
                        mle.num_vars,
                    ))
            })
            .collect();

        let V_x = V_x?;
        let V_x_gamma = Self::lin_comb_V_s(&gamma, &V_x) * proof.v;
        if V_x_gamma != e_y {
            return Err(SpartanError::VerificationError(
                "linear combination of powers of gamma and V_x != e_y".to_string(),
            ));
        }

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
