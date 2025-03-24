use ark_ff::Zero;

use crate::{
    ccs::{
        ccs_f::{Statement_F, CCS_F},
        ccs_z::{Statement_Z, CCS_Z},
    },
    field::{conversion::FieldMap, RandomField},
    poly_f::mle::DenseMultilinearExtension,
    sumcheck::{utils::eq_eval, MLSumcheck, SumCheckError::SumCheckFailed, SumcheckProof},
    transcript::KeccakTranscript,
    zip::{code::ZipSpec, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};

use super::{
    errors::{LookupError, MleEvaluationError, SpartanError, ZincError},
    structs::{LookupProof, SpartanProof, ZincProof, ZincVerifier},
    utils::{draw_random_field, SqueezeBeta, SqueezeGamma},
};

pub trait Verifier<const N: usize> {
    fn verify(
        &self,
        cm_i: &Statement_Z,
        proof: ZincProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z,
    ) -> Result<(), ZincError<N>>;
}

impl<const N: usize, S: ZipSpec> Verifier<N> for ZincVerifier<N, S> {
    fn verify(
        &self,
        cm_i: &Statement_Z,
        proof: ZincProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z,
    ) -> Result<(), ZincError<N>> {
        let field_config = draw_random_field::<N>(&cm_i.public_input, transcript);
        // TODO: Write functionality to let the verifier know that there are no denominators that can be divided by q(As an honest prover)
        let ccs_F = ccs.map_to_field(field_config);
        let statement_f = cm_i.map_to_field(field_config);

        let (rx_ry, e_y, gamma) =
            SpartanVerifier::<N>::verify(self, &proof.spartan_proof, transcript, &ccs_F)
                .map_err(ZincError::SpartanError)?;

        self.verify_pcs_proof(&statement_f, &proof, &rx_ry, e_y, gamma, &ccs_F, transcript)?;

        LookupVerifier::<N>::verify(self, proof.lookup_proof).map_err(ZincError::LookupError)?;

        Ok(())
    }
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
        proof: &SpartanProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(Vec<RandomField<N>>, RandomField<N>, RandomField<N>), SpartanError<N>>;
}

impl<const N: usize, S: ZipSpec> SpartanVerifier<N> for ZincVerifier<N, S> {
    fn verify(
        &self,
        proof: &SpartanProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(Vec<RandomField<N>>, RandomField<N>, RandomField<N>), SpartanError<N>> {
        // Step 1: Generate the beta challenges.
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, unsafe { *ccs.config.as_ptr() });

        //Step 2: The sumcheck.
        let (r_x, s) =
            self.verify_linearization_proof(&proof.linearization_sumcheck, transcript, ccs)?;

        // Step 3. Check V_s is congruent to s
        Self::verify_linearization_claim(&beta_s, &r_x, s, proof, ccs)?;

        let gamma = transcript.squeeze_gamma_challenge(unsafe { *ccs.config.as_ptr() });

        let second_sumcheck_claimed_sum = Self::lin_comb_V_s(&gamma, &proof.V_s);

        let (r_y, e_y) = self.verify_second_sumcheck_proof(
            &proof.second_sumcheck,
            transcript,
            ccs,
            second_sumcheck_claimed_sum,
        )?;

        Ok(([r_x, r_y].concat(), e_y, gamma))
    }
}

pub trait LookupVerifier<const N: usize> {
    fn verify(&self, proof: LookupProof<N>) -> Result<(), LookupError>;
}

impl<const N: usize, S: ZipSpec> LookupVerifier<N> for ZincVerifier<N, S> {
    fn verify(&self, _proof: LookupProof<N>) -> Result<(), LookupError> {
        todo!()
    }
}

impl<const N: usize, S: ZipSpec> ZincVerifier<N, S> {
    fn verify_linearization_proof(
        &self,
        proof: &SumcheckProof<N>,
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
            unsafe { *ccs.config.as_ptr() },
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
        proof: &SumcheckProof<N>,
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
            unsafe { *ccs.config.as_ptr() },
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

    fn verify_pcs_proof(
        &self,
        cm_i: &Statement_F<N>,
        proof: &ZincProof<N>,
        rx_ry: &[RandomField<N>],
        e_y: RandomField<N>,
        gamma: RandomField<N>,
        ccs: &CCS_F<N>,
        transcript: &mut KeccakTranscript,
    ) -> Result<(), SpartanError<N>> {
        let param = MultilinearZip::<N, S, _>::setup(ccs.m, transcript);
        let mut pcs_transcript = PcsTranscript::from_proof(&proof.pcs_proof);
        let r_y = &rx_ry[ccs.s..];

        MultilinearZip::<N, S, KeccakTranscript>::verify_f(
            &param,
            &proof.z_comm,
            r_y,
            &proof.v,
            &mut pcs_transcript,
            unsafe { *ccs.config.as_ptr() },
        )?;

        // Evaluate constraints at rx_ry point
        let V_x = cm_i.constraints.iter().map(|M| {
            let mle = DenseMultilinearExtension::from_matrix(M, unsafe { *ccs.config.as_ptr() });
            mle.evaluate(rx_ry, unsafe { *ccs.config.as_ptr() })
                .ok_or(MleEvaluationError::IncorrectLength(rx_ry.len(), mle.num_vars))
        }).collect::<Result<Vec<_>, _>>()?;

        // Check final verification equation
        let V_x_gamma = Self::lin_comb_V_s(&gamma, &V_x) * proof.v;
        if V_x_gamma != e_y {
            return Err(SpartanError::VerificationError(
                "linear combination of powers of gamma and V_x != e_y".to_string(),
            ));
        }

        Ok(())
    }
}
