use ark_ff::Zero;

use super::{
    errors::{MleEvaluationError, SpartanError, ZincError},
    structs::{SpartanProof, ZincProof, ZincVerifier, ZipProof},
    utils::{draw_random_field, SqueezeBeta, SqueezeGamma},
};
use crate::field_config::ConfigRef;
use crate::{
    ccs::{
        ccs_f::{Statement_F, CCS_F},
        ccs_z::{Statement_Z, CCS_Z},
    },
    field::{conversion::FieldMap, RandomField},
    field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension,
    sumcheck::{utils::eq_eval, MLSumcheck, SumCheckError::SumCheckFailed, SumcheckProof},
    transcript::KeccakTranscript,
    zip::{code::ZipSpec, pcs::structs::MultilinearZip, pcs_transcript::PcsTranscript},
};

pub trait Verifier<const I: usize, const N: usize> {
    fn verify<'cfg>(
        &self,
        cm_i: &Statement_Z<I>,
        proof: ZincProof<'cfg, I, N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<I>,
        config: &'cfg FieldConfig<N>,
    ) -> Result<(), ZincError<N>>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:;
}

impl<const I: usize, const N: usize, S: ZipSpec> Verifier<I, N> for ZincVerifier<I, N, S> {
    fn verify<'cfg>(
        &self,
        statement: &Statement_Z<I>,
        proof: ZincProof<'cfg, I, N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<I>,
        config: &'cfg FieldConfig<N>,
    ) -> Result<(), ZincError<N>>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:,
    {
        if draw_random_field::<I, N>(&statement.public_input, transcript) != *config {
            return Err(ZincError::FieldConfigError);
        }
        // TODO: Write functionality to let the verifier know that there are no denominators that can be divided by q(As an honest prover)
        let ccs_F = ccs.map_to_field(ConfigRef::from(config));
        let statement_f = statement.map_to_field(ConfigRef::from(config));

        let verification_points =
            SpartanVerifier::<N>::verify(self, &proof.spartan_proof, &ccs_F, transcript, config)
                .map_err(ZincError::SpartanError)?;

        self.verify_pcs_proof(
            &statement_f,
            &proof.zip_proof,
            &verification_points,
            &ccs_F,
            transcript,
            config,
        )?;

        Ok(())
    }
}

/// Verifier for the Linearization subprotocol.
pub trait SpartanVerifier<'cfg, const N: usize> {
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
        proof: &SpartanProof<'cfg, N>,
        ccs: &CCS_F<'cfg, N>,
        transcript: &mut KeccakTranscript,
        config: &'cfg FieldConfig<N>,
    ) -> Result<VerificationPoints<'cfg, N>, SpartanError<N>>;
}

impl<'cfg, const I: usize, const N: usize, S: ZipSpec> SpartanVerifier<'cfg, N>
    for ZincVerifier<I, N, S>
{
    fn verify(
        &self,
        proof: &SpartanProof<'cfg, N>,
        ccs: &CCS_F<'cfg, N>,
        transcript: &mut KeccakTranscript,
        config: &'cfg FieldConfig<N>,
    ) -> Result<VerificationPoints<'cfg, N>, SpartanError<N>> {
        // Step 1: Generate the beta challenges.
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, ConfigRef::from(config));

        //Step 2: The sumcheck.
        let (r_x, s) =
            self.verify_linearization_proof(&proof.linearization_sumcheck, transcript, ccs)?;

        // Step 3. Check V_s is congruent to s
        Self::verify_linearization_claim(&beta_s, &r_x, s, proof, ccs)?;

        let gamma: RandomField<'cfg, N> =
            transcript.squeeze_gamma_challenge(ConfigRef::from(config));

        let second_sumcheck_claimed_sum = Self::lin_comb_V_s(&gamma, &proof.V_s);

        let (r_y, e_y) = self.verify_second_sumcheck_proof(
            &proof.second_sumcheck,
            transcript,
            ccs,
            second_sumcheck_claimed_sum,
        )?;

        Ok(VerificationPoints {
            rx_ry: [r_x, r_y].concat(),
            e_y,
            gamma,
        })
    }
}

impl<const I: usize, const N: usize, S: ZipSpec> ZincVerifier<I, N, S> {
    fn verify_linearization_proof<'cfg>(
        &self,
        proof: &SumcheckProof<'cfg, N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(Vec<RandomField<'cfg, N>>, RandomField<'cfg, N>), SpartanError<N>> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let nvars = ccs.s;
        let degree = ccs.d + 1;

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            nvars,
            degree,
            RandomField::zero(),
            proof,
            unsafe { ConfigRef::new(*ccs.config.as_ptr()) },
        )?;

        Ok((subclaim.point, subclaim.expected_evaluation))
    }

    fn verify_linearization_claim<'cfg>(
        beta_s: &[RandomField<'cfg, N>],
        point_r: &[RandomField<'cfg, N>],
        s: RandomField<'cfg, N>,
        proof: &SpartanProof<'cfg, N>,
        ccs: &CCS_F<'cfg, N>,
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
                should_equal_s.into(),
                s.into(),
            )));
        }

        Ok(())
    }

    fn verify_second_sumcheck_proof<'cfg>(
        &self,
        proof: &SumcheckProof<'cfg, N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<'cfg, N>,
        claimed_sum: RandomField<'cfg, N>,
    ) -> Result<(Vec<RandomField<'cfg, N>>, RandomField<'cfg, N>), SpartanError<N>> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let nvars = ccs.s_prime;
        let degree = 2;

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            nvars,
            degree,
            claimed_sum,
            proof,
            unsafe { ConfigRef::new(*ccs.config.as_ptr()) },
        )?;

        Ok((subclaim.point, subclaim.expected_evaluation))
    }

    fn lin_comb_V_s<'cfg>(
        gamma: &RandomField<'cfg, N>,
        V_s: &[RandomField<'cfg, N>],
    ) -> RandomField<'cfg, N> {
        let mut res = RandomField::zero();
        for V_i in V_s.iter().rev() {
            res *= gamma;
            res += V_i;
        }
        res
    }

    fn verify_pcs_proof<'cfg>(
        &self,
        cm_i: &Statement_F<'cfg, N>,
        zip_proof: &ZipProof<'cfg, I, N>,
        verification_points: &VerificationPoints<'cfg, N>,
        ccs: &CCS_F<'cfg, N>,
        transcript: &mut KeccakTranscript,
        config: &'cfg FieldConfig<N>,
    ) -> Result<(), SpartanError<N>>
    where
        [(); 2 * I]:,
        [(); 4 * I]:,
        [(); 8 * I]:,
    {
        let param =
            MultilinearZip::<I, { 2 * I }, { 4 * I }, { 8 * I }, S, KeccakTranscript>::setup(
                ccs.m, transcript,
            );
        let mut pcs_transcript = PcsTranscript::from_proof(&zip_proof.pcs_proof);
        let r_y = &verification_points.rx_ry[ccs.s..];

        MultilinearZip::<I, { 2 * I }, { 4 * I }, { 8 * I }, S, KeccakTranscript>::verify(
            &param,
            &zip_proof.z_comm,
            r_y,
            zip_proof.v,
            &mut pcs_transcript,
            config,
        )?;

        // Evaluate constraints at rx_ry point
        let V_xy = cm_i
            .constraints
            .iter()
            .map(|M| {
                let mle = DenseMultilinearExtension::from_matrix(M, ConfigRef::from(config));
                mle.evaluate(&verification_points.rx_ry, ConfigRef::from(config))
                    .ok_or(MleEvaluationError::IncorrectLength(
                        verification_points.rx_ry.len(),
                        mle.num_vars,
                    ))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Check final verification equation
        let V_x_gamma = Self::lin_comb_V_s(&verification_points.gamma, &V_xy) * zip_proof.v;
        if V_x_gamma != verification_points.e_y {
            return Err(SpartanError::PCSVerificationError(
                "linear combination of powers of gamma and V_x != e_y".to_string(),
            ));
        }

        Ok(())
    }
}

pub struct VerificationPoints<'cfg, const N: usize> {
    pub rx_ry: Vec<RandomField<'cfg, N>>,
    pub e_y: RandomField<'cfg, N>,
    pub gamma: RandomField<'cfg, N>,
}
