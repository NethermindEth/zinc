use ark_std::{boxed::Box, vec::Vec};
use num_traits::Zero;

use super::{
    errors::{MleEvaluationError, SpartanError, ZincError},
    structs::{SpartanProof, ZincProof, ZincVerifier, ZipProof},
    utils::{SqueezeBeta, SqueezeGamma, draw_random_field},
};
use crate::{
    ccs::{
        ccs_f::{CCS_F, Statement_F},
        ccs_z::{CCS_Z, Statement_Z},
    },
    field::RandomField,
    poly_f::mle::DenseMultilinearExtension,
    sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed, SumcheckProof, utils::eq_eval},
    traits::{ConfigReference, FieldMap, FromRef, Integer, MapsToField, ZipTypes},
    transcript::KeccakTranscript,
    zip::{
        code::LinearCodeSpec, code_raa::RaaCode, pcs::structs::MultilinearZip,
        pcs_transcript::PcsTranscript,
    },
};

pub trait Verifier<I: Integer, C: ConfigReference, S: LinearCodeSpec> {
    fn verify(
        &self,
        cm_i: &Statement_Z<I>,
        proof: ZincProof<C>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<I>,
        config: C,
    ) -> Result<(), ZincError>;
}

// TODO
impl<ZT: ZipTypes, C: ConfigReference, S: LinearCodeSpec> Verifier<ZT::N, C, S>
    for ZincVerifier<ZT, C, S>
where
    ZT::L: MapsToField<C>,
    ZT::K: MapsToField<C>,
    ZT::N: MapsToField<C>,
    ZT::N: FromRef<C::I>,
    C::I: FromRef<ZT::N>,
    C::I: FromRef<<ZT::N as Integer>::I>,
    Self: SpartanVerifier<C>,
{
    fn verify(
        &self,
        statement: &Statement_Z<ZT::N>,
        proof: ZincProof<C>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_Z<ZT::N>,
        config: C,
    ) -> Result<(), ZincError> {
        if draw_random_field::<ZT::N, C>(&statement.public_input, transcript)
            != *config.reference().unwrap()
        {
            return Err(ZincError::FieldConfigError);
        }
        // TODO: Write functionality to let the verifier know that there are no denominators that can be divided by q(As an honest prover)
        let ccs_F = ccs.map_to_field(config);
        let statement_f = statement.map_to_field(config);

        let verification_points =
            SpartanVerifier::verify(self, &proof.spartan_proof, &ccs_F, transcript, config)
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
pub trait SpartanVerifier<C: ConfigReference> {
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
        proof: &SpartanProof<C>,
        ccs: &CCS_F<C>,
        transcript: &mut KeccakTranscript,
        config: C,
    ) -> Result<VerificationPoints<C>, SpartanError>;
}

impl<ZT: ZipTypes, C: ConfigReference, S: LinearCodeSpec> SpartanVerifier<C>
    for ZincVerifier<ZT, C, S>
{
    fn verify(
        &self,
        proof: &SpartanProof<C>,
        ccs: &CCS_F<C>,
        transcript: &mut KeccakTranscript,
        config: C,
    ) -> Result<VerificationPoints<C>, SpartanError> {
        // Step 1: Generate the beta challenges.
        let beta_s = transcript.squeeze_beta_challenges(ccs.s, config);

        //Step 2: The sumcheck.
        let (r_x, s) =
            self.verify_linearization_proof(&proof.linearization_sumcheck, transcript, ccs)?;

        // Step 3. Check V_s is congruent to s
        Self::verify_linearization_claim(&beta_s, &r_x, s, proof, ccs)?;

        let gamma = transcript.squeeze_gamma_challenge(config);

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

impl<ZT: ZipTypes, C: ConfigReference, S: LinearCodeSpec> ZincVerifier<ZT, C, S> {
    fn verify_linearization_proof(
        &self,
        proof: &SumcheckProof<C>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<C>,
    ) -> Result<(Vec<RandomField<C>>, RandomField<C>), SpartanError> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let nvars = ccs.s;
        let degree = ccs.d + 1;

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            nvars,
            degree,
            RandomField::zero(),
            proof,
            unsafe { C::new(*ccs.config.as_ptr()) },
        )?;

        Ok((subclaim.point, subclaim.expected_evaluation))
    }

    fn verify_linearization_claim(
        beta_s: &[RandomField<C>],
        point_r: &[RandomField<C>],
        s: RandomField<C>,
        proof: &SpartanProof<C>,
        ccs: &CCS_F<C>,
    ) -> Result<(), SpartanError> {
        let e = eq_eval(point_r, beta_s)?;
        let should_equal_s = e * ccs // e * (\sum c_i * \Pi_{j \in S_i} u_j)
            .c
            .iter()
            .enumerate()
            .map(|(i, c)| {
                c.clone()
                    * ccs.S[i]
                        .iter()
                        .map(|&j| &proof.V_s[j])
                        .product::<RandomField<_>>()
            }) // c_i * \Pi_{j \in S_i} u_j
            .sum::<RandomField<_>>(); // \sum c_i * \Pi_{j \in S_i} u_j

        if should_equal_s != s {
            return Err(SpartanError::SumCheckError(SumCheckFailed(
                Box::new(should_equal_s.into()),
                Box::new(s.into()),
            )));
        }

        Ok(())
    }

    fn verify_second_sumcheck_proof(
        &self,
        proof: &SumcheckProof<C>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<C>,
        claimed_sum: RandomField<C>,
    ) -> Result<(Vec<RandomField<C>>, RandomField<C>), SpartanError> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let nvars = ccs.s_prime;
        let degree = 2;

        let subclaim = MLSumcheck::verify_as_subprotocol(
            transcript,
            nvars,
            degree,
            claimed_sum,
            proof,
            unsafe { C::new(*ccs.config.as_ptr()) },
        )?;

        Ok((subclaim.point, subclaim.expected_evaluation))
    }

    fn lin_comb_V_s(gamma: &RandomField<C>, V_s: &[RandomField<C>]) -> RandomField<C> {
        let mut res = RandomField::zero();
        for V_i in V_s.iter().rev() {
            res *= gamma;
            res += V_i;
        }
        res
    }

    fn verify_pcs_proof(
        &self,
        cm_i: &Statement_F<C>,
        zip_proof: &ZipProof<C>,
        verification_points: &VerificationPoints<C>,
        ccs: &CCS_F<C>,
        transcript: &mut KeccakTranscript,
        config: C,
    ) -> Result<(), SpartanError>
    where
        ZT::L: MapsToField<C>,
        ZT::K: MapsToField<C>,
    {
        let linear_code = RaaCode::<ZT>::new(&self.lc_spec, ccs.m, transcript);
        let param = MultilinearZip::<ZT, _>::setup(ccs.m, linear_code);
        let mut pcs_transcript = PcsTranscript::from_proof(&zip_proof.pcs_proof);
        let r_y = &verification_points.rx_ry[ccs.s..];

        MultilinearZip::<ZT, _>::verify(
            &param,
            &zip_proof.z_comm,
            r_y,
            zip_proof.v.clone(),
            &mut pcs_transcript,
            config,
        )?;

        // Evaluate constraints at rx_ry point
        let V_xy = cm_i
            .constraints
            .iter()
            .map(|M| {
                let mle = DenseMultilinearExtension::from_matrix(M, config);
                mle.evaluate(&verification_points.rx_ry, config).ok_or(
                    MleEvaluationError::IncorrectLength(
                        verification_points.rx_ry.len(),
                        mle.num_vars,
                    ),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Check final verification equation
        let V_x_gamma = Self::lin_comb_V_s(&verification_points.gamma, &V_xy) * &zip_proof.v;
        if V_x_gamma != verification_points.e_y {
            return Err(SpartanError::PCSVerificationError(
                "linear combination of powers of gamma and V_x != e_y".into(),
            ));
        }

        Ok(())
    }
}

pub struct VerificationPoints<C: ConfigReference> {
    pub rx_ry: Vec<RandomField<C>>,
    pub e_y: RandomField<C>,
    pub gamma: RandomField<C>,
}
