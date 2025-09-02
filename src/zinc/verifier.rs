use ark_std::{boxed::Box, vec::Vec};
use crypto_bigint::Int;

use super::{
    errors::{MleEvaluationError, SpartanError, ZincError},
    structs::{SpartanProof, ZincProof, ZincVerifier, ZipProof},
    utils::{SqueezeBeta, SqueezeGamma},
};
use crate::{
    ccs::{CcsF, CcsZ, Statement},
    poly::dense::DenseMultilinearExtension,
    sumcheck::{MLSumcheck, SumCheckError::SumCheckFailed, SumcheckProof},
    traits::Field,
    transcript::KeccakTranscript,
    zip::{
        code::LinearCodeSpec, code_raa::RaaCode, pcs::structs::MultilinearZip,
        pcs_transcript::PcsTranscript,
    },
};

pub trait Verifier<const I: usize, F: Field<LIMBS>, const LIMBS: usize, S: LinearCodeSpec> {
    fn verify(
        &self,
        cm_i: &Statement<Int<I>>,
        proof: ZincProof<F>,
        transcript: &mut KeccakTranscript,
        ccs: &CcsZ<Int<I>>,
    ) -> Result<(), ZincError<F>>;
}

// TODO
impl<
    const N: usize,
    const L: usize,
    const K: usize,
    const M: usize,
    F: Field<LIMBS>,
    const LIMBS: usize,
    S: LinearCodeSpec,
> Verifier<N, F, LIMBS, S> for ZincVerifier<N, L, K, M, F, LIMBS, S>
where
    Self: SpartanVerifier<F, LIMBS>,
{
    fn verify(
        &self,
        statement: &Statement<Int<N>>,
        proof: ZincProof<F>,
        transcript: &mut KeccakTranscript,
        ccs: &CcsZ<Int<N>>,
    ) -> Result<(), ZincError<F>> {
        // TODO: Write functionality to let the verifier know that there are no denominators that can be divided by q(As an honest prover)
        let ccs_F = ccs.map_to_field();
        let statement_f = statement.map_to_field();

        let verification_points =
            SpartanVerifier::<F, LIMBS>::verify(self, &proof.spartan_proof, &ccs_F, transcript)
                .map_err(ZincError::SpartanError)?;

        self.verify_pcs_proof(
            &statement_f,
            &proof.zip_proof,
            &verification_points,
            &ccs_F,
            transcript,
        )?;

        Ok(())
    }
}

/// Verifier for the Linearization subprotocol.
pub trait SpartanVerifier<F: Field<LIMBS>, const LIMBS: usize> {
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
        proof: &SpartanProof<F>,
        ccs: &CcsF<F>,
        transcript: &mut KeccakTranscript,
    ) -> Result<VerificationPoints<F>, SpartanError<F>>;
}

impl<
    const N: usize,
    const L: usize,
    const K: usize,
    const M: usize,
    F: Field<LIMBS>,
    const LIMBS: usize,
    S: LinearCodeSpec,
> SpartanVerifier<F, LIMBS> for ZincVerifier<N, L, K, M, F, LIMBS, S>
{
    fn verify(
        &self,
        proof: &SpartanProof<F>,
        ccs: &CcsF<F>,
        transcript: &mut KeccakTranscript,
    ) -> Result<VerificationPoints<F>, SpartanError<F>> {
        // Step 1: Generate the beta challenges.
        let beta_s = transcript.squeeze_beta_challenges(ccs.s);

        //Step 2: The sumcheck.
        let (r_x, s) =
            self.verify_linearization_proof(&proof.linearization_sumcheck, transcript, ccs)?;

        // Step 3. Check V_s is congruent to s
        Self::verify_linearization_claim(&beta_s, &r_x, s, proof, ccs)?;

        let gamma: F = transcript.squeeze_gamma_challenge();

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

impl<
    const N: usize,
    const L: usize,
    const K: usize,
    const M: usize,
    F: Field<LIMBS>,
    const LIMBS: usize,
    S: LinearCodeSpec,
> ZincVerifier<N, L, K, M, F, LIMBS, S>
{
    fn verify_linearization_proof(
        &self,
        proof: &SumcheckProof<F>,
        transcript: &mut KeccakTranscript,
        ccs: &CcsF<F>,
    ) -> Result<(Vec<F>, F), SpartanError<F>> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let nvars = ccs.s;
        let degree = ccs.d + 1;

        let claimed_sum = MLSumcheck::extract_sum(proof);
        let subclaim =
            MLSumcheck::verify_as_subprotocol(transcript, nvars, degree, claimed_sum, proof)?;

        Ok((subclaim.point, subclaim.expected_evaluation))
    }

    fn verify_linearization_claim(
        beta_s: &[F],
        point_r: &[F],
        s: F,
        proof: &SpartanProof<F>,
        ccs: &CcsF<F>,
    ) -> Result<(), SpartanError<F>> {
        let e = {
            // Evaluate eq(beta_s, x) at x = point_r using the same routine
            // used to build the eq MLE, to avoid any subtle inconsistencies.
            let eq_mle = crate::sumcheck::utils::build_eq_x_r(beta_s)?;
            eq_mle
                .evaluate(point_r)
                .ok_or(MleEvaluationError::IncorrectLength(
                    point_r.len(),
                    eq_mle.num_vars,
                ))?
        };
        let should_equal_s = e * ccs // e * (\sum c_i * \Pi_{j \in S_i} u_j)
            .c
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let prod = ccs.S[i].iter().fold(F::one(), |acc, &j| acc * proof.V_s[j]);
                *c * prod
            }) // c_i * \Pi_{j \in S_i} u_j
            .fold(F::zero(), |acc, term| acc + term); // \sum c_i * \Pi_{j \in S_i} u_j

        if should_equal_s != s {
            return Err(SpartanError::SumCheckError(SumCheckFailed(
                Box::new(should_equal_s),
                Box::new(s),
            )));
        }

        Ok(())
    }

    fn verify_second_sumcheck_proof(
        &self,
        proof: &SumcheckProof<F>,
        transcript: &mut KeccakTranscript,
        ccs: &CcsF<F>,
        claimed_sum: F,
    ) -> Result<(Vec<F>, F), SpartanError<F>> {
        // The polynomial has degree <= ccs.d + 1 and log_m (ccs.s) vars.
        let nvars = ccs.s_prime;
        let degree = 2;

        let subclaim =
            MLSumcheck::verify_as_subprotocol(transcript, nvars, degree, claimed_sum, proof)?;

        Ok((subclaim.point, subclaim.expected_evaluation))
    }

    fn lin_comb_V_s(gamma: &F, V_s: &[F]) -> F {
        let mut res = F::zero();
        for V_i in V_s.iter().rev() {
            res = res * gamma;
            res += V_i;
        }
        res
    }

    fn verify_pcs_proof(
        &self,
        cm_i: &Statement<F>,
        zip_proof: &ZipProof<F>,
        verification_points: &VerificationPoints<F>,
        ccs: &CcsF<F>,
        transcript: &mut KeccakTranscript,
    ) -> Result<(), SpartanError<F>> {
        let linear_code = RaaCode::<N, L, K, M>::new(&self.lc_spec, ccs.m, transcript);
        let param = MultilinearZip::<N, L, K, M, _>::setup(ccs.m, linear_code);
        let mut pcs_transcript = PcsTranscript::from_proof(&zip_proof.pcs_proof);
        let r_y = &verification_points.rx_ry[ccs.s..];

        MultilinearZip::<N, L, K, M, _>::verify(
            &param,
            &zip_proof.z_comm,
            r_y,
            zip_proof.v,
            &mut pcs_transcript,
        )?;

        // Evaluate constraints at rx_ry point
        let V_xy = cm_i
            .constraints
            .iter()
            .map(|mat| {
                let mle = DenseMultilinearExtension::from_matrix(mat);
                mle.evaluate(&verification_points.rx_ry)
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
                "linear combination of powers of gamma and V_x != e_y".into(),
            ));
        }

        Ok(())
    }
}

#[derive(Debug)]
pub struct VerificationPoints<F> {
    pub rx_ry: Vec<F>,
    pub e_y: F,
    pub gamma: F,
}
