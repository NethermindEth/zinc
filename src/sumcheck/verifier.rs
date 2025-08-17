//! Verifier

use crypto_bigint::NonZero;

use crate::{
    sumcheck::{IPForMLSumcheck, SumCheckError, prover::ProverMsg},
    traits::{Field, Ring},
    transcript::KeccakTranscript,
};

pub const SQUEEZE_NATIVE_ELEMENTS_NUM: usize = 1;

/// Verifier Message
#[derive(Clone)]
pub struct VerifierMsg<F> {
    /// randomness sampled by verifier
    pub randomness: F,
}

/// Verifier State
pub struct VerifierState<F> {
    round: usize,
    nv: usize,
    max_multiplicands: usize,
    finished: bool,
    /// a list storing the univariate polynomial in evaluation form sent by the prover at each round
    polynomials_received: Vec<Vec<F>>,
    /// a list storing the randomness sampled by the verifier at each round
    randomness: Vec<F>,
}

/// Subclaim when verifier is convinced
#[derive(Debug)]
pub struct SubClaim<F> {
    /// the multi-dimensional point that this multilinear extension is evaluated to
    pub point: Vec<F>,
    /// the expected evaluation
    pub expected_evaluation: F,
}

impl<F: Field<LIMBS>, const LIMBS: usize> IPForMLSumcheck<F, LIMBS> {
    /// initialize the verifier
    pub fn verifier_init(nvars: usize, degree: usize) -> VerifierState<F> {
        VerifierState {
            round: 1,
            nv: nvars,
            max_multiplicands: degree,
            finished: false,
            polynomials_received: Vec::with_capacity(nvars),
            randomness: Vec::with_capacity(nvars),
        }
    }

    /// Run verifier at current round, given prover message
    ///
    /// Normally, this function should perform actual verification. Instead, `verify_round` only samples
    /// and stores randomness and perform verifications altogether in `check_and_generate_subclaim` at
    /// the last step.
    pub fn verify_round(
        prover_msg: &ProverMsg<F>,
        verifier_state: &mut VerifierState<F>,
        transcript: &mut KeccakTranscript,
    ) -> VerifierMsg<F> {
        if verifier_state.finished {
            panic!("Incorrect verifier state: Verifier is already finished.");
        }

        // Now, verifier should check if the received P(0) + P(1) = expected. The check is moved to
        // `check_and_generate_subclaim`, and will be done after the last round.

        let msg = Self::sample_round(transcript);
        verifier_state.randomness.push(msg.randomness);
        verifier_state
            .polynomials_received
            .push(prover_msg.evaluations.clone());

        // Now, verifier should set `expected` to P(r).
        // This operation is also moved to `check_and_generate_subclaim`,
        // and will be done after the last round.

        if verifier_state.round == verifier_state.nv {
            // accept and close
            verifier_state.finished = true;
        } else {
            verifier_state.round += 1;
        }
        msg
    }

    /// verify the sumcheck phase, and generate the subclaim
    ///
    /// If the asserted sum is correct, then the multilinear polynomial evaluated at `subclaim.point`
    /// is `subclaim.expected_evaluation`. Otherwise, it is highly unlikely that those two will be equal.
    /// Larger field size guarantees smaller soundness error.
    pub fn check_and_generate_subclaim(
        verifier_state: VerifierState<F>,
        asserted_sum: F,
    ) -> Result<SubClaim<F>, SumCheckError<F>> {
        if !verifier_state.finished {
            panic!("Verifier has not finished.");
        }

        let mut expected = asserted_sum;
        if verifier_state.polynomials_received.len() != verifier_state.nv {
            panic!("insufficient rounds");
        }
        for i in 0..verifier_state.nv {
            let evaluations = &verifier_state.polynomials_received[i];
            if evaluations.len() != verifier_state.max_multiplicands + 1 {
                return Err(SumCheckError::MaxDegreeExceeded);
            }

            let p0 = &evaluations[0];
            if verifier_state.max_multiplicands > 0 {
                let p1 = &evaluations[1];
                if *p0 + *p1 != expected {
                    return Err(SumCheckError::SumCheckFailed(
                        Box::new(*p0 + *p1),
                        Box::new(expected),
                    ));
                }
            } else {
                // Degree 0, constant polynomial
                if *p0 != expected {
                    return Err(SumCheckError::SumCheckFailed(
                        Box::new(*p0),
                        Box::new(expected),
                    ));
                }
            }

            expected = interpolate_uni_poly(evaluations, verifier_state.randomness[i]);
        }

        Ok(SubClaim {
            point: verifier_state.randomness,
            expected_evaluation: expected,
        })
    }

    /// simulate a verifier message without doing verification
    ///
    /// Given the same calling context, `transcript_round` output exactly the same message as
    /// `verify_round`
    #[inline]
    pub fn sample_round(transcript: &mut KeccakTranscript) -> VerifierMsg<F> {
        VerifierMsg {
            randomness: transcript.get_challenge(),
        }
    }
}

/// interpolate the *unique* univariate polynomial of degree *at most*
/// p_i.len()-1 passing through the y-values in p_i at x = 0,..., p_i.len()-1
/// and evaluate this  polynomial at `eval_at`. In other words, efficiently compute
///  \sum_{i=0}^{len p_i - 1} p_i[i] * (\prod_{j!=i} (eval_at - j)/(i-j))
pub(crate) fn interpolate_uni_poly<F: Ring + From<u64> + From<i64> + From<i128> + From<u128>>(
    p_i: &[F],
    x: F,
) -> F {
    // We will need these a few times
    let zero = F::ZERO;
    let one = F::one();

    let len = p_i.len();

    // Handle edge cases first
    if len == 0 {
        return zero;
    }
    if len == 1 {
        return p_i[0];
    }

    let mut evals = Vec::with_capacity(len);

    // Build evals = [x-0, x-1, ..., x-(len-1)] and prod = prod_j (x - j)
    let mut prod = one;
    let mut j = zero;
    for coeff in p_i {
        let tmp = x - j;
        // Early return if x equals some integer j in [0, len-1]
        if tmp.is_zero() {
            return *coeff;
        }
        evals.push(tmp);
        prod = prod * tmp;
        j += &one;
    }

    let mut res = zero;
    // we want to compute \prod (j!=i) (i-j) for a given i
    //
    // we start from the last step, which is
    //  denom[len-1] = (len-1) * (len-2) *... * 2 * 1
    // the step before that is
    //  denom[len-2] = (len-2) * (len-3) * ... * 2 * 1 * -1
    // and the step before that is
    //  denom[len-3] = (len-3) * (len-4) * ... * 2 * 1 * -1 * -2
    //
    // i.e., for any i, the one before this will be derived from
    //  denom[i-1] = - denom[i] * (len-i) / i
    //
    // that is, we only need to store
    // - the last denom for i = len-1, and
    // - the ratio between the current step and the last step, which is the
    //   product of -(len-i) / i from all previous steps and we store
    //   this product as a fraction number to reduce field divisions.

    // We know
    //  - 2^61 < factorial(20) < 2^62
    //  - 2^122 < factorial(33) < 2^123
    // so we will be able to compute the ratio
    //  - for len <= 20 with i64
    //  - for len <= 33 with i128
    //  - for len >  33 with BigInt
    if p_i.len() <= 20 {
        let last_denom = F::from(u64_factorial(len - 1));

        let mut ratio_numerator = 1i64;
        let mut ratio_enumerator = 1u64;

        for i in (0..len).rev() {
            let ratio_numerator_f = F::from(ratio_numerator);
            let ratio_enumerator_f = F::from(ratio_enumerator);

            let x = prod * ratio_enumerator_f
                / NonZero::new(last_denom * ratio_numerator_f * evals[i]).unwrap();

            res += &(p_i[i] * x.unwrap());

            // compute ratio for the next step which is current_ratio * -(len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i64 - i as i64);
                ratio_enumerator *= i as u64;
            }
        }
    } else if p_i.len() <= 33 {
        let last_denom = F::from(u128_factorial(len - 1));
        let mut ratio_numerator = 1i128;
        let mut ratio_enumerator = 1u128;

        for i in (0..len).rev() {
            let ratio_numerator_f = F::from(ratio_numerator);

            let ratio_enumerator_f = F::from(ratio_enumerator);

            let x = prod * ratio_enumerator_f
                / NonZero::new(last_denom * ratio_numerator_f * evals[i]).unwrap();
            res += &(p_i[i] * x.unwrap());

            // compute ratio for the next step which is current_ratio * -(len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i128 - i as i128);
                ratio_enumerator *= i as u128;
            }
        }
    } else {
        // since we are using field operations, we can merge
        // `last_denom` and `ratio_numerator` into a single field element.
        let mut denom_up = field_factorial::<F>(len - 1);
        let mut denom_down = one;

        for i in (0..len).rev() {
            let x = prod * denom_down / NonZero::new(denom_up * evals[i]).unwrap();
            res += &(p_i[i] * x.unwrap());

            // compute denom for the next step is -current_denom * (len-i)/i
            if i != 0 {
                let denom_up_factor = F::from((len - i) as u64);
                denom_up = denom_up * denom_up_factor.negate();

                let denom_down_factor = F::from(i as u64);
                denom_down = denom_down * denom_down_factor;
            }
        }
    }

    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn field_factorial<F: Ring + From<u64>>(a: usize) -> F {
    let mut res: F = F::one();
    for i in 1..=(a as u64) {
        res = res * F::from(i);
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u128_factorial(a: usize) -> u128 {
    let mut res = 1u128;
    for i in 1..=a {
        res *= i as u128;
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u64_factorial(a: usize) -> u64 {
    let mut res = 1u64;
    for i in 1..=a {
        res *= i as u64;
    }
    res
}
