//! Verifier
use crate::field_config::ConfigRef;
use ark_ff::One;
use ark_std::vec::Vec;

use super::{prover::ProverMsg, IPForMLSumcheck, SumCheckError};
use crate::{
    field::{conversion::FieldMap, RandomField},
    transcript::KeccakTranscript as Transcript,
};

pub const SQUEEZE_NATIVE_ELEMENTS_NUM: usize = 1;

/// Verifier Message
#[derive(Clone)]
pub struct VerifierMsg<'cfg, const N: usize> {
    /// randomness sampled by verifier
    pub randomness: RandomField<'cfg, N>,
}

/// Verifier State
pub struct VerifierState<'cfg, const N: usize> {
    round: usize,
    nv: usize,
    max_multiplicands: usize,
    finished: bool,
    /// a list storing the univariate polynomial in evaluation form sent by the prover at each round
    polynomials_received: Vec<Vec<RandomField<'cfg, N>>>,
    /// a list storing the randomness sampled by the verifier at each round
    randomness: Vec<RandomField<'cfg, N>>,
    /// The configuration of the field that the sumcheck protocol is working in
    config: ConfigRef<'cfg, N>,
}

/// Subclaim when verifier is convinced
#[derive(Debug)]
pub struct SubClaim<'cfg, const N: usize> {
    /// the multi-dimensional point that this multilinear extension is evaluated to
    pub point: Vec<RandomField<'cfg, N>>,
    /// the expected evaluation
    pub expected_evaluation: RandomField<'cfg, N>,
}

impl<const N: usize> IPForMLSumcheck<N> {
    /// initialize the verifier
    pub fn verifier_init(nvars: usize, degree: usize, config: ConfigRef<N>) -> VerifierState<N> {
        VerifierState {
            round: 1,
            nv: nvars,
            max_multiplicands: degree,
            finished: false,
            polynomials_received: Vec::with_capacity(nvars),
            randomness: Vec::with_capacity(nvars),
            config,
        }
    }

    /// Run verifier at current round, given prover message
    ///
    /// Normally, this function should perform actual verification. Instead, `verify_round` only samples
    /// and stores randomness and perform verifications altogether in `check_and_generate_subclaim` at
    /// the last step.
    pub fn verify_round<'cfg>(
        prover_msg: ProverMsg<'cfg, N>,
        verifier_state: &mut VerifierState<'cfg, N>,
        transcript: &mut Transcript,
    ) -> VerifierMsg<'cfg, N> {
        if verifier_state.finished {
            panic!("Incorrect verifier state: Verifier is already finished.");
        }

        // Now, verifier should check if the received P(0) + P(1) = expected. The check is moved to
        // `check_and_generate_subclaim`, and will be done after the last round.

        let msg = Self::sample_round(transcript, verifier_state.config);
        verifier_state.randomness.push(msg.randomness);
        verifier_state
            .polynomials_received
            .push(prover_msg.evaluations);

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
    pub fn check_and_generate_subclaim<'cfg>(
        verifier_state: VerifierState<'cfg, N>,
        asserted_sum: RandomField<'cfg, N>,
        config: ConfigRef<'cfg, N>,
    ) -> Result<SubClaim<'cfg, N>, SumCheckError<N>> {
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
                panic!("incorrect number of evaluations");
            }

            let p0 = evaluations[0];
            let p1 = evaluations[1];
            if p0 + p1 != expected {
                return Err(SumCheckError::SumCheckFailed(
                    (p0 + p1).into(),
                    expected.into(),
                ));
            }
            expected = interpolate_uni_poly(evaluations, verifier_state.randomness[i], config);
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
    pub fn sample_round<'cfg>(
        transcript: &mut Transcript,
        config: ConfigRef<'cfg, N>,
    ) -> VerifierMsg<'cfg, N> {
        VerifierMsg {
            randomness: transcript.get_challenge(config),
        }
    }
}

/// interpolate the *unique* univariate polynomial of degree *at most*
/// p_i.len()-1 passing through the y-values in p_i at x = 0,..., p_i.len()-1
/// and evaluate this  polynomial at `eval_at`. In other words, efficiently compute
///  \sum_{i=0}^{len p_i - 1} p_i[i] * (\prod_{j!=i} (eval_at - j)/(i-j))
pub(crate) fn interpolate_uni_poly<'cfg, const N: usize>(
    p_i: &[RandomField<'cfg, N>],
    x: RandomField<'cfg, N>,
    config: ConfigRef<'cfg, N>,
) -> RandomField<'cfg, N> {
    // We will need these a few times
    let zero = 0u64.map_to_field(config);
    let one = 1u64.map_to_field(config);

    let len = p_i.len();

    let mut evals = vec![];

    let mut prod = x;
    evals.push(x);

    //`prod = \prod_{j} (x - j)`
    // we return early if 0 <= x < len, i.e. if the desired value has been passed
    let mut j = zero;
    for i in 1..len {
        if x == j {
            return p_i[i - 1];
        }
        j += &one;

        let tmp = x - j;
        evals.push(tmp);
        prod *= tmp;
    }

    if x == j {
        return p_i[len - 1];
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
        let mut last_denom = u64_factorial(len - 1).map_to_field(config);

        last_denom.set_config(config);

        let mut ratio_numerator = 1i64;
        let mut ratio_enumerator = 1u64;

        for i in (0..len).rev() {
            let ratio_numerator_f = if ratio_numerator < 0 {
                let mut res = (-ratio_numerator as u64).map_to_field(config);
                res.set_config(config);
                -res
            } else {
                let mut res = (ratio_numerator as u64).map_to_field(config);
                res.set_config(config);
                res
            };

            let mut ratio_enumerator_f = ratio_enumerator.map_to_field(config);
            ratio_enumerator_f.set_config(config);

            let x = prod * ratio_enumerator_f / (last_denom * ratio_numerator_f * evals[i]);

            res += &(p_i[i] * x);

            // compute ratio for the next step which is current_ratio * -(len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i64 - i as i64);
                ratio_enumerator *= i as u64;
            }
        }
    } else if p_i.len() <= 33 {
        let last_denom = u128_factorial(len - 1).map_to_field(config);
        let mut ratio_numerator = 1i128;
        let mut ratio_enumerator = 1u128;

        for i in (0..len).rev() {
            let ratio_numerator_f = if ratio_numerator < 0 {
                let mut res = (-ratio_numerator as u128).map_to_field(config);
                res.set_config(config);
                -res
            } else {
                let mut res = (ratio_numerator as u128).map_to_field(config);
                res.set_config(config);
                res
            };

            let mut ratio_enumerator_f = ratio_enumerator.map_to_field(config);
            ratio_enumerator_f.set_config(config);

            let x: RandomField<N> =
                prod * ratio_enumerator_f / (last_denom * ratio_numerator_f * evals[i]);
            res += &(p_i[i] * x);

            // compute ratio for the next step which is current_ratio * -(len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i128 - i as i128);
                ratio_enumerator *= i as u128;
            }
        }
    } else {
        // since we are using field operations, we can merge
        // `last_denom` and `ratio_numerator` into a single field element.
        let mut denom_up = field_factorial::<N>(len - 1, config);
        let mut denom_down = one;

        for i in (0..len).rev() {
            let x = prod * denom_down / (denom_up * evals[i]);
            res += &(p_i[i] * x);

            // compute denom for the next step is -current_denom * (len-i)/i
            if i != 0 {
                let mut denom_up_factor = ((len - i) as u64).map_to_field(config);
                denom_up_factor.set_config(config);
                denom_up *= -denom_up_factor;

                let mut denom_down_factor = (i as u64).map_to_field(config);
                denom_down_factor.set_config(config);
                denom_down *= denom_down_factor;
            }
        }
    }

    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn field_factorial<const N: usize>(a: usize, config: ConfigRef<N>) -> RandomField<N> {
    let mut res = RandomField::one();
    for i in 1..=a {
        res *= (i as u64).map_to_field(config);
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

// #[cfg(test)]
// mod test {
//     use super::*;

//     use ark_poly::univariate::DensePolynomial;
//     use ark_poly::DenseUVPolynomial;
//     use ark_poly::Polynomial;
//     use ark_std::vec::Vec;
//     use ark_std::UniformRand;

//     type F = ark_test_curves::bls12_381::Fr;

//     #[test]
//     fn test_interpolation() {
//         let mut prng = ark_std::test_rng();

//         // test a polynomial with 20 known points, i.e., with degree 19
//         let poly = DensePolynomial::<F>::rand(20 - 1, &mut prng);
//         let evals = (0..20)
//             .map(|i| poly.evaluate(&F::from(i)))
//             .collect::<Vec<F>>();
//         let query = F::rand(&mut prng);

//         assert_eq!(poly.evaluate(&query), interpolate_uni_poly(&evals, query));

//         // test a polynomial with 33 known points, i.e., with degree 32
//         let poly = DensePolynomial::<F>::rand(33 - 1, &mut prng);
//         let evals = (0..33)
//             .map(|i| poly.evaluate(&F::from(i)))
//             .collect::<Vec<F>>();
//         let query = F::rand(&mut prng);

//         assert_eq!(poly.evaluate(&query), interpolate_uni_poly(&evals, query));

//         // test a polynomial with 64 known points, i.e., with degree 63
//         let poly = DensePolynomial::<F>::rand(64 - 1, &mut prng);
//         let evals = (0..64)
//             .map(|i| poly.evaluate(&F::from(i)))
//             .collect::<Vec<F>>();
//         let query = F::rand(&mut prng);

//         assert_eq!(poly.evaluate(&query), interpolate_uni_poly(&evals, query));

//         // test interpolation when we ask for the value at an x-cordinate
//         // we are already passing, i.e. in the range 0 <= x < len(values) - 1
//         let evals = vec![0, 1, 4, 9]
//             .into_iter()
//             .map(F::from)
//             .collect::<Vec<F>>();
//         assert_eq!(interpolate_uni_poly(&evals, F::from(3)), F::from(9));
//     }
// }
