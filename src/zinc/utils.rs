#![allow(non_snake_case)]

use ark_std::vec::Vec;

use super::errors::{MleEvaluationError, SpartanError};
use crate::{
    ccs::{CSError, CcsF, mat_vec_mul},
    poly::dense::DenseMultilinearExtension,
    sparse_matrix::SparseMatrix,
    sumcheck::utils::build_eq_x_r,
    traits::Field,
    transcript::KeccakTranscript,
};

/// Prepare the main linearization polynomial.
///
/// $$ g(\vec{\mathbf{x}}) := eq(\vec{\beta}, \vec{\mathbf{x}}) \cdot
/// \left(
/// \sum\_{i=1}^{n\_s} c\_i \cdot
/// \left[
/// \prod\_{j \in S\_i}
/// \left(
/// \sum\_{\vec{\mathbf{b}} \in \\{0,1\\}^{\log n\_c}}
/// \text{mle}[M\_j](\vec{\mathbf{x}}, \vec{\mathbf{b}}) \cdot \text{mle}\[\mathbf{z}\_{ccs}\](\vec{b})
/// \right)
/// \right]
/// \right) $$
///
/// # Parameters:
///
/// * `c` (`&[NTT]`): The second multiplicand of the polynomial is a linear combination of products of lists of MLEs, c is the coefficients of the lists
///
/// * `M_mles` (`&[DenseMultilinearExtension<NTT>]`): MLEs that the polynomial is constructed from
///
/// * `S` (`&[Vec<usize>]`): ] indices for the MLE lists
///
/// * `beta_s` (`&[NTT]`): Randomness
///
/// # Returns:
///
/// * The MLEs which form the polynomial
/// * The max degree of the polynomial
///
/// # Errors:
/// * Will return an error if any of the MLEs are of the wrong size
///
pub fn prepare_lin_sumcheck_polynomial<F: Field<LIMBS>, const LIMBS: usize>(
    c: &[F],
    d: &usize,
    M_mles: &[DenseMultilinearExtension<F>],
    S: &[Vec<usize>],
    beta_s: &[F],
) -> Result<(Vec<DenseMultilinearExtension<F>>, usize), SpartanError<F>> {
    // Flatten M_mles following S over non-zero c entries; eq(beta, x) is appended last.
    let len = 1 + c
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, _)| S[i].len())
        .sum::<usize>();

    let mut mles = Vec::with_capacity(len);

    for (i, _) in c.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
        for &j in &S[i] {
            mles.push(M_mles[j].clone());
        }
    }

    mles.push(build_eq_x_r(beta_s)?);

    Ok((mles, d + 1))
}

pub(crate) fn sumcheck_polynomial_comb_fn_1<F: Field<LIMBS>, const LIMBS: usize>(
    vals: &[F],
    ccs: &CcsF<F>,
) -> F {
    // vals layout must match the flattened order constructed in prepare_lin_sumcheck_polynomial
    let mut result = F::zero();
    let mut k = 0usize;

    'outer: for (i, c) in ccs.c.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let mut term = *c;
        for _ in &ccs.S[i] {
            let v = vals[k];
            k += 1;
            if v.is_zero() {
                continue 'outer;
            }
            term = term * v;
        }
        result += term;
    }

    // eq() is the last term added
    result * vals[vals.len() - 1]
}

pub(crate) trait SqueezeBeta<F: Field<LIMBS>, const LIMBS: usize> {
    fn squeeze_beta_challenges(&mut self, n: usize) -> Vec<F>;
}

impl<F: Field<LIMBS>, const LIMBS: usize> SqueezeBeta<F, LIMBS> for KeccakTranscript {
    fn squeeze_beta_challenges(&mut self, n: usize) -> Vec<F> {
        self.absorb(b"beta_s");

        self.get_challenges(n)
    }
}

pub(crate) trait SqueezeGamma<F: Field<LIMBS>, const LIMBS: usize> {
    fn squeeze_gamma_challenge(&mut self) -> F;
}

impl<F: Field<LIMBS>, const LIMBS: usize> SqueezeGamma<F, LIMBS> for KeccakTranscript {
    fn squeeze_gamma_challenge(&mut self) -> F {
        self.absorb(b"gamma");

        self.get_challenge()
    }
}

// Prepare MLE's of the form mle[M_i \cdot z_ccs](x), a.k.a. \sum mle[M_i](x, b) * mle[z_ccs](b).
pub(super) fn calculate_Mz_mles<E, F: Field<LIMBS>, const LIMBS: usize>(
    constraints: &[SparseMatrix<F>],
    ccs_s: usize,
    z_ccs: &[F],
) -> Result<Vec<DenseMultilinearExtension<F>>, E>
where
    E: From<MleEvaluationError> + From<CSError> + Sync + Send,
{
    to_mles_err::<F, LIMBS, _, E, CSError>(ccs_s, constraints.iter().map(|M| mat_vec_mul(M, z_ccs)))
}

fn to_mles_err<F: Field<LIMBS>, const LIMBS: usize, I, E, E1>(
    n_vars: usize,
    mle_s: I,
) -> Result<Vec<DenseMultilinearExtension<F>>, E>
where
    I: IntoIterator<Item = Result<Vec<F>, E1>>,
    E: From<MleEvaluationError> + From<E1>,
{
    mle_s
        .into_iter()
        .map(|m| {
            let m = m?;
            if 1 << n_vars < m.len() {
                Err(MleEvaluationError::IncorrectLength(1 << n_vars, m.len()).into())
            } else {
                Ok(DenseMultilinearExtension::from_evaluations_vec(n_vars, m))
            }
        })
        .collect::<Result<_, E>>()
}

// pub fn draw_random_field<const I: usize, F: Field<LIMBS>, const LIMBS: usize>(
//     public_inputs: &[Int<I>],
//     transcript: &mut KeccakTranscript,
// ) -> F::C {
//     for input in public_inputs {
//         transcript.absorb(cast_slice(input.as_words()));
//     }
//     // Method for efficient random prime sampling not yet implemented
//     // Fixing the random prime q for now
//     F::C::new(get_prime::<F>(transcript))
// }
