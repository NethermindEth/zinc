#![allow(dead_code, non_snake_case)]
use ark_ff::Zero;

use crate::{
    ccs::{ccs_f::CCS_F, error::CSError, utils::mat_vec_mul},
    field::RandomField,
    field_config::FieldConfig,
    poly::mle::DenseMultilinearExtension,
    sparse_matrix::SparseMatrix,
    sumcheck::utils::build_eq_x_r,
    transcript::KeccakTranscript,
};

use super::errors::{MleEvaluationError, SpartanError};

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
pub fn prepare_lin_sumcheck_polynomial<const N: usize>(
    c: &[RandomField<N>],
    d: &usize,
    M_mles: &[DenseMultilinearExtension<N>],
    S: &[Vec<usize>],
    beta_s: &[RandomField<N>],
    config: *const FieldConfig<N>,
) -> Result<(Vec<DenseMultilinearExtension<N>>, usize), SpartanError<N>> {
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

    mles.push(build_eq_x_r(beta_s, config)?);

    Ok((mles, d + 1))
}

pub(crate) fn sumcheck_polynomial_comb_fn_1<const N: usize>(
    vals: &[RandomField<N>],
    ccs: &CCS_F<N>,
) -> RandomField<N> {
    let mut result = RandomField::zero();
    'outer: for (i, &c) in ccs.c.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let mut term = c;
        for &j in &ccs.S[i] {
            if vals[j].is_zero() {
                continue 'outer;
            }
            term *= vals[j];
        }
        result += &term;
    }
    // eq() is the last term added
    result * vals[vals.len() - 1]
}

pub(crate) trait SqueezeBeta<const N: usize> {
    fn squeeze_beta_challenges(
        &mut self,
        n: usize,
        config: *const FieldConfig<N>,
    ) -> Vec<RandomField<N>>;
}

impl<const N: usize> SqueezeBeta<N> for KeccakTranscript {
    fn squeeze_beta_challenges(
        &mut self,
        n: usize,
        config: *const FieldConfig<N>,
    ) -> Vec<RandomField<N>> {
        self.absorb(b"beta_s");

        self.get_challenges(n, config)
    }
}

pub(crate) trait SqueezeGamma<const N: usize> {
    fn squeeze_gamma_challenge(&mut self, config: *const FieldConfig<N>) -> RandomField<N>;
}

impl<const N: usize> SqueezeGamma<N> for KeccakTranscript {
    fn squeeze_gamma_challenge(&mut self, config: *const FieldConfig<N>) -> RandomField<N> {
        self.absorb(b"gamma");

        self.get_challenge(config)
    }
}

// Prepare MLE's of the form mle[M_i \cdot z_ccs](x), a.k.a. \sum mle[M_i](x, b) * mle[z_ccs](b).
pub(super) fn calculate_Mz_mles<E, const N: usize>(
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
        constraints.iter().map(|M| mat_vec_mul(M, z_ccs)),
        config,
    )
}

fn to_mles_err<const N: usize, I, E, E1>(
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
