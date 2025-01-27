//! Prover

use ark_ff::Zero;
use ark_std::{cfg_into_iter, cfg_iter_mut, vec::Vec};
use lattirust_poly::{
    mle::MultilinearExtension,
    polynomials::{DenseMultilinearExtension, VirtualPolynomial},
};

use crate::field::RandomField;

use super::{verifier::VerifierMsg, IPForMLSumcheck};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Prover Message
#[derive(Clone, Debug, PartialEq)]
pub struct ProverMsg<const N: usize> {
    /// evaluations on P(0), P(1), P(2), ...
    pub(crate) evaluations: Vec<RandomField<N>>,
}

/// Prover State
pub struct ProverState<const N: usize> {
    /// sampled randomness given by the verifier
    pub randomness: Vec<RandomField<N>>,
    /// Stores the list of products that is meant to be added together. Each multiplicand is represented by
    /// the index in flattened_ml_extensions
    pub list_of_products: Vec<(RandomField<N>, Vec<usize>)>,
    /// Stores a list of multilinear extensions in which `self.list_of_products` points to
    pub flattened_ml_extensions: Vec<DenseMultilinearExtension<RandomField<N>>>,
    /// Number of variables
    pub num_vars: usize,
    /// Max number of multiplicands in a product
    pub max_multiplicands: usize,
    /// The current round number
    pub round: usize,
}

impl<const N: usize> IPForMLSumcheck<N> {
    /// initialize the prover to argue for the sum of polynomial over {0,1}^`num_vars`
    ///
    /// The polynomial is represented by a list of products of polynomials along with its coefficient that is meant to be added together.
    ///
    /// This data structure of the polynomial is a list of list of `(coefficient, DenseMultilinearExtension)`.
    /// * Number of products n = `polynomial.products.len()`,
    /// * Number of multiplicands of ith product m_i = `polynomial.products[i].1.len()`,
    /// * Coefficient of ith product c_i = `polynomial.products[i].0`
    ///
    /// The resulting polynomial is
    ///
    /// $$\sum_{i=0}^{n}C_i\cdot\prod_{j=0}^{m_i}P_{ij}$$
    ///
    pub fn prover_init(polynomial: &VirtualPolynomial<R>) -> ProverState<N> {
        if polynomial.aux_info.num_variables == 0 {
            panic!("Attempt to prove a constant.")
        }

        // create a deep copy of all unique MLExtensions
        let flattened_ml_extensions = ark_std::cfg_iter!(polynomial.flattened_ml_extensions)
            .map(|x| x.as_ref().clone())
            .collect();

        ProverState {
            randomness: Vec::with_capacity(polynomial.aux_info.num_variables),
            list_of_products: polynomial.products.clone(),
            flattened_ml_extensions,
            num_vars: polynomial.aux_info.num_variables,
            max_multiplicands: polynomial.aux_info.max_degree,
            round: 0,
        }
    }

    /// receive message from verifier, generate prover message, and proceed to next round
    ///
    /// Main algorithm used is from section 3.2 of [XZZPS19](https://eprint.iacr.org/2019/317.pdf#subsection.3.2).
    pub fn prove_round(
        prover_state: &mut ProverState<N>,
        v_msg: &Option<VerifierMsg<N>>,
    ) -> ProverMsg<N> {
        if let Some(msg) = v_msg {
            if prover_state.round == 0 {
                panic!("first round should be prover first.");
            }
            prover_state.randomness.push(msg.randomness);

            // fix argument
            let i = prover_state.round;
            let r = prover_state.randomness[i - 1];
            cfg_iter_mut!(prover_state.flattened_ml_extensions).for_each(|multiplicand| {
                *multiplicand = multiplicand.fix_variables(&[r.into()]);
            });
        } else if prover_state.round > 0 {
            panic!("verifier message is empty");
        }

        prover_state.round += 1;

        if prover_state.round > prover_state.num_vars {
            panic!("Prover is not active");
        }

        let i = prover_state.round;
        let nv = prover_state.num_vars;
        let degree = prover_state.max_multiplicands; // the degree of univariate polynomial sent by prover at this round

        #[cfg(not(feature = "parallel"))]
        let zeros = (
            vec![RandomField::<N>::zero(); degree + 1],
            vec![RandomField::<N>::zero(); degree + 1],
        );
        #[cfg(feature = "parallel")]
        let zeros = || (vec![R::zero(); degree + 1], vec![R::zero(); degree + 1]);

        // generate sum
        let fold_result =
            cfg_into_iter!(0..1 << (nv - i)).fold(zeros, |(mut products_sum, mut product), b| {
                // In effect, this fold is essentially doing simply:
                // for b in 0..1 << (nv - i) {
                let index = b << 1;
                for (coefficient, products) in &prover_state.list_of_products {
                    product.fill(*coefficient);
                    for &jth_product in products {
                        let table = &prover_state.flattened_ml_extensions[jth_product];
                        let mut start = table[index];
                        let step = table[index + 1] - start;
                        for p in product.iter_mut() {
                            *p *= start;
                            start += step;
                        }
                    }
                    for (sum, prod) in products_sum.iter_mut().zip(product.iter()) {
                        *sum += prod;
                    }
                }
                (products_sum, product)
            });

        #[cfg(not(feature = "parallel"))]
        let products_sum = fold_result.0;

        // When rayon is used, the `fold` operation results in a iterator of `Vec<F>` rather than a single `Vec<F>`. In this case, we simply need to sum them.
        #[cfg(feature = "parallel")]
        let products_sum = fold_result.map(|scratch| scratch.0).reduce(
            || vec![R::zero(); degree + 1],
            |mut overall_products_sum, sublist_sum| {
                overall_products_sum
                    .iter_mut()
                    .zip(sublist_sum.iter())
                    .for_each(|(f, s)| *f += s);
                overall_products_sum
            },
        );

        ProverMsg {
            evaluations: products_sum,
        }
    }
}
