// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// Adapted for rings by Nethermind

//! This module defines our main mathematical object `DensePolynomial`; and
//! various functions associated with it.
use ark_std::{
    cfg_iter_mut, end_timer,
    rand::{Rng, RngCore},
    start_timer, vec,
    vec::Vec,
};
#[cfg(feature = "parallel")]
use rayon::iter::*;

use crate::{
    poly::ArithErrors,
    poly_f::{
        mle::DenseMultilinearExtension,
        polynomials::{random_mle_list, RefCounter},
    },
    traits::{Field, FieldMap},
};

#[allow(clippy::type_complexity)]
pub fn rand_poly<F: Field>(
    nv: usize,
    num_multiplicands_range: (usize, usize),
    num_products: usize,
    config: F::R,
    rng: &mut impl RngCore,
) -> Result<
    (
        (Vec<DenseMultilinearExtension<F>>, usize),
        Vec<(F, Vec<usize>)>,
        F,
    ),
    ArithErrors,
> {
    let mut sum = 0u64.map_to_field(config);
    let mut mles = vec![];
    let mut products = Vec::with_capacity(num_products);
    let mut degree = 0;
    let mut current_mle_index = 0;
    for _ in 0..num_products {
        let num_multiplicands = rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
        degree = num_multiplicands.max(degree);
        let (product, product_sum) = random_mle_list(nv, num_multiplicands, rng, config);
        let product = product
            .into_iter()
            .map(|p| RefCounter::into_inner(p).unwrap())
            .collect::<Vec<_>>();

        let coefficient = F::rand_with_config(rng, config);
        mles.extend(product);
        sum += &(product_sum * &coefficient);

        let indices: Vec<usize> =
            (current_mle_index..current_mle_index + num_multiplicands).collect();
        products.push((coefficient, indices));
        current_mle_index += num_multiplicands;
    }

    Ok(((mles, degree), products, sum))
}

pub fn rand_poly_comb_fn<F: Field>(vals: &[F], products: &[(F, Vec<usize>)], config: F::R) -> F {
    let mut result = 0u64.map_to_field(config);
    for (coef, indices) in products {
        let mut term = coef.clone();
        for &i in indices {
            term *= &vals[i];
        }
        result += &term;
    }

    result
}

/// Evaluate eq polynomial.
pub fn eq_eval<F: Field>(x: &[F], y: &[F]) -> Result<F, ArithErrors> {
    if x.len() != y.len() {
        return Err(ArithErrors::InvalidParameters(
            "x and y have different length".into(),
        ));
    }
    let start = start_timer!(|| "eq_eval");
    let mut res = F::one();
    for (xi, yi) in x.iter().zip(y.iter()) {
        let xi_yi = xi.clone() * yi;
        res *= xi_yi.clone() + xi_yi - xi.clone() - yi.clone() + F::one();
    }
    end_timer!(start);
    Ok(res)
}

/// This function build the eq(x, r) polynomial for any given r.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r<F: Field>(
    r: &[F],
    config: F::R,
) -> Result<DenseMultilinearExtension<F>, ArithErrors> {
    let evals = build_eq_x_r_vec(r)?;
    let mle = DenseMultilinearExtension::from_evaluations_vec(r.len(), evals, config);

    Ok(mle)
}
/// This function build the eq(x, r) polynomial for any given r, and output the
/// evaluation of eq(x, r) in its vector form.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r_vec<F: Field>(r: &[F]) -> Result<Vec<F>, ArithErrors> {
    // we build eq(x,r) from its evaluations
    // we want to evaluate eq(x,r) over x \in {0, 1}^num_vars
    // for example, with num_vars = 4, x is a binary vector of 4, then
    //  0 0 0 0 -> (1-r0)   * (1-r1)    * (1-r2)    * (1-r3)
    //  1 0 0 0 -> r0       * (1-r1)    * (1-r2)    * (1-r3)
    //  0 1 0 0 -> (1-r0)   * r1        * (1-r2)    * (1-r3)
    //  1 1 0 0 -> r0       * r1        * (1-r2)    * (1-r3)
    //  ....
    //  1 1 1 1 -> r0       * r1        * r2        * r3
    // we will need 2^num_var evaluations

    let mut eval = Vec::new();
    build_eq_x_r_helper(r, &mut eval)?;

    Ok(eval)
}

/// A helper function to build eq(x, r) recursively.
/// This function takes `r.len()` steps, and for each step it requires a maximum
/// `r.len()-1` multiplications.
fn build_eq_x_r_helper<F: Field>(r: &[F], buf: &mut Vec<F>) -> Result<(), ArithErrors> {
    if r.is_empty() {
        return Err(ArithErrors::InvalidParameters("r length is 0".into()));
    } else if r.len() == 1 {
        // initializing the buffer with [1-r_0, r_0]
        buf.push(F::one() - r[0].clone());
        buf.push(r[0].clone());
    } else {
        build_eq_x_r_helper(&r[1..], buf)?;

        // suppose at the previous step we received [b_1, ..., b_k]
        // for the current step we will need
        // if x_0 = 0:   (1-r0) * [b_1, ..., b_k]
        // if x_0 = 1:   r0 * [b_1, ..., b_k]
        // let mut res = vec![];
        // for &b_i in buf.iter() {
        //     let tmp = r[0] * b_i;
        //     res.push(b_i - tmp);
        //     res.push(tmp);
        // }
        // *buf = res;

        let mut res = vec![F::zero(); buf.len() << 1];
        cfg_iter_mut!(res).enumerate().for_each(|(i, val)| {
            let bi = buf[i >> 1].clone();
            let tmp = r[0].clone() * bi.clone();
            if (i & 1) == 0 {
                *val = bi - tmp;
            } else {
                *val = tmp;
            }
        });
        *buf = res;
    }

    Ok(())
}
