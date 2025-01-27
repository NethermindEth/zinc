// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// Adapted for rings by Nethermind

use ark_std::{end_timer, rand::RngCore, start_timer, string::ToString, vec::*};
#[cfg(feature = "parallel")]
use rayon::prelude::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

use super::{util::get_batched_nv, ArithErrors, RefCounter};
pub use crate::mle::DenseMultilinearExtension;
use crate::mle::MultilinearExtension;
use lattirust_ring::Ring;

/// Sample a random list of multilinear polynomials.
/// Returns
/// - the list of polynomials,
/// - its sum of polynomial evaluations over the boolean hypercube.
pub fn random_mle_list<R: Ring, Rn: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut Rn,
) -> (Vec<RefCounter<DenseMultilinearExtension<R>>>, R) {
    let start = start_timer!(|| "sample random mle list");
    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv));
    }
    let mut sum = R::zero();

    for _ in 0..1 << nv {
        let mut product = R::one();

        for e in multiplicands.iter_mut() {
            let val = R::rand(rng);
            e.push(val);
            product *= val;
        }
        sum += product;
    }

    let list = multiplicands
        .into_iter()
        .map(|x| RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    (list, sum)
}

// Build a randomize list of mle-s whose sum is zero.
pub fn random_zero_mle_list<R: Ring, Rn: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut Rn,
) -> Vec<RefCounter<DenseMultilinearExtension<R>>> {
    let start = start_timer!(|| "sample random zero mle list");

    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv));
    }
    for _ in 0..1 << nv {
        multiplicands[0].push(R::zero());
        for e in multiplicands.iter_mut().skip(1) {
            e.push(R::rand(rng));
        }
    }

    let list = multiplicands
        .into_iter()
        .map(|x| RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    list
}

pub fn identity_permutation<R: Ring>(num_vars: usize, num_chunks: usize) -> Vec<R> {
    let len = (num_chunks as u64) * (1u64 << num_vars);
    (0..len).map(R::from).collect()
}

/// A list of MLEs that represents an identity permutation
pub fn identity_permutation_mles<R: Ring>(
    num_vars: usize,
    num_chunks: usize,
) -> Vec<RefCounter<DenseMultilinearExtension<R>>> {
    let mut res = vec![];
    for i in 0..num_chunks {
        let shift = (i * (1 << num_vars)) as u64;
        let s_id_vec = (shift..shift + (1u64 << num_vars)).map(R::from).collect();
        res.push(RefCounter::new(
            DenseMultilinearExtension::from_evaluations_vec(num_vars, s_id_vec),
        ));
    }
    res
}

pub fn random_permutation<R: Ring, Rn: RngCore>(
    num_vars: usize,
    num_chunks: usize,
    rng: &mut Rn,
) -> Vec<R> {
    let len = (num_chunks as u64) * (1u64 << num_vars);
    let mut s_id_vec: Vec<R> = (0..len).map(R::from).collect();
    let mut s_perm_vec = vec![];
    for _ in 0..len {
        let index = (rng.next_u64() as usize) % s_id_vec.len();
        s_perm_vec.push(s_id_vec.remove(index));
    }
    s_perm_vec
}

/// A list of MLEs that represent a random permutation
pub fn random_permutation_mles<R: Ring, Rn: RngCore>(
    num_vars: usize,
    num_chunks: usize,
    rng: &mut Rn,
) -> Vec<RefCounter<DenseMultilinearExtension<R>>> {
    let s_perm_vec = random_permutation(num_vars, num_chunks, rng);
    let mut res = vec![];
    let n = 1 << num_vars;
    for i in 0..num_chunks {
        res.push(RefCounter::new(
            DenseMultilinearExtension::from_evaluations_vec(
                num_vars,
                s_perm_vec[i * n..i * n + n].to_vec(),
            ),
        ));
    }
    res
}

pub fn evaluate_opt<R: Ring>(poly: &DenseMultilinearExtension<R>, point: &[R]) -> R {
    assert_eq!(poly.num_vars, point.len());
    fix_variables(poly, point).evaluations[0]
}

pub fn fix_variables<R: Ring>(
    poly: &DenseMultilinearExtension<R>,
    partial_point: &[R],
) -> DenseMultilinearExtension<R> {
    assert!(
        partial_point.len() <= poly.num_vars,
        "invalid size of partial point"
    );
    let nv = poly.num_vars;
    let mut poly = poly.evaluations.to_vec();
    let dim = partial_point.len();
    // evaluate single variable of partial point from left to right
    for (i, point) in partial_point.iter().enumerate().take(dim) {
        poly = fix_one_variable_helper(&poly, nv - i, point);
    }

    DenseMultilinearExtension::<R>::from_evaluations_slice(nv - dim, &poly[..1 << (nv - dim)])
}

fn fix_one_variable_helper<R: Ring>(data: &[R], nv: usize, point: &R) -> Vec<R> {
    let mut res = vec![R::zero(); 1 << (nv - 1)];

    // evaluate single variable of partial point from left to right
    #[cfg(not(feature = "parallel"))]
    for i in 0..1 << (nv - 1) {
        res[i] = data[i] + (data[(i << 1) + 1] - data[i << 1]) * point;
    }

    #[cfg(feature = "parallel")]
    res.par_iter_mut().enumerate().for_each(|(i, x)| {
        *x = data[i << 1] + (data[(i << 1) + 1] - data[i << 1]) * point;
    });

    res
}

pub fn evaluate_no_par<R: Ring>(poly: &DenseMultilinearExtension<R>, point: &[R]) -> R {
    assert_eq!(poly.num_vars, point.len());
    fix_variables_no_par(poly, point).evaluations[0]
}

fn fix_variables_no_par<R: Ring>(
    poly: &DenseMultilinearExtension<R>,
    partial_point: &[R],
) -> DenseMultilinearExtension<R> {
    assert!(
        partial_point.len() <= poly.num_vars,
        "invalid size of partial point"
    );
    let nv = poly.num_vars;
    let mut poly = poly.evaluations.to_vec();
    let dim = partial_point.len();
    // evaluate single variable of partial point from left to right
    for i in 1..dim + 1 {
        let r = partial_point[i - 1];
        for b in 0..1 << (nv - i) {
            poly[b] = poly[b << 1] + (poly[(b << 1) + 1] - poly[b << 1]) * r;
        }
    }
    DenseMultilinearExtension::from_evaluations_slice(nv - dim, &poly[..1 << (nv - dim)])
}

/// merge a set of polynomials. Returns an error if the
/// polynomials do not share a same number of nvs.
pub fn merge_polynomials<R: Ring>(
    polynomials: &[RefCounter<DenseMultilinearExtension<R>>],
) -> Result<RefCounter<DenseMultilinearExtension<R>>, ArithErrors> {
    let nv = polynomials[0].num_vars();
    for poly in polynomials.iter() {
        if nv != poly.num_vars() {
            return Err(ArithErrors::InvalidParameters(
                "num_vars do not match for polynomials".to_string(),
            ));
        }
    }

    let merged_nv = get_batched_nv(nv, polynomials.len());
    let mut scalars = vec![];
    for poly in polynomials.iter() {
        scalars.extend_from_slice(poly.to_evaluations().as_slice());
    }
    scalars.extend_from_slice(vec![R::zero(); (1 << merged_nv) - scalars.len()].as_ref());
    Ok(RefCounter::new(
        DenseMultilinearExtension::from_evaluations_vec(merged_nv, scalars),
    ))
}

pub fn fix_last_variables_no_par<R: Ring>(
    poly: &DenseMultilinearExtension<R>,
    partial_point: &[R],
) -> DenseMultilinearExtension<R> {
    let mut res = fix_last_variable_no_par(poly, partial_point.last().unwrap());
    for p in partial_point.iter().rev().skip(1) {
        res = fix_last_variable_no_par(&res, p);
    }
    res
}

fn fix_last_variable_no_par<R: Ring>(
    poly: &DenseMultilinearExtension<R>,
    partial_point: &R,
) -> DenseMultilinearExtension<R> {
    let nv = poly.num_vars();
    let half_len = 1 << (nv - 1);
    let mut res = vec![R::zero(); half_len];
    for (i, e) in res.iter_mut().enumerate().take(half_len) {
        *e = poly.evaluations[i]
            + *partial_point * (poly.evaluations[i + half_len] - poly.evaluations[i]);
    }
    DenseMultilinearExtension::from_evaluations_vec(nv - 1, res)
}
pub fn fix_last_variables<R: Ring>(
    poly: &DenseMultilinearExtension<R>,
    partial_point: &[R],
) -> DenseMultilinearExtension<R> {
    assert!(
        partial_point.len() <= poly.num_vars,
        "invalid size of partial point"
    );
    let nv = poly.num_vars;
    let mut poly = poly.evaluations.to_vec();
    let dim = partial_point.len();
    // evaluate single variable of partial point from left to right
    for (i, point) in partial_point.iter().rev().enumerate().take(dim) {
        poly = fix_last_variable_helper(&poly, nv - i, point);
    }

    DenseMultilinearExtension::<R>::from_evaluations_slice(nv - dim, &poly[..1 << (nv - dim)])
}

fn fix_last_variable_helper<R: Ring>(data: &[R], nv: usize, point: &R) -> Vec<R> {
    let half_len = 1 << (nv - 1);
    let mut res = vec![R::zero(); half_len];

    // evaluate single variable of partial point from left to right
    #[cfg(not(feature = "parallel"))]
    for b in 0..half_len {
        res[b] = data[b] + (data[b + half_len] - data[b]) * point;
    }

    #[cfg(feature = "parallel")]
    res.par_iter_mut().enumerate().for_each(|(i, x)| {
        *x = data[i] + (data[i + half_len] - data[i]) * point;
    });

    res
}
