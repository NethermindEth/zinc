// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// Adapted for rings by Nethermind

use ark_std::{end_timer, rand::RngCore, start_timer, vec, vec::Vec};
use crypto_bigint::Random;
use num_traits::{One, Zero};

use super::RefCounter;
pub use crate::poly_f::mle::DenseMultilinearExtension;
use crate::{
    field::RandomField,
    poly::{ArithErrors, get_batched_nv},
    poly_f::mle::MultilinearExtension,
    traits::{ConfigReference, FieldMap},
};

/// Sample a random list of multilinear polynomials.
/// Returns
/// - the list of polynomials,
/// - its sum of polynomial evaluations over the boolean hypercube.
pub fn random_mle_list<C: ConfigReference, Rn: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut Rn,
    config: C,
) -> (
    Vec<RefCounter<DenseMultilinearExtension<C>>>,
    RandomField<C>,
) {
    let start = start_timer!(|| "sample random mle list");
    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv));
    }
    let mut sum = RandomField::zero();

    for _ in 0..1 << nv {
        let mut product = RandomField::one();

        for e in multiplicands.iter_mut() {
            let val = RandomField::rand_with_config(rng, config);
            e.push(val.clone());
            product *= &val;
        }
        sum += &product;
    }

    let list = multiplicands
        .into_iter()
        .map(|x| {
            RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(
                nv,
                x,
                Some(config),
            ))
        })
        .collect();

    end_timer!(start);
    (list, sum)
}

// Build a randomize list of mle-s whose sum is zero.
pub fn random_zero_mle_list<C: ConfigReference, Rn: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut Rn,
    config: C,
) -> Vec<RefCounter<DenseMultilinearExtension<C>>> {
    let start = start_timer!(|| "sample random zero mle list");

    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv));
    }
    for _ in 0..1 << nv {
        multiplicands[0].push(RandomField::zero());
        for e in multiplicands.iter_mut().skip(1) {
            e.push(RandomField::random(rng));
        }
    }

    let list = multiplicands
        .into_iter()
        .map(|x| {
            RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(
                nv,
                x,
                Some(config),
            ))
        })
        .collect();

    end_timer!(start);
    list
}

pub fn identity_permutation<C: ConfigReference>(
    num_vars: usize,
    num_chunks: usize,
    config: C,
) -> Vec<RandomField<C>> {
    let len = (num_chunks as u64) * (1u64 << num_vars);
    (0..len).map(|i| i.map_to_field(config)).collect()
}

/// A list of MLEs that represents an identity permutation
pub fn identity_permutation_mles<C: ConfigReference>(
    num_vars: usize,
    num_chunks: usize,
    config: C,
) -> Vec<RefCounter<DenseMultilinearExtension<C>>> {
    let mut res = vec![];
    for i in 0..num_chunks {
        let shift = (i * (1 << num_vars)) as u64;
        let s_id_vec = (shift..shift + (1u64 << num_vars))
            .map(|i| i.map_to_field(config))
            .collect();
        res.push(RefCounter::new(
            DenseMultilinearExtension::from_evaluations_vec(num_vars, s_id_vec, Some(config)),
        ));
    }
    res
}

pub fn random_permutation<C: ConfigReference, Rn: RngCore>(
    num_vars: usize,
    num_chunks: usize,
    rng: &mut Rn,
    config: C,
) -> Vec<RandomField<C>> {
    let len = (num_chunks as u64) * (1u64 << num_vars);
    let mut s_id_vec: Vec<RandomField<C>> = (0..len).map(|i| i.map_to_field(config)).collect();
    let mut s_perm_vec = vec![];
    for _ in 0..len {
        let index = (rng.next_u64() as usize) % s_id_vec.len();
        s_perm_vec.push(s_id_vec.remove(index));
    }
    s_perm_vec
}

/// A list of MLEs that represent a random permutation
pub fn random_permutation_mles<C: ConfigReference, Rn: RngCore>(
    num_vars: usize,
    num_chunks: usize,
    rng: &mut Rn,
    config: C,
) -> Vec<RefCounter<DenseMultilinearExtension<C>>> {
    let s_perm_vec = random_permutation(num_vars, num_chunks, rng, config);
    let mut res = vec![];
    let n = 1 << num_vars;
    for i in 0..num_chunks {
        res.push(RefCounter::new(
            DenseMultilinearExtension::from_evaluations_vec(
                num_vars,
                s_perm_vec[i * n..i * n + n].to_vec(),
                Some(config),
            ),
        ));
    }
    res
}

pub fn evaluate_opt<C: ConfigReference>(
    poly: &DenseMultilinearExtension<C>,
    point: &[RandomField<C>],
    config: C,
) -> RandomField<C> {
    assert_eq!(poly.num_vars, point.len());
    fix_variables(poly, point, config).evaluations[0].clone()
}

pub fn fix_variables<C: ConfigReference>(
    poly: &DenseMultilinearExtension<C>,
    partial_point: &[RandomField<C>],
    config: C,
) -> DenseMultilinearExtension<C> {
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

    DenseMultilinearExtension::<C>::from_evaluations_slice(
        nv - dim,
        &poly[..1 << (nv - dim)],
        config,
    )
}

fn fix_one_variable_helper<C: ConfigReference>(
    data: &[RandomField<C>],
    nv: usize,
    point: &RandomField<C>,
) -> Vec<RandomField<C>> {
    let mut res = vec![RandomField::zero(); 1 << (nv - 1)];

    // evaluate single variable of partial point from left to right

    for i in 0..1 << (nv - 1) {
        res[i] = data[i].clone() + (data[(i << 1) + 1].clone() - data[i << 1].clone()) * point;
    }

    res
}

pub fn evaluate_no_par<C: ConfigReference>(
    poly: &DenseMultilinearExtension<C>,
    point: &[RandomField<C>],
    config: C,
) -> RandomField<C> {
    assert_eq!(poly.num_vars, point.len());
    fix_variables_no_par(poly, point, config).evaluations[0].clone()
}

fn fix_variables_no_par<C: ConfigReference>(
    poly: &DenseMultilinearExtension<C>,
    partial_point: &[RandomField<C>],
    config: C,
) -> DenseMultilinearExtension<C> {
    assert!(
        partial_point.len() <= poly.num_vars,
        "invalid size of partial point"
    );
    let nv = poly.num_vars;
    let mut poly = poly.evaluations.to_vec();
    let dim = partial_point.len();
    // evaluate single variable of partial point from left to right
    for i in 1..dim + 1 {
        let r = partial_point[i - 1].clone();
        for b in 0..1 << (nv - i) {
            poly[b] =
                poly[b << 1].clone() + (poly[(b << 1) + 1].clone() - poly[b << 1].clone()) * &r;
        }
    }
    DenseMultilinearExtension::from_evaluations_slice(nv - dim, &poly[..1 << (nv - dim)], config)
}

/// merge a set of polynomials. Returns an error if the
/// polynomials do not share a same number of nvs.
pub fn merge_polynomials<C: ConfigReference>(
    polynomials: &[RefCounter<DenseMultilinearExtension<C>>],
    config: C,
) -> Result<RefCounter<DenseMultilinearExtension<C>>, ArithErrors> {
    let nv = polynomials[0].num_vars();
    for poly in polynomials.iter() {
        if nv != poly.num_vars() {
            return Err(ArithErrors::InvalidParameters(
                "num_vars do not match for polynomials".into(),
            ));
        }
    }

    let merged_nv = get_batched_nv(nv, polynomials.len());
    let mut scalars = vec![];
    for poly in polynomials.iter() {
        scalars.extend_from_slice(poly.to_evaluations().as_slice());
    }
    scalars.extend_from_slice(vec![RandomField::zero(); (1 << merged_nv) - scalars.len()].as_ref());
    Ok(RefCounter::new(
        DenseMultilinearExtension::from_evaluations_vec(merged_nv, scalars, Some(config)),
    ))
}

pub fn fix_last_variables_no_par<C: ConfigReference>(
    poly: &DenseMultilinearExtension<C>,
    partial_point: &[RandomField<C>],
    config: C,
) -> DenseMultilinearExtension<C> {
    let mut res = fix_last_variable_no_par(poly, partial_point.last().unwrap(), config);
    for p in partial_point.iter().rev().skip(1) {
        res = fix_last_variable_no_par(&res, p, config);
    }
    res
}

fn fix_last_variable_no_par<C: ConfigReference>(
    poly: &DenseMultilinearExtension<C>,
    partial_point: &RandomField<C>,
    config: C,
) -> DenseMultilinearExtension<C> {
    let nv = poly.num_vars();
    let half_len = 1 << (nv - 1);
    let mut res = vec![RandomField::zero(); half_len];
    for (i, e) in res.iter_mut().enumerate().take(half_len) {
        *e = poly.evaluations[i].clone()
            + partial_point.clone()
                * (poly.evaluations[i + half_len].clone() - poly.evaluations[i].clone());
    }
    DenseMultilinearExtension::from_evaluations_vec(nv - 1, res, Some(config))
}
pub fn fix_last_variables<C: ConfigReference>(
    poly: &DenseMultilinearExtension<C>,
    partial_point: &[RandomField<C>],
    config: C,
) -> DenseMultilinearExtension<C> {
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

    DenseMultilinearExtension::from_evaluations_slice(nv - dim, &poly[..1 << (nv - dim)], config)
}

fn fix_last_variable_helper<C: ConfigReference>(
    data: &[RandomField<C>],
    nv: usize,
    point: &RandomField<C>,
) -> Vec<RandomField<C>> {
    let half_len = 1 << (nv - 1);
    let mut res = vec![RandomField::zero(); half_len];

    // evaluate single variable of partial point from left to right

    for b in 0..half_len {
        res[b] = data[b].clone() + (data[b + half_len].clone() - data[b].clone()) * point;
    }

    res
}
