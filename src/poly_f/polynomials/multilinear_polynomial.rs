// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// Adapted for rings by Nethermind

use ark_ff::{One, Zero};
use ark_std::{end_timer, rand::RngCore, start_timer, vec::*};

use super::RefCounter;
use crate::poly_f::mle::dense::DenseMultilinearExtension;
use crate::{
    field::{rand_with_config, RandomField},
    field_config::FieldConfig,
};

/// Sample a random list of multilinear polynomials.
/// Returns
/// - the list of polynomials,
/// - its sum of polynomial evaluations over the boolean hypercube.
pub fn random_mle_list<const N: usize, Rn: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut Rn,
    config: *const FieldConfig<N>,
) -> (
    Vec<RefCounter<DenseMultilinearExtension<N>>>,
    RandomField<N>,
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
            let val = rand_with_config(rng, config);
            e.push(val);
            product *= val;
        }
        sum += &product;
    }

    let list = multiplicands
        .into_iter()
        .map(|x| {
            RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(
                nv, x, config,
            ))
        })
        .collect();

    end_timer!(start);
    (list, sum)
}
