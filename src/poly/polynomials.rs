// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// Adapted for rings by Nethermind

use ark_std::{end_timer, start_timer};
use rand_core::RngCore;

use crate::{
    poly::{
        ArithErrors, RefCounter, dense::DenseMultilinearExtension, mle::MultilinearExtension,
        utils::get_batched_nv,
    },
    traits::Ring,
};

/// Sample a random list of multilinear polynomials.
/// Returns
/// - the list of polynomials,
/// - its sum of polynomial evaluations over the boolean hypercube.
pub fn random_mle_list<F: Ring, Rn: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut Rn,
) -> (Vec<RefCounter<DenseMultilinearExtension<F>>>, F) {
    let start = start_timer!(|| "sample random mle list");
    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv));
    }
    let mut sum = F::zero();

    for _ in 0..1 << nv {
        let mut product = F::one();

        for e in multiplicands.iter_mut() {
            let val = F::random(rng);
            e.push(val);
            product = product * val;
        }
        sum += &product;
    }

    let list = multiplicands
        .into_iter()
        .map(|x| RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    (list, sum)
}

// Build a randomize list of mle-s whose sum is zero.
pub fn random_zero_mle_list<F: Ring, Rn: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut Rn,
) -> Vec<RefCounter<DenseMultilinearExtension<F>>> {
    let start = start_timer!(|| "sample random zero mle list");

    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv));
    }
    for _ in 0..1 << nv {
        multiplicands[0].push(F::zero());
        for e in multiplicands.iter_mut().skip(1) {
            e.push(F::random(rng));
        }
    }

    let list = multiplicands
        .into_iter()
        .map(|x| RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    list
}

pub fn identity_permutation<F: Ring + From<u64>>(num_vars: usize, num_chunks: usize) -> Vec<F> {
    let len = (num_chunks as u64) * (1u64 << num_vars);
    (0..len).map(F::from).collect()
}

/// A list of MLEs that represents an identity permutation
pub fn identity_permutation_mles<F: Ring + From<u64>>(
    num_vars: usize,
    num_chunks: usize,
) -> Vec<RefCounter<DenseMultilinearExtension<F>>> {
    let mut res = vec![];
    for i in 0..num_chunks {
        let shift = (i * (1 << num_vars)) as u64;
        let s_id_vec = (shift..shift + (1u64 << num_vars)).map(F::from).collect();
        res.push(RefCounter::new(
            DenseMultilinearExtension::from_evaluations_vec(num_vars, s_id_vec),
        ));
    }
    res
}

pub fn random_permutation<F: Ring + From<u64>, Rn: RngCore>(
    num_vars: usize,
    num_chunks: usize,
    rng: &mut Rn,
) -> Vec<F> {
    let len = (num_chunks as u64) * (1u64 << num_vars);
    let mut s_id_vec: Vec<F> = (0..len).map(F::from).collect();
    let mut s_perm_vec = vec![];
    for _ in 0..len {
        let index = (rng.next_u64() as usize) % s_id_vec.len();
        s_perm_vec.push(s_id_vec.remove(index));
    }
    s_perm_vec
}

/// A list of MLEs that represent a random permutation
pub fn random_permutation_mles<F: Ring + From<u64>, Rn: RngCore>(
    num_vars: usize,
    num_chunks: usize,
    rng: &mut Rn,
) -> Vec<RefCounter<DenseMultilinearExtension<F>>> {
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

pub fn evaluate_opt<F: Ring>(poly: &DenseMultilinearExtension<F>, point: &[F]) -> F {
    assert_eq!(poly.num_vars, point.len());
    fix_variables(poly, point).evaluations[0]
}

pub fn fix_variables<F: Ring>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &[F],
) -> DenseMultilinearExtension<F> {
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

    DenseMultilinearExtension::<F>::from_evaluations_slice(nv - dim, &poly[..1 << (nv - dim)])
}

fn fix_one_variable_helper<F: Ring>(data: &[F], nv: usize, point: &F) -> Vec<F> {
    let mut res = vec![F::zero(); 1 << (nv - 1)];

    // evaluate single variable of partial point from left to right

    for i in 0..1 << (nv - 1) {
        res[i] = data[i] + (data[(i << 1) + 1] - data[i << 1]) * point;
    }

    res
}

pub fn evaluate_no_par<F: Ring>(poly: &DenseMultilinearExtension<F>, point: &[F]) -> F {
    assert_eq!(poly.num_vars, point.len());
    fix_variables_no_par(poly, point).evaluations[0]
}

fn fix_variables_no_par<F: Ring>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &[F],
) -> DenseMultilinearExtension<F> {
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
pub fn merge_polynomials<F: Ring>(
    polynomials: &[RefCounter<DenseMultilinearExtension<F>>],
) -> Result<RefCounter<DenseMultilinearExtension<F>>, ArithErrors> {
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
    scalars.extend_from_slice(vec![F::zero(); (1 << merged_nv) - scalars.len()].as_ref());
    Ok(RefCounter::new(
        DenseMultilinearExtension::from_evaluations_vec(merged_nv, scalars),
    ))
}

pub fn fix_last_variables_no_par<F: Ring>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &[F],
) -> DenseMultilinearExtension<F> {
    let mut res = fix_last_variable_no_par(poly, partial_point.last().unwrap());
    for p in partial_point.iter().rev().skip(1) {
        res = fix_last_variable_no_par(&res, p);
    }
    res
}

fn fix_last_variable_no_par<F: Ring>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &F,
) -> DenseMultilinearExtension<F> {
    let nv = poly.num_vars();
    let half_len = 1 << (nv - 1);
    let mut res = vec![F::zero(); half_len];
    for (i, e) in res.iter_mut().enumerate().take(half_len) {
        *e = poly.evaluations[i]
            + *partial_point * (poly.evaluations[i + half_len] - poly.evaluations[i]);
    }
    DenseMultilinearExtension::from_evaluations_vec(nv - 1, res)
}

pub fn fix_last_variables<F: Ring>(
    poly: &DenseMultilinearExtension<F>,
    partial_point: &[F],
) -> DenseMultilinearExtension<F> {
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

    DenseMultilinearExtension::from_evaluations_slice(nv - dim, &poly[..1 << (nv - dim)])
}

fn fix_last_variable_helper<F: Ring>(data: &[F], nv: usize, point: &F) -> Vec<F> {
    let half_len = 1 << (nv - 1);
    let mut res = vec![F::zero(); half_len];

    // evaluate single variable of partial point from left to right

    for b in 0..half_len {
        res[b] = data[b] + (data[b + half_len] - data[b]) * point;
    }

    res
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{U128, const_monty_params};
    use num_traits::{One, Zero};

    use super::*;
    use crate::{field::RandomField, poly::mle::MultilinearExtension};

    const_monty_params!(ModP, U128, "0076F668F4274572E39A3EA8285319B5");
    type Fp = RandomField<ModP, { U128::LIMBS }>;

    #[test]
    fn test_identity_permutation_values_and_mles() {
        let v = identity_permutation::<Fp>(3, 2);
        let expected: Vec<Fp> = (0u64..(2 * (1u64 << 3))).map(Fp::from).collect();
        assert_eq!(v, expected);

        let mles = identity_permutation_mles::<Fp>(2, 3);
        assert_eq!(mles.len(), 3);
        for (i, m) in mles.iter().enumerate() {
            assert_eq!(m.num_vars(), 2);
            let shift = (i * (1 << 2)) as u64;
            let evals = m.to_evaluations();
            let expected_chunk: Vec<Fp> = (shift..shift + (1u64 << 2)).map(Fp::from).collect();
            assert_eq!(evals, expected_chunk);
        }
    }

    #[test]
    fn test_evaluate_vs_trait_evaluate() {
        let evals = vec![
            Fp::from(1u64),
            Fp::from(2u64),
            Fp::from(3u64),
            Fp::from(4u64),
        ];
        let mle = DenseMultilinearExtension::from_evaluations_vec(2, evals);
        for p in [
            vec![Fp::zero(), Fp::zero()],
            vec![Fp::one(), Fp::zero()],
            vec![Fp::zero(), Fp::one()],
            vec![Fp::one(), Fp::one()],
        ] {
            let b = evaluate_no_par(&mle, &p);
            let c = mle.evaluate(&p).unwrap();
            assert_eq!(b, c);
        }
    }

    #[test]
    fn test_fix_variables_and_fix_last_variables_match_no_par() {
        let evals: Vec<Fp> = (0u64..8u64).map(Fp::from).collect();
        let mle = DenseMultilinearExtension::from_evaluations_vec(3, evals);

        let r = [Fp::from(7u64), Fp::from(11u64)];
        let by_trait = mle.fixed_variables(&r);
        let slow = fix_variables_no_par(&mle, &r);
        assert_eq!(by_trait.num_vars(), slow.num_vars());
        assert_eq!(by_trait.to_evaluations(), slow.to_evaluations());

        let r2 = [Fp::from(5u64), Fp::from(3u64)];
        let reduced = fix_last_variables_no_par(&mle, &r2);
        assert_eq!(reduced.num_vars(), 1);
        for s in [Fp::zero(), Fp::one()] {
            let lhs = reduced.evaluate(&[s]).unwrap();
            let rhs = mle.evaluate(&[s, r2[0], r2[1]]).unwrap();
            assert_eq!(lhs, rhs);
        }
    }

    #[test]
    fn test_merge_polynomials_padding_and_error() {
        let p1 = RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![Fp::from(1), Fp::from(2), Fp::from(3), Fp::from(4)],
        ));
        let p2 = RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![Fp::from(5), Fp::from(6), Fp::from(7), Fp::from(8)],
        ));
        let p3 = RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![Fp::from(9), Fp::from(10), Fp::from(11), Fp::from(12)],
        ));

        let merged = merge_polynomials(&[p1.clone(), p2.clone(), p3.clone()]).unwrap();
        // merged_nv should be enough to fit 3 * 4 = 12 elements => next power of two: 16 -> nv = 4
        assert_eq!(merged.num_vars(), 4);
        let mut expected = vec![];
        expected.extend_from_slice(&p1.to_evaluations());
        expected.extend_from_slice(&p2.to_evaluations());
        expected.extend_from_slice(&p3.to_evaluations());
        expected.resize(1 << merged.num_vars(), Fp::zero());
        assert_eq!(merged.to_evaluations(), expected);

        let bad = RefCounter::new(DenseMultilinearExtension::from_evaluations_vec(
            1,
            vec![Fp::from(13), Fp::from(17)],
        ));
        let err = merge_polynomials(&[p1, bad, p3]).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("num_vars do not match"));
    }

    #[test]
    fn test_random_zero_mle_list_sum_of_products_is_zero() {
        let nv = 3usize;
        let degree = 3usize;
        struct Dummy(u64);
        impl RngCore for Dummy {
            fn next_u32(&mut self) -> u32 {
                let x = self.0 as u32;
                self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
                x
            }
            fn next_u64(&mut self) -> u64 {
                let x = self.0;
                self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
                x
            }
            fn fill_bytes(&mut self, dest: &mut [u8]) {
                for chunk in dest.chunks_mut(8) {
                    let r = self.next_u64().to_le_bytes();
                    let len = chunk.len();
                    chunk.copy_from_slice(&r[..len]);
                }
            }
        }
        let mut rng = Dummy(123456789);
        let list = random_zero_mle_list::<Fp, _>(nv, degree, &mut rng);
        assert_eq!(list.len(), degree);
        let n = 1usize << nv;
        assert!(list[0].to_evaluations().iter().all(|v| v.is_zero()));
        let mut total = Fp::zero();
        for idx in 0..n {
            let mut prod = Fp::one();
            for poly in list.iter() {
                prod *= &poly[idx];
            }
            total += &prod;
        }
        assert!(total.is_zero());
    }
}
