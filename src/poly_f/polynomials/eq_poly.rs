use ark_ff::One;

use crate::{field::RandomField as F, lookup::Math};

use super::util::unsafe_allocate_zero_vec;

pub struct EqPolynomial<const N: usize> {
    r: Vec<F<N>>,
}

pub struct EqPlusOnePolynomial<const N: usize> {
    x: Vec<F<N>>,
}

const PARALLEL_THRESHOLD: usize = 16;

impl<const N: usize> EqPolynomial<N> {
    pub fn new(r: Vec<F<N>>) -> Self {
        EqPolynomial { r }
    }

    pub fn evaluate(&self, rx: &[F<N>]) -> F<N> {
        assert_eq!(self.r.len(), rx.len());
        (0..rx.len())
            .map(|i| self.r[i] * rx[i] + (F::one() - self.r[i]) * (F::one() - rx[i]))
            .product()
    }

    pub fn evals(r: &[F<N>]) -> Vec<F<N>> {
        match r.len() {
            0..=PARALLEL_THRESHOLD => Self::evals_serial(r, None),
            _ => Self::evals_parallel(r, None),
        }
    }

    /// Computes the table of coefficients:
    ///     scaling_factor * eq(r, x) for all x in {0, 1}^n
    /// serially. More efficient for short `r`.
    fn evals_serial(r: &[F<N>], scaling_factor: Option<F<N>>) -> Vec<F<N>> {
        let mut evals: Vec<F<N>> = vec![scaling_factor.unwrap_or(F::one()); r.len().pow2()];
        let mut size = 1;
        for j in 0..r.len() {
            // in each iteration, we double the size of chis
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                // copy each element from the prior iteration twice
                let scalar = evals[i / 2];
                evals[i] = scalar * r[j];
                evals[i - 1] = scalar - evals[i];
            }
        }
        evals
    }

    /// Computes the table of coefficients:
    ///     scaling_factor * eq(r, x) for all x in {0, 1}^n
    /// computing biggest layers of the dynamic programming tree in parallel.

    pub fn evals_parallel(r: &[F<N>], scaling_factor: Option<F<N>>) -> Vec<F<N>> {
        let final_size = r.len().pow2();
        let mut evals: Vec<F<N>> = unsafe_allocate_zero_vec(final_size);
        let mut size = 1;
        evals[0] = scaling_factor.unwrap_or(F::one());

        for r in r.iter().rev() {
            let (evals_left, evals_right) = evals.split_at_mut(size);
            let (evals_right, _) = evals_right.split_at_mut(size);

            evals_left
                .iter_mut()
                .zip(evals_right.iter_mut())
                .for_each(|(x, y)| {
                    *y = *x * *r;
                    *x -= *y;
                });

            size *= 2;
        }

        evals
    }
}

impl<const N: usize> EqPlusOnePolynomial<N> {
    pub fn new(x: Vec<F<N>>) -> Self {
        EqPlusOnePolynomial { x }
    }

    /* This MLE is 1 if y = x + 1 for x in the range [0... 2^l-2].
    That is, it ignores the case where x is all 1s, outputting 0.
    Assumes x and y are provided big-endian. */
    pub fn evaluate(&self, y: &[F<N>]) -> F<N> {
        let l = self.x.len();
        let x = &self.x;
        assert!(y.len() == l);
        let one = F::one();

        /* If y+1 = x, then the two bit vectors are of the following form.
            Let k be the longest suffix of 1s in x.
            In y, those k bits are 0.
            Then, the next bit in x is 0 and the next bit in y is 1.
            The remaining higher bits are the same in x and y.
        */
        (0..l)
            .into_iter()
            .map(|k| {
                let lower_bits_product = (0..k)
                    .map(|i| x[l - 1 - i] * (F::one() - y[l - 1 - i]))
                    .product::<F<N>>();
                let kth_bit_product = (F::one() - x[l - 1 - k]) * y[l - 1 - k];
                let higher_bits_product = ((k + 1)..l)
                    .map(|i| {
                        x[l - 1 - i] * y[l - 1 - i] + (one - x[l - 1 - i]) * (one - y[l - 1 - i])
                    })
                    .product::<F<N>>();
                lower_bits_product * kth_bit_product * higher_bits_product
            })
            .sum()
    }

    pub fn evals(r: &[F<N>], scaling_factor: Option<F<N>>) -> (Vec<F<N>>, Vec<F<N>>) {
        let ell = r.len();
        let mut eq_evals: Vec<F<N>> = unsafe_allocate_zero_vec(ell.pow2());
        eq_evals[0] = scaling_factor.unwrap_or(F::one());
        let mut eq_plus_one_evals: Vec<F<N>> = unsafe_allocate_zero_vec(ell.pow2());

        // i indicates the LENGTH of the prefix of r for which the eq_table is calculated
        let eq_evals_helper = |eq_evals: &mut Vec<F<N>>, r: &[F<N>], i: usize| {
            debug_assert!(i != 0);
            let step = 1 << (ell - i); // step = (full / size)/2

            let mut selected: Vec<_> = eq_evals.iter_mut().step_by(step).collect();

            selected.chunks_mut(2).for_each(|chunk| {
                *chunk[1] = *chunk[0] * r[i - 1];
                *chunk[0] -= *chunk[1];
            });
        };

        for i in 0..ell {
            let step = 1 << (ell - i);
            let half_step = step / 2;

            let r_lower_product =
                (F::one() - r[i]) * r.iter().skip(i + 1).copied().product::<F<N>>();

            eq_plus_one_evals
                .iter_mut()
                .enumerate()
                .skip(half_step)
                .step_by(step)
                .for_each(|(index, v)| {
                    *v = eq_evals[index - half_step] * r_lower_product;
                });

            eq_evals_helper(&mut eq_evals, r, i + 1);
        }

        (eq_evals, eq_plus_one_evals)
    }
}
