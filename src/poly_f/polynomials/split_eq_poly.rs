//! Implements the Dao-Thaler optimization for EQ polynomial evaluations
//! https://eprint.iacr.org/2024/1210.pdf

use crate::field::RandomField as F;
use crate::field_config::FieldConfig;
use crate::poly_f::polynomials::DenseMultilinearExtension;

use super::eq_poly::EqPolynomial;
#[derive(Debug, Clone, PartialEq)]
pub struct SplitEqPolynomial<const N: usize> {
    num_vars: usize,
    config: *const FieldConfig<N>,
    pub(crate) E1: Vec<F<N>>,
    pub(crate) E1_len: usize,
    pub(crate) E2: Vec<F<N>>,
    pub(crate) E2_len: usize,
}

impl<const N: usize> SplitEqPolynomial<N> {
    pub fn new(w: &[F<N>], config: *const FieldConfig<N>) -> Self {
        let m = w.len() / 2;
        let (w2, w1) = w.split_at(m);
        let (E2, E1) = (EqPolynomial::evals(w2), EqPolynomial::evals(w1));
        let E1_len = E1.len();
        let E2_len = E2.len();
        Self {
            num_vars: w.len(),
            config,
            E1,
            E1_len,
            E2,
            E2_len,
        }
    }

    pub fn get_num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn len(&self) -> usize {
        if self.E1_len == 1 {
            self.E2_len
        } else {
            self.E1_len * self.E2_len
        }
    }

    pub fn bind(&mut self, r: F<N>) {
        if self.E1_len == 1 {
            // E_1 is already completely bound, so we bind E_2
            let n = self.E2_len / 2;
            for i in 0..n {
                self.E2[i] = self.E2[2 * i] + r * (self.E2[2 * i + 1] - self.E2[2 * i]);
            }
            self.E2_len = n;
        } else {
            // Bind E_1
            let n = self.E1_len / 2;
            for i in 0..n {
                self.E1[i] = self.E1[2 * i] + r * (self.E1[2 * i + 1] - self.E1[2 * i]);
            }
            self.E1_len = n;

            // If E_1 is now completely bound, we will be switching over to the
            // linear-time sumcheck prover, using E_1 * E_2:
            if self.E1_len == 1 {
                self.E2[..self.E2_len]
                    .iter_mut()
                    .for_each(|eval| *eval *= self.E1[0]);
            }
        }
    }

    #[cfg(test)]
    pub fn merge(&self) -> DenseMultilinearExtension<N> {
        if self.E1_len == 1 {
            DenseMultilinearExtension::from_evaluations_vec(
                self.num_vars,
                self.E2[..self.E2_len].to_vec(),
                self.config,
            )
        } else {
            let mut merged = vec![];
            for i in 0..self.E2_len {
                for j in 0..self.E1_len {
                    merged.push(self.E2[i] * self.E1[j])
                }
            }
            DenseMultilinearExtension::from_evaluations_vec(self.num_vars, merged, self.config)
        }
    }
}
