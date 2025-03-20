#![allow(non_snake_case)]
use ark_ff::Zero;
use i256::{I256, I512};

use crate::field::RandomField as F;
use crate::field_config::FieldConfig;

use ark_ff::UniformRand;
use ark_std::fmt::Debug;
use ark_std::rand::distributions::Uniform;
use ark_std::rand::Rng;
use ark_std::rand::RngCore;
use itertools::Itertools;
use std::collections::BTreeSet;
use std::iter;
#[allow(dead_code)]
const PROB_MULTIPLIER: usize = 18;
#[allow(dead_code)]
const INVERSE_RATE: usize = 2;
pub trait LinearCodes<const N: usize>: Sync + Send {
    fn row_len(&self) -> usize;

    fn codeword_len(&self) -> usize;

    fn num_column_opening(&self) -> usize;

    fn num_proximity_testing(&self) -> usize;

    fn encode(&self, input: &[I256]) -> Vec<I512>;
}

#[derive(Clone, Debug)]
pub struct Zip<const N: usize> {
    row_len: usize,
    codeword_len: usize,
    num_column_opening: usize,
    num_proximity_testing: usize,
    a: SparseMatrixZ,
    b: SparseMatrixZ,
}

impl<const N: usize> Zip<N> {
    pub fn proof_size<S: ZipSpec>(n_0: usize, c: usize, r: usize) -> usize {
        let log2_q = N;
        let num_ldt = S::num_proximity_testing(log2_q, c, n_0);
        (1 + num_ldt) * c + S::num_column_opening() * r
    }

    pub fn new_multilinear<S: ZipSpec>(num_vars: usize, n_0: usize, rng: impl RngCore) -> Self {
        assert!(1 << num_vars > n_0);

        let log2_q = N;

        let row_len = num_vars.pow(2).isqrt().next_power_of_two();

        let codeword_len = S::codeword_len(row_len);

        let num_column_opening = S::num_column_opening();
        let num_proximity_testing = S::num_proximity_testing(log2_q, row_len, n_0);

        let (a, b) = S::matrices(codeword_len / 2, row_len, row_len / 2, rng);
        Self {
            row_len,
            codeword_len,
            num_column_opening,
            num_proximity_testing,
            a,
            b,
        }
    }
    pub fn encode_i64(&self, row: &[i64]) -> Vec<I512> {
        let wider_row: Vec<_> = row.iter().map(|i| I256::from(*i)).collect();
        Self::encode(self, &wider_row)
    }
    pub fn encode_f(&self, row: &[F<N>], field: *const FieldConfig<N>) -> Vec<F<N>> {
        let mut code = Vec::with_capacity(self.codeword_len);
        let a_f = SparseMatrixF::new(&self.a, field);
        let b_f = SparseMatrixF::new(&self.b, field);
        code.extend(SparseMatrixF::mat_vec_mul(&a_f, row));
        code.extend(SparseMatrixF::mat_vec_mul(&b_f, row));

        code
    }
}

impl<const N: usize> LinearCodes<N> for Zip<N> {
    fn row_len(&self) -> usize {
        self.row_len
    }

    fn codeword_len(&self) -> usize {
        self.codeword_len
    }

    fn num_column_opening(&self) -> usize {
        self.num_column_opening
    }

    fn num_proximity_testing(&self) -> usize {
        self.num_proximity_testing
    }

    fn encode(&self, row: &[I256]) -> Vec<I512> {
        let mut code = Vec::with_capacity(self.codeword_len);
        code.extend(SparseMatrixZ::mat_vec_mul(&self.a, row));
        code.extend(SparseMatrixZ::mat_vec_mul(&self.b, row));

        code
    }
}

pub trait ZipSpec: Debug {
    const LAMBDA: f64;
    const ALPHA: f64;
    const BETA: f64;
    const R: f64;

    fn delta() -> f64 {
        Self::BETA / Self::R
    }

    fn mu() -> f64 {
        Self::R - 1f64 - Self::R * Self::ALPHA
    }

    fn nu() -> f64 {
        Self::BETA + Self::ALPHA * Self::BETA + 0.03
    }

    fn c_n(n: usize) -> usize {
        let alpha = Self::ALPHA;
        let beta = Self::BETA;
        let n = n as f64;
        std::cmp::min(
            std::cmp::max(ceil(1.28 * beta * n), ceil(beta * n) + 4),
            ceil(
                ((110.0 / n) + h(beta) + alpha * h(1.28 * beta / alpha))
                    / (beta * (alpha / (1.28 * beta)).log2()),
            ),
        )
    }

    fn d_n(log2_q: usize, n: usize) -> usize {
        let alpha = Self::ALPHA;
        let beta = Self::BETA;
        let r = Self::R;
        let mu = Self::mu();
        let nu = Self::nu();
        let log2_q = log2_q as f64;
        let n = n as f64;
        std::cmp::min(
            ceil((2.0 * beta + ((r - 1.0) + 110.0 / n) / log2_q) * n),
            ceil(
                (r * alpha * h(beta / r) + mu * h(nu / mu) + 110.0 / n)
                    / (alpha * beta * (mu / nu).log2()),
            ),
        )
    }

    fn num_column_opening() -> usize {
        1000
    }

    fn num_proximity_testing(_log2_q: usize, _n: usize, _n_0: usize) -> usize {
        1
    }

    fn codeword_len(n: usize) -> usize {
        n * INVERSE_RATE
    }

    fn matrices(
        rows: usize,
        cols: usize,
        density: usize,
        mut rng: impl RngCore,
    ) -> (SparseMatrixZ, SparseMatrixZ) {
        let dim = SparseMatrixDimension::new(rows, cols, density);
        (
            SparseMatrixZ::new(dim, &mut rng),
            SparseMatrixZ::new(dim, &mut rng),
        )
    }
}

macro_rules! impl_spec_128 {
    ($(($name:ident, $alpha:literal, $beta:literal, $r:literal)),*) => {
        $(
            #[derive(Debug)]
            pub struct $name;
            impl ZipSpec for $name {
                const LAMBDA: f64 = 128.0;
                const ALPHA: f64 = $alpha;
                const BETA: f64 = $beta;
                const R: f64 = $r;
            }
        )*
    };
}

// Figure 2 in [GLSTW21](https://eprint.iacr.org/2021/1043.pdf).
impl_spec_128!(
    (ZipSpec1, 0.1195, 0.0284, 1.420),
    (ZipSpec2, 0.1380, 0.0444, 1.470),
    (ZipSpec3, 0.1780, 0.0610, 1.521),
    (ZipSpec4, 0.2000, 0.0820, 1.640),
    (ZipSpec5, 0.2110, 0.0970, 1.616),
    (ZipSpec6, 0.2380, 0.1205, 1.720)
);

#[derive(Clone, Copy, Debug)]
pub struct SparseMatrixDimension {
    n: usize,
    m: usize,
    d: usize,
}

impl SparseMatrixDimension {
    fn new(n: usize, m: usize, d: usize) -> Self {
        Self { n, m, d }
    }
}

#[derive(Clone, Debug)]
pub struct SparseMatrixZ {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, i128)>,
}

impl SparseMatrixZ {
    fn new(dimension: SparseMatrixDimension, mut rng: impl RngCore) -> Self {
        let cells = iter::repeat_with(|| {
            let mut columns = BTreeSet::<usize>::new();
            (&mut rng)
                .sample_iter(&Uniform::new(0, dimension.m))
                .filter(|column| columns.insert(*column))
                .take(dimension.d)
                .count();
            columns
                .into_iter()
                .map(|column| (column, i128::rand(&mut rng)))
                .collect_vec()
        })
        .take(dimension.n)
        .flatten()
        .collect();
        Self { dimension, cells }
    }

    fn rows(&self) -> impl Iterator<Item = &[(usize, i128)]> {
        self.cells.chunks(self.dimension.d)
    }

    fn mat_vec_mul(&self, vector: &[I256]) -> Vec<I512> {
        assert_eq!(
            self.dimension.m,
            vector.len(),
            "Vector length must match matrix column dimension"
        );

        let mut result = vec![I512::from(0); self.dimension.n];

        self.rows().enumerate().for_each(|(row_idx, cells)| {
            let mut sum = I512::from(0);
            for (column, coeff) in cells.iter() {
                sum += I512::from(*coeff) * I256_to_I512(vector[*column]);
            }
            result[row_idx] = sum;
        });

        result
    }
}

#[derive(Clone, Debug)]
pub struct SparseMatrixF<const N: usize> {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, F<N>)>,
}

impl<const N: usize> SparseMatrixF<N> {
    fn new(sparse_matrix: &SparseMatrixZ, config: *const FieldConfig<N>) -> Self {
        let cells_f: Vec<(usize, F<N>)> = sparse_matrix
            .cells
            .iter()
            .map(|(col_index, val)| (*col_index, F::from_i128(*val, config)))
            .collect();
        Self {
            dimension: sparse_matrix.dimension,
            cells: cells_f,
        }
    }

    fn rows(&self) -> impl Iterator<Item = &[(usize, F<N>)]> {
        self.cells.chunks(self.dimension.d)
    }

    fn mat_vec_mul(&self, vector: &[F<N>]) -> Vec<F<N>> {
        assert_eq!(
            self.dimension.m,
            vector.len(),
            "Vector length must match matrix column dimension"
        );

        let mut result = vec![F::zero(); self.dimension.n];

        self.rows().enumerate().for_each(|(row_idx, cells)| {
            let mut sum = F::zero();
            for (column, coeff) in cells.iter() {
                sum += &(*coeff * vector[*column]);
            }
            result[row_idx] = sum;
        });

        result
    }
}

pub fn steps(start: i64) -> impl Iterator<Item = i64> {
    steps_by(start, 1i64)
}

pub fn steps_by(start: i64, step: i64) -> impl Iterator<Item = i64> {
    iter::successors(Some(start), move |state| Some(step + *state))
}

// H(p) = -p \log_2(p) - (1 - p) \log_2(1 - p)
fn h(p: f64) -> f64 {
    assert!(0f64 < p && p < 1f64);
    let one_minus_p = 1f64 - p;
    -p * p.log2() - one_minus_p * one_minus_p.log2()
}

fn ceil(v: f64) -> usize {
    v.ceil() as usize
}

pub(super) fn I256_to_I512(i: I256) -> I512 {
    let mut bytes = [0u8; 64];

    // Preserve sign extension for negative values
    let sign_byte = if i.is_negative() { 0xFF } else { 0x00 };
    bytes[..32].fill(sign_byte); // Fill the upper bytes with the sign
    bytes[32..].copy_from_slice(&i.to_be_bytes());

    I512::from_be_bytes(bytes)
}

#[cfg(test)]
mod tests {
    use i256::{I256, I512};

    use crate::zip::code::I256_to_I512;

    #[test]
    fn test_I256_to_I512_positive() {
        let i = I256::from(123456);
        let result = I256_to_I512(i);
        assert_eq!(result, I512::from(123456));
    }

    #[test]
    fn test_I256_to_I512_negative() {
        let i = I256::from(-987654);
        let result = I256_to_I512(i);
        assert_eq!(result, I512::from(-987654));
    }
}
