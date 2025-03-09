use i256::I256;
use itertools::Itertools;
use std::collections::BTreeSet;
use std::iter;

use ark_ff::UniformRand;
use ark_std::fmt::Debug;
use ark_std::rand::distributions::Uniform;
use ark_std::rand::Rng;
use ark_std::rand::RngCore;
#[allow(dead_code)]
const PROB_MULTIPLIER: usize = 18;
#[allow(dead_code)]
const INVERSE_RATE: usize = 2;
pub trait LinearCodes<const N: usize>: Sync + Send {
    fn row_len(&self) -> usize;

    fn codeword_len(&self) -> usize;

    fn num_column_opening(&self) -> usize;

    fn num_proximity_testing(&self) -> usize;

    fn encode(&self, input: &[i64]) -> Vec<I256>;
}

#[derive(Clone, Debug)]
pub struct Zip<const N: usize> {
    row_len: usize,
    codeword_len: usize,
    num_column_opening: usize,
    num_proximity_testing: usize,
    a: SparseMatrix,
    b: SparseMatrix,
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
        let min_log2_n = (n_0 + 1).next_power_of_two().ilog2() as usize;

        let (_, row_len) =
            (min_log2_n..=num_vars).fold((usize::MAX, 0), |(min_proof_size, row_len), log2_n| {
                let proof_size = Self::proof_size::<S>(n_0, 1 << log2_n, 1 << (num_vars - log2_n));
                if proof_size < min_proof_size {
                    (proof_size, 1 << log2_n)
                } else {
                    (min_proof_size, row_len)
                }
            });

        let codeword_len = S::codeword_len(row_len);
        let num_column_opening = S::num_column_opening();
        let num_proximity_testing = S::num_proximity_testing(log2_q, row_len, n_0);
        let (a, b) = S::matrices(log2_q, row_len, n_0, rng);

        Self {
            row_len,
            codeword_len,
            num_column_opening,
            num_proximity_testing,
            a,
            b,
        }
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

    fn encode(&self, row: &[i64]) -> Vec<I256> {
        let mut code = Vec::with_capacity(self.codeword_len);
        code.extend(SparseMatrix::mat_vec_mul(&self.a, row));
        code.extend(SparseMatrix::mat_vec_mul(&self.b, row));
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
        2617
    }

    fn num_proximity_testing(log2_q: usize, n: usize, _n_0: usize) -> usize {
        ceil(Self::LAMBDA / (log2_q as f64 - (Self::codeword_len(n) as f64).log2()))
    }

    fn codeword_len(n: usize) -> usize {
        n * INVERSE_RATE
    }

    fn matrices(
        rows: usize,
        cols: usize,
        density: usize,
        mut rng: impl RngCore,
    ) -> (SparseMatrix, SparseMatrix) {
        let dim = SparseMatrixDimension::new(rows, cols * INVERSE_RATE, density);
        (
            SparseMatrix::new(dim, &mut rng),
            SparseMatrix::new(dim, &mut rng),
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
pub struct SparseMatrix {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, i128)>,
}

impl SparseMatrix {
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

    // fn dot(&self, array: &[i64]) -> Vec<i64> {
    //     let mut target = vec![0i64; self.dimension.m];
    //     self.dot_into(array, &mut target);
    //     target
    // }

    fn mat_vec_mul(&self, vector: &[i64]) -> Vec<I256> {
        assert_eq!(
            self.dimension.m,
            vector.len(),
            "Vector length must match matrix column dimension"
        );

        let mut result = vec![I256::from(0); self.dimension.n];

        self.rows().enumerate().for_each(|(row_idx, cells)| {
            let mut sum = I256::from(0);
            for (column, coeff) in cells.iter() {
                sum += I256::from(*coeff) * (I256::from(vector[*column]));
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
