#![allow(non_snake_case)]
use ark_ff::Zero;
use crypto_bigint::Int;

use crate::field::conversion::FieldMap;
use crate::field::RandomField as F;
use crate::field_config::FieldConfig;
use crate::zip::utils::expand;

use ark_std::fmt::Debug;

use itertools::Itertools;
use std::collections::BTreeSet;
use std::iter;

use super::pcs::structs::ZipTranscript;

const INVERSE_RATE: usize = 2;
pub trait LinearCodes<const N: usize, const M: usize>: Sync + Send {
    fn row_len(&self) -> usize;

    fn codeword_len(&self) -> usize;

    fn num_column_opening(&self) -> usize;

    fn num_proximity_testing(&self) -> usize;

    fn encode(&self, input: &[Int<N>]) -> Vec<Int<M>>;
}

#[derive(Clone, Debug)]
pub struct Zip<const N: usize, const L: usize> {
    row_len: usize,
    codeword_len: usize,
    num_column_opening: usize,
    num_proximity_testing: usize,
    a: SparseMatrixZ<L>,
    b: SparseMatrixZ<L>,
}

impl<const N: usize, const L: usize> Zip<N, L> {
    pub fn proof_size<S: ZipSpec>(n_0: usize, c: usize, r: usize) -> usize {
        let log2_q = N;
        let num_ldt = S::num_proximity_testing(log2_q, c, n_0);
        (1 + num_ldt) * c + S::num_column_opening() * r
    }

    pub fn new_multilinear<S: ZipSpec, T: ZipTranscript<L>>(
        num_vars: usize,
        n_0: usize,
        transcript: &mut T,
    ) -> Self {
        assert!(1 << num_vars > n_0);

        let log2_q = N;

        let row_len = ((1 << num_vars) as u64).isqrt().next_power_of_two() as usize;

        let codeword_len = S::codeword_len(row_len);

        let num_column_opening = S::num_column_opening();
        let num_proximity_testing = S::num_proximity_testing(log2_q, row_len, n_0);

        let (a, b) = S::matrices(codeword_len / 2, row_len, row_len / 2, transcript);
        Self {
            row_len,
            codeword_len,
            num_column_opening,
            num_proximity_testing,
            a,
            b,
        }
    }

    pub fn encode_f(&self, row: &[F<N>], field: *const FieldConfig<N>) -> Vec<F<N>> {
        let mut code = Vec::with_capacity(self.codeword_len);
        let a_f = SparseMatrixF::new(&self.a, field);
        let b_f = SparseMatrixF::new(&self.b, field);
        code.extend(SparseMatrixF::mat_vec_mul(&a_f, row));
        code.extend(SparseMatrixF::mat_vec_mul(&b_f, row));

        code
    }

    pub(crate) fn encode_wide<const M: usize>(&self, row: &[Int<M>]) -> Vec<Int<M>> {
        let mut code = Vec::with_capacity(self.codeword_len);
        code.extend(SparseMatrixZ::mat_vec_mul(&self.a, row));
        code.extend(SparseMatrixZ::mat_vec_mul(&self.b, row));

        code
    }
}

impl<const N: usize, const M: usize, const L: usize> LinearCodes<N, M> for Zip<N, L> {
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

    fn encode(&self, row: &[Int<N>]) -> Vec<Int<M>> {
        let mut code = Vec::with_capacity(self.codeword_len);
        code.extend(SparseMatrixZ::mat_vec_mul(&self.a, row));
        code.extend(SparseMatrixZ::mat_vec_mul(&self.b, row));

        code
    }
}

pub trait ZipSpec: Debug {
    fn num_column_opening() -> usize {
        1000
    }

    fn num_proximity_testing(_log2_q: usize, _n: usize, _n_0: usize) -> usize {
        1
    }

    fn codeword_len(n: usize) -> usize {
        n * INVERSE_RATE
    }

    fn matrices<const L: usize, T: ZipTranscript<L>>(
        rows: usize,
        cols: usize,
        density: usize,
        transcript: &mut T,
    ) -> (SparseMatrixZ<L>, SparseMatrixZ<L>) {
        let dim = SparseMatrixDimension::new(rows, cols, density);
        (
            SparseMatrixZ::new(dim, transcript),
            SparseMatrixZ::new(dim, transcript),
        )
    }
}

macro_rules! impl_spec_128 {
    ($(($name:ident,)),*) => {
        $(
            #[derive(Debug)]
            pub struct $name;
            impl ZipSpec for $name {

            }
        )*
    };
}

// Figure 2 in [GLSTW21](https://eprint.iacr.org/2021/1043.pdf).
impl_spec_128!((ZipSpec1,));

#[derive(Clone, Copy, Debug)]
pub struct SparseMatrixDimension {
    n: usize, // number of rows
    m: usize, // number of columns
    d: usize, // number of non-zero elements per row
}

impl std::fmt::Display for SparseMatrixDimension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}x{} matrix with {} non-zero elements per row",
            self.n, self.m, self.d
        )
    }
}

impl SparseMatrixDimension {
    fn new(n: usize, m: usize, d: usize) -> Self {
        Self { n, m, d }
    }
}

#[derive(Clone, Debug)]
pub struct SparseMatrixZ<const L: usize> {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, Int<L>)>,
}

impl<const L: usize> SparseMatrixZ<L> {
    fn new<T: ZipTranscript<L>>(dimension: SparseMatrixDimension, transcript: &mut T) -> Self {
        let cells = iter::repeat_with(|| {
            let mut columns = BTreeSet::<usize>::new();
            transcript.sample_unique_columns(0..dimension.m, &mut columns, dimension.d);
            columns
                .into_iter()
                .map(|column| (column, transcript.get_encoding_element()))
                .collect_vec()
        })
        .take(dimension.n)
        .flatten()
        .collect();
        Self { dimension, cells }
    }

    fn rows(&self) -> impl Iterator<Item = &[(usize, Int<L>)]> {
        self.cells.chunks(self.dimension.d)
    }

    fn mat_vec_mul<const N: usize, const M: usize>(&self, vector: &[Int<N>]) -> Vec<Int<M>> {
        assert_eq!(
            self.dimension.m,
            vector.len(),
            "Vector length must match matrix column dimension"
        );

        let mut result = vec![Int::from(0); self.dimension.n];

        self.rows().enumerate().for_each(|(row_idx, cells)| {
            let mut sum = Int::<M>::zero();
            for (column, coeff) in cells.iter() {
                sum += expand::<L, M>(coeff) * expand::<N, M>(&vector[*column]);
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
    fn new<const L: usize>(
        sparse_matrix: &SparseMatrixZ<L>,
        config: *const FieldConfig<N>,
    ) -> Self {
        let cells_f: Vec<(usize, F<N>)> = sparse_matrix
            .cells
            .iter()
            .map(|(col_index, val)| (*col_index, val.map_to_field(config)))
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
