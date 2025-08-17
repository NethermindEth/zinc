#![allow(non_snake_case)]

use std::{collections::BTreeSet, fmt::Debug, iter};

use crypto_bigint::Int;
use itertools::Itertools;

use crate::{
    traits::{Field, Transcript},
    zip::utils::expand,
};

pub trait LinearCode<const N: usize, const L: usize, const K: usize, const M: usize>:
    Sync + Send
{
    /// Length of each input row before encoding
    fn row_len(&self) -> usize;

    /// Length of each encoded codeword (output length after encoding)
    fn codeword_len(&self) -> usize;

    /// Number of columns to open during verification (security parameter)
    fn num_column_opening(&self) -> usize;

    /// Number of proximity tests to perform (security parameter)
    fn num_proximity_testing(&self) -> usize;

    /// Encodes a row of cryptographic integers using this linear encoding scheme.
    ///
    /// This function is optimized for the prover's context where we work with cryptographic integers.
    /// It's more efficient than `encode_f` as it avoids field conversions.
    ///
    /// # Parameters
    /// - `row`: Slice of cryptographic integers to encode
    ///
    /// # Returns
    /// A vector of cryptographic integers representing the encoded row
    fn encode(&self, row: &[Int<N>]) -> Vec<Int<M>> {
        self.encode_wide(row)
    }

    /// Encodes a row of cryptographic integers using this linear encoding scheme.
    ///
    /// This function is optimized for the prover's context where we work with cryptographic integers.
    /// It's more efficient than `encode_f` as it avoids field conversions.
    ///
    /// # Parameters
    /// - `row`: Slice of cryptographic integers to encode
    ///
    /// # Returns
    /// A vector of cryptographic integers representing the encoded row
    fn encode_wide<const IN: usize, const OUT: usize>(&self, row: &[Int<IN>]) -> Vec<Int<OUT>>;

    /// Encodes a row of field elements using this linear encoding scheme.
    ///
    /// This function is used when working with field elements directly and performs the encoding
    /// by first converting the sparse matrices to field elements.
    ///
    /// # Parameters
    /// - `row`: Slice of field elements to encode
    /// - `field`: Field configuration for the conversion
    ///
    /// # Returns
    /// A vector of field elements representing the encoded row
    fn encode_f<F: Field<LIMBS>, const LIMBS: usize>(&self, row: &[F]) -> Vec<F>;
}

/// A linear code implementation used for the Zip PCS.
///
/// # Type Parameters
/// - `I`: The input cryptographic integer type. Represents the field elements being encoded.
/// - `L`: The matrix element type. A larger cryptographic integer type used for sparse matrix
///   operations to prevent overflow during encoding. Must be at least as large as `I`.
#[derive(Clone, Debug)]
pub struct ZipLinearCode<const N: usize, const L: usize, const K: usize, const M: usize> {
    /// Length of each input row before encoding
    row_len: usize,

    /// Length of each encoded codeword (output length after encoding)
    codeword_len: usize,

    /// Number of columns to open during verification (security parameter)
    num_column_opening: usize,

    /// Number of proximity tests to perform (security parameter)
    num_proximity_testing: usize,

    /// First sparse matrix used in the encoding process
    a: SparseMatrixZ<L>,

    /// Second sparse matrix used in the encoding process
    b: SparseMatrixZ<L>,
}

impl<const N: usize, const L: usize, const K: usize, const M: usize> ZipLinearCode<N, L, K, M> {
    pub fn new<S: LinearCodeSpec, T: Transcript<L>>(
        spec: &S,
        poly_size: usize,
        transcript: &mut T,
    ) -> Self {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        Self::new_multilinear::<S, T>(spec, num_vars, 20.min((1 << num_vars) - 1), transcript)
    }

    /// Creates a new linear code instance for multilinear polynomials.
    ///
    /// # Parameters
    /// - `num_vars`: Number of variables in the multilinear polynomial
    /// - `n_0`: Number of rows in the matrix representation of the polynomial
    /// - `transcript`: Reference to a transcript for generating random challenges
    fn new_multilinear<S: LinearCodeSpec, T: Transcript<L>>(
        spec: &S,
        num_vars: usize,
        n_0: usize,
        transcript: &mut T,
    ) -> Self {
        assert!(1 << num_vars > n_0);

        let log2_q = N;

        let row_len = ((1 << num_vars) as u64).isqrt().next_power_of_two() as usize;

        let codeword_len = row_len * spec.repetition_factor();

        let num_column_opening = spec.num_column_opening();
        let num_proximity_testing = spec.num_proximity_testing(log2_q, row_len, n_0);

        let (a, b) = Self::matrices(codeword_len / 2, row_len, row_len / 2, transcript);
        Self {
            row_len,
            codeword_len,
            num_column_opening,
            num_proximity_testing,
            a,
            b,
        }
    }

    pub fn proof_size<S: LinearCodeSpec>(spec: S, n_0: usize, c: usize, r: usize) -> usize {
        let log2_q = N;
        // Number of low-degree tests
        let num_ldt = spec.num_proximity_testing(log2_q, c, n_0);
        (1 + num_ldt) * c + spec.num_column_opening() * r
    }

    fn matrices<const L1: usize, T: Transcript<L1>>(
        rows: usize,
        cols: usize,
        density: usize,
        transcript: &mut T,
    ) -> (SparseMatrixZ<L1>, SparseMatrixZ<L1>) {
        let dim = SparseMatrixDimension::new(rows, cols, density);
        (
            SparseMatrixZ::sample_new(dim, transcript),
            SparseMatrixZ::sample_new(dim, transcript),
        )
    }
}

impl<const N: usize, const L: usize, const K: usize, const M: usize> LinearCode<N, L, K, M>
    for ZipLinearCode<N, L, K, M>
{
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

    fn encode_wide<const IN: usize, const OUT: usize>(&self, row: &[Int<IN>]) -> Vec<Int<OUT>> {
        debug_assert_eq!(
            row.len(),
            self.row_len,
            "Row length must match the code's row length"
        );
        let mut code = Vec::with_capacity(self.codeword_len);
        code.extend(self.a.mat_vec_mul::<IN, OUT>(row));
        code.extend(self.b.mat_vec_mul::<IN, OUT>(row));
        code
    }

    fn encode_f<F: Field<LIMBS>, const LIMBS: usize>(&self, row: &[F]) -> Vec<F> {
        debug_assert_eq!(
            row.len(),
            self.row_len,
            "Row length must match the code's row length"
        );
        let mut code = Vec::with_capacity(self.codeword_len);
        let a_f = SparseMatrixF::new(&self.a);
        let b_f = SparseMatrixF::new(&self.b);
        code.extend(a_f.mat_vec_mul(row));
        code.extend(b_f.mat_vec_mul(row));
        code
    }
}

pub trait LinearCodeSpec: Debug {
    fn num_column_opening(&self) -> usize;

    /// A.k.a. inverse rate, the ratio of codeword length to input row length.
    /// Has to be at a power of 2.
    fn repetition_factor(&self) -> usize;

    fn num_proximity_testing(&self, _log2_q: usize, _n: usize, _n_0: usize) -> usize;
}

// Figure 2 in [GLSTW21](https://eprint.iacr.org/2021/1043.pdf).
#[derive(Debug)]
pub struct DefaultLinearCodeSpec;
impl LinearCodeSpec for DefaultLinearCodeSpec {
    fn num_column_opening(&self) -> usize {
        1000
    }

    fn repetition_factor(&self) -> usize {
        2
    }

    fn num_proximity_testing(&self, _log2_q: usize, _n: usize, _n_0: usize) -> usize {
        1
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SparseMatrixDimension {
    /// Number of rows
    n: usize,
    /// Number of columns
    m: usize,
    /// Number of non-zero elements per row
    d: usize,
}

impl ark_std::fmt::Display for SparseMatrixDimension {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
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

/// Sparse matrix over a ring of integers.
#[derive(Clone, Debug)]
pub struct SparseMatrixZ<const I: usize> {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, Int<I>)>,
}

impl<const L: usize> SparseMatrixZ<L> {
    /// Creates a new sparse matrix with the given dimension and samples its cells using the
    /// provided transcript.
    fn sample_new<T: Transcript<L>>(dimension: SparseMatrixDimension, transcript: &mut T) -> Self {
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

    pub fn rows(&self) -> impl Iterator<Item = &[(usize, Int<L>)]> {
        self.cells.chunks(self.dimension.d)
    }

    /// Multiplies the sparse matrix by a vector of cryptographic integers.
    pub fn mat_vec_mul<const N: usize, const M: usize>(&self, vector: &[Int<N>]) -> Vec<Int<M>> {
        assert_eq!(
            self.dimension.m,
            vector.len(),
            "Vector length must match matrix column dimension"
        );

        let mut result = vec![Int::<M>::from_i64(0i64); self.dimension.n];

        self.rows().enumerate().for_each(|(row_idx, cells)| {
            let mut sum = Int::<M>::ZERO;
            for (column, coeff) in cells.iter() {
                sum += &(expand::<L, M>(coeff) * expand::<N, M>(&vector[*column]));
            }
            result[row_idx] = sum;
        });

        result
    }

    pub fn to_dense(&self) -> Vec<Vec<Int<L>>> {
        let mut r: Vec<Vec<Int<L>>> =
            vec![vec![Int::<L>::ZERO; self.dimension.m]; self.dimension.n];
        for (row_i, (col_i, value)) in self.cells.iter().enumerate() {
            r[row_i][*col_i] = *value;
        }
        r
    }
}

/// Sparse matrix over a field.
#[derive(Clone, Debug)]
pub struct SparseMatrixF<F: Field<LIMBS>, const LIMBS: usize> {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, F)>,
}

impl<F: Field<LIMBS>, const LIMBS: usize> SparseMatrixF<F, LIMBS> {
    pub fn new<const L: usize>(sparse_matrix: &SparseMatrixZ<L>) -> Self {
        let cells_f: Vec<(usize, F)> = sparse_matrix
            .cells
            .iter()
            .map(|(col_index, val)| (*col_index, val.resize().into()))
            .collect();
        Self {
            dimension: sparse_matrix.dimension,
            cells: cells_f,
        }
    }

    fn rows(&self) -> impl Iterator<Item = &[(usize, F)]> {
        self.cells.chunks(self.dimension.d)
    }

    /// Multiplies the sparse matrix by a vector of cryptographic integers.
    pub fn mat_vec_mul(&self, vector: &[F]) -> Vec<F> {
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
