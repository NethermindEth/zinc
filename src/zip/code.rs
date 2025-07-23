#![allow(non_snake_case)]

use ark_std::{collections::BTreeSet, fmt::Debug, iter, marker::PhantomData, vec, vec::Vec};
use itertools::Itertools;

use super::pcs::structs::ZipTranscript;
use crate::{
    traits::{Field, FieldMap, Integer, Words, ZipTypes},
    zip::utils::expand,
};

const INVERSE_RATE: usize = 2;
pub trait LinearCode<ZT: ZipTypes>: Sync + Send {
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
    fn encode(&self, row: &[ZT::N]) -> Vec<ZT::M> {
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
    fn encode_wide<In, Out>(&self, row: &[In]) -> Vec<Out>
    where
        In: Integer,
        Out: Integer + for<'a> From<&'a In> + for<'a> From<&'a ZT::L>;

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
    fn encode_f<F: Field>(&self, row: &[F], field: F::R) -> Vec<F>
    where
        ZT::L: FieldMap<F, Output = F>;
}

/// A linear code implementation used for the Zip PCS.
///
/// # Type Parameters
/// - `I`: The input cryptographic integer type. Represents the field elements being encoded.
/// - `L`: The matrix element type. A larger cryptographic integer type used for sparse matrix
///   operations to prevent overflow during encoding. Must be at least as large as `I`.
#[derive(Clone, Debug)]
pub struct ZipLinearCode<ZT: ZipTypes> {
    /// Length of each input row before encoding
    row_len: usize,

    /// Length of each encoded codeword (output length after encoding)
    codeword_len: usize,

    /// Number of columns to open during verification (security parameter)
    num_column_opening: usize,

    /// Number of proximity tests to perform (security parameter)
    num_proximity_testing: usize,

    /// First sparse matrix used in the encoding process
    a: SparseMatrixZ<ZT::L>,

    /// Second sparse matrix used in the encoding process
    b: SparseMatrixZ<ZT::L>,

    phantom: PhantomData<ZT>,
}

impl<ZT: ZipTypes> ZipLinearCode<ZT> {
    pub fn new<S: ZipLinearCodeSpec, T: ZipTranscript<ZT::L>>(
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
    fn new_multilinear<S: ZipLinearCodeSpec, T: ZipTranscript<ZT::L>>(
        spec: &S,
        num_vars: usize,
        n_0: usize,
        transcript: &mut T,
    ) -> Self {
        assert!(1 << num_vars > n_0);

        let log2_q = <ZT::N as Integer>::W::num_words();

        let row_len = ((1 << num_vars) as u64).isqrt().next_power_of_two() as usize;

        let codeword_len = spec.codeword_len(row_len);

        let num_column_opening = spec.num_column_opening();
        let num_proximity_testing = spec.num_proximity_testing(log2_q, row_len, n_0);

        let (a, b) = spec.matrices(codeword_len / 2, row_len, row_len / 2, transcript);
        Self {
            row_len,
            codeword_len,
            num_column_opening,
            num_proximity_testing,
            a,
            b,
            phantom: PhantomData,
        }
    }

    pub fn proof_size<S: ZipLinearCodeSpec>(spec: S, n_0: usize, c: usize, r: usize) -> usize {
        let log2_q = <ZT::N as Integer>::W::num_words();
        // Number of low-degree tests
        let num_ldt = spec.num_proximity_testing(log2_q, c, n_0);
        (1 + num_ldt) * c + spec.num_column_opening() * r
    }
}

impl<ZT: ZipTypes> LinearCode<ZT> for ZipLinearCode<ZT> {
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

    fn encode_wide<In, Out>(&self, row: &[In]) -> Vec<Out>
    where
        In: Integer,
        Out: Integer + for<'a> From<&'a In> + for<'a> From<&'a ZT::L>,
    {
        let mut code = Vec::with_capacity(self.codeword_len);
        code.extend(self.a.mat_vec_mul(row));
        code.extend(self.b.mat_vec_mul(row));
        code
    }

    fn encode_f<F: Field>(&self, row: &[F], field: F::R) -> Vec<F>
    where
        ZT::L: FieldMap<F, Output = F>,
    {
        let mut code = Vec::with_capacity(self.codeword_len);
        let a_f = SparseMatrixF::new(&self.a, field);
        let b_f = SparseMatrixF::new(&self.b, field);
        code.extend(a_f.mat_vec_mul(row));
        code.extend(b_f.mat_vec_mul(row));
        code
    }
}

pub trait ZipLinearCodeSpec: Debug {
    fn num_column_opening(&self) -> usize {
        1000
    }

    fn num_proximity_testing(&self, _log2_q: usize, _n: usize, _n_0: usize) -> usize {
        1
    }

    fn codeword_len(&self, n: usize) -> usize {
        n * INVERSE_RATE
    }

    fn matrices<L: Integer, T: ZipTranscript<L>>(
        &self,
        rows: usize,
        cols: usize,
        density: usize,
        transcript: &mut T,
    ) -> (SparseMatrixZ<L>, SparseMatrixZ<L>) {
        let dim = SparseMatrixDimension::new(rows, cols, density);
        (
            SparseMatrixZ::sample_new(dim, transcript),
            SparseMatrixZ::sample_new(dim, transcript),
        )
    }
}

macro_rules! impl_spec_128 {
    ($(($name:ident,)),*) => {
        $(
            #[derive(Debug)]
            pub struct $name;
            impl ZipLinearCodeSpec for $name {

            }
        )*
    };
}

// Figure 2 in [GLSTW21](https://eprint.iacr.org/2021/1043.pdf).
impl_spec_128!((ZipLinearCodeSpec1,));

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
pub struct SparseMatrixZ<I: Integer> {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, I)>,
}

impl<L: Integer> SparseMatrixZ<L> {
    /// Creates a new sparse matrix with the given dimension and samples its cells using the
    /// provided transcript.
    fn sample_new<T: ZipTranscript<L>>(
        dimension: SparseMatrixDimension,
        transcript: &mut T,
    ) -> Self {
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

    fn rows(&self) -> impl Iterator<Item = &[(usize, L)]> {
        self.cells.chunks(self.dimension.d)
    }

    /// Multiplies the sparse matrix by a vector of cryptographic integers.
    fn mat_vec_mul<N: Integer, M: Integer + for<'a> From<&'a N> + for<'a> From<&'a L>>(
        &self,
        vector: &[N],
    ) -> Vec<M> {
        assert_eq!(
            self.dimension.m,
            vector.len(),
            "Vector length must match matrix column dimension"
        );

        let mut result = vec![M::from_i64(0i64); self.dimension.n];

        self.rows().enumerate().for_each(|(row_idx, cells)| {
            let mut sum = M::ZERO;
            for (column, coeff) in cells.iter() {
                sum += &(expand::<L, M>(coeff) * expand::<N, M>(&vector[*column]));
            }
            result[row_idx] = sum;
        });

        result
    }
}

/// Sparse matrix over a field.
#[derive(Clone, Debug)]
pub struct SparseMatrixF<F: Field> {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, F)>,
}

impl<F: Field> SparseMatrixF<F> {
    fn new<L: Integer + FieldMap<F, Output = F>>(
        sparse_matrix: &SparseMatrixZ<L>,
        config: F::R,
    ) -> Self {
        let cells_f: Vec<(usize, F)> = sparse_matrix
            .cells
            .iter()
            .map(|(col_index, val)| (*col_index, val.map_to_field(config)))
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
    fn mat_vec_mul(&self, vector: &[F]) -> Vec<F> {
        assert_eq!(
            self.dimension.m,
            vector.len(),
            "Vector length must match matrix column dimension"
        );

        let mut result = vec![F::zero(); self.dimension.n];

        self.rows().enumerate().for_each(|(row_idx, cells)| {
            let mut sum = F::zero();
            for (column, coeff) in cells.iter() {
                sum += &(coeff.clone() * &vector[*column]);
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
