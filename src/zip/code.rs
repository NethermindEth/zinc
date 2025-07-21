#![allow(non_snake_case)]

use ark_std::{collections::BTreeSet, fmt::Debug, iter, marker::PhantomData, vec, vec::Vec};
use itertools::Itertools;

use super::pcs::structs::ZipTranscript;
use crate::{
    traits::{Field, FieldMap, Integer, Words},
    zip::utils::expand,
};

const INVERSE_RATE: usize = 2;
pub trait LinearCodes<In: Integer, Im: Integer>: Sync + Send {
    fn row_len(&self) -> usize;

    fn codeword_len(&self) -> usize;

    fn num_column_opening(&self) -> usize;

    fn num_proximity_testing(&self) -> usize;

    fn encode(&self, input: &[In]) -> Vec<Im>;
}

#[derive(Clone, Debug)]
pub struct Zip<I: Integer, L: Integer> {
    row_len: usize,
    codeword_len: usize,
    num_column_opening: usize,
    num_proximity_testing: usize,
    a: SparseMatrixZ<L>,
    b: SparseMatrixZ<L>,
    phantom: PhantomData<I>,
}

impl<I: Integer, L: Integer> Zip<I, L> {
    pub fn proof_size<S: ZipSpec>(n_0: usize, c: usize, r: usize) -> usize {
        let log2_q = I::W::num_words();
        let num_ldt = S::num_proximity_testing(log2_q, c, n_0);
        (1 + num_ldt) * c + S::num_column_opening() * r
    }

    pub fn new_multilinear<S: ZipSpec, T: ZipTranscript<L>>(
        num_vars: usize,
        n_0: usize,
        transcript: &mut T,
    ) -> Self {
        assert!(1 << num_vars > n_0);

        let log2_q = I::W::num_words();

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
            phantom: PhantomData,
        }
    }

    pub fn encode_f<F: Field>(&self, row: &[F], field: F::R) -> Vec<F>
    where
        L: FieldMap<F, Output = F>,
    {
        let mut code = Vec::with_capacity(self.codeword_len);
        let a_f = SparseMatrixF::new(&self.a, field);
        let b_f = SparseMatrixF::new(&self.b, field);
        code.extend(SparseMatrixF::mat_vec_mul(&a_f, row));
        code.extend(SparseMatrixF::mat_vec_mul(&b_f, row));

        code
    }

    pub(crate) fn encode_wide<M: Integer + for<'a> From<&'a L> + for<'a> From<&'a M>>(
        &self,
        row: &[M],
    ) -> Vec<M> {
        let mut code = Vec::with_capacity(self.codeword_len);
        code.extend(SparseMatrixZ::mat_vec_mul(&self.a, row));
        code.extend(SparseMatrixZ::mat_vec_mul(&self.b, row));

        code
    }
}

impl<N: Integer, M: Integer + for<'a> From<&'a N> + for<'a> From<&'a L>, L: Integer>
    LinearCodes<N, M> for Zip<N, L>
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

    fn encode(&self, row: &[N]) -> Vec<M> {
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

    fn matrices<L: Integer, T: ZipTranscript<L>>(
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

#[derive(Clone, Debug)]
pub struct SparseMatrixZ<I: Integer> {
    dimension: SparseMatrixDimension,
    cells: Vec<(usize, I)>,
}

impl<L: Integer> SparseMatrixZ<L> {
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

    fn rows(&self) -> impl Iterator<Item = &[(usize, L)]> {
        self.cells.chunks(self.dimension.d)
    }

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
