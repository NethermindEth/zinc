use crate::traits::{Integer, Field, FieldMap, Words, ZipTypes};
use crate::zip::code::{LinearCode, LinearCodeSpec, SparseMatrixF, SparseMatrixZ};
use crate::zip::pcs::structs::ZipTranscript;
use num_traits::Zero;
use rand::{Rng, SeedableRng};
use ark_std::fmt::Debug;
use ark_std::ops::{Add, AddAssign, Mul};
use ark_std::vec;
use ark_std::vec::Vec;

/// Implementation of a repeat-accumulate-accumulate (RAA) codes over the binary field,
/// as defined by the Blaze paper (https://eprint.iacr.org/2024/1609)
#[derive(Debug, Clone)]
pub struct RaaCode<ZT: ZipTypes> {
    row_len: usize,

    repetition_factor: usize,

    num_column_opening: usize,

    num_proximity_testing: usize,

    /// Permutation matrix (M_{\pi_1})
    perm_matrix_1: SparseMatrixZ<ZT::L>,

    /// Permutation matrix (M_{\pi_2})
    perm_matrix_2: SparseMatrixZ<ZT::L>,

    /// Matrix that repeats input vector N times, denoted as `F_{r}` in Blaze paper
    repetition_matrix: SparseMatrixZ<ZT::L>,

    encoding_matrix: DenseMatrixZ<ZT::L>,
}

impl<ZT: ZipTypes> RaaCode<ZT> {
    pub fn new<S: LinearCodeSpec, T: ZipTranscript<ZT::L>>(
        spec: &S,
        poly_size: usize,
        transcript: &mut T,
    ) -> Self {
        // Taken from original Zip codes
        let num_vars = poly_size.ilog2() as usize;
        let row_len = ((1 << num_vars) as u64).isqrt().next_power_of_two() as usize;
        let repetition_factor = spec.repetition_factor();
        let codeword_len = row_len * repetition_factor;

        // Note: Could be made more efficient by using u8 instead of ZT::L (all matrices are binary),
        // but since this is a setup phase, we don't really care.

        let acc_matrix: DenseMatrixZ<ZT::L> = gen_accumulation_matrix(codeword_len);

        let perm_matrix_1: SparseMatrixZ<ZT::L> =
            SparseMatrixZ::sample_permutation(codeword_len, transcript);

        let perm_matrix_2: SparseMatrixZ<ZT::L> =
            SparseMatrixZ::sample_permutation(codeword_len, transcript);

        let num_column_opening = spec.num_column_opening();
        let log2_q = <ZT::N as Integer>::W::num_words();
        let n_0 = 20.min((1 << num_vars) - 1);
        let num_proximity_testing = spec.num_proximity_testing(log2_q, row_len, n_0);

        let repetition_matrix: SparseMatrixZ<ZT::L> =
            SparseMatrixZ::repetition(row_len, repetition_factor);

        // Pre-multiply the matrices to create the encoding matrix
        let encoding_matrix = &acc_matrix;
        let encoding_matrix = encoding_matrix.mul_mod_2(&perm_matrix_2.to_dense());
        let encoding_matrix = encoding_matrix.mul_mod_2(&acc_matrix.rows);
        let encoding_matrix = encoding_matrix.mul_mod_2(&perm_matrix_1.to_dense());
        let mut encoding_matrix = encoding_matrix.mul_mod_2(&repetition_matrix.to_dense());

        // Apply a trick to randomize the signs of the encoding matrix.
        // This does not change the code's properties, but produces codewords with
        // quadratically smaller entries on average.
        encoding_matrix.randomize_sign(transcript);

        Self {
            row_len,
            repetition_factor,
            num_column_opening,
            num_proximity_testing,
            perm_matrix_1,
            perm_matrix_2,
            repetition_matrix,
            encoding_matrix,
        }
    }
}

/// Generate an accumulation matrix for the RAA code. It is a lower triangular matrix:
/// ```text
/// 1 0 0 0
/// 1 1 0 0
/// 1 1 1 0
/// 1 1 1 1
/// ```
fn gen_accumulation_matrix<I: Integer>(dim: usize) -> DenseMatrixZ<I> {
    let mut rows = vec![vec![I::ZERO; dim]; dim];
    for i in 0..dim {
        for j in 0..=i {
            rows[i][j] = I::ONE;
        }
    }
    DenseMatrixZ { rows }
}

/// Multiply each element of the input vector by the lower triangular matrix of the appropriate size.
fn mul_by_acc_matrix<I>(input: &[I]) -> Vec<I>
where
    I: Zero + Add<I> + Clone,
{
    let mut result = vec![I::zero(); input.len()];
    result[0] = input[0].clone();
    for i in 1..input.len() {
        result[i] = result[i - 1].clone() + input[i].clone();
    }
    result
}

impl<ZT: ZipTypes> LinearCode<ZT> for RaaCode<ZT> {
    fn row_len(&self) -> usize {
        self.row_len
    }

    fn codeword_len(&self) -> usize {
        self.row_len * self.repetition_factor
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
        debug_assert_eq!(
            row.len(),
            self.row_len,
            "Row length must match the code's row length"
        );
        let result: Vec<Out> = self.repetition_matrix.mat_vec_mul(row);
        let result: Vec<Out> = self.perm_matrix_1.mat_vec_mul(&result);
        let result: Vec<Out> = mul_by_acc_matrix(&result);
        let result: Vec<Out> = self.perm_matrix_2.mat_vec_mul(&result);
        let result: Vec<Out> = mul_by_acc_matrix(&result);
        debug_assert_eq!(result.len(), self.codeword_len());
        result

        // let alt_result: Vec<Out> = self.encoding_matrix.mat_vec_mul(row);
        // assert_eq!(
        //     result, alt_result,
        //     "Encoding matrices should produce the same result"
        // );
        // debug_assert_eq!(alt_result.len(), self.codeword_len());
        // alt_result
    }

    fn encode_f<F: Field>(&self, row: &[F], field: F::R) -> Vec<F>
    where
        ZT::L: FieldMap<F, Output = F>,
    {
        debug_assert_eq!(
            row.len(),
            self.row_len,
            "Row length must match the code's row length"
        );

        let repetition_matrix = SparseMatrixF::new(&self.repetition_matrix, field);
        let perm_matrix_1 = SparseMatrixF::new(&self.perm_matrix_1, field);
        let perm_matrix_2 = SparseMatrixF::new(&self.perm_matrix_2, field);

        let result: Vec<F> = repetition_matrix.mat_vec_mul(row);
        let result: Vec<F> = perm_matrix_1.mat_vec_mul(&result);
        let result: Vec<F> = mul_by_acc_matrix(&result);
        let result: Vec<F> = perm_matrix_2.mat_vec_mul(&result);
        let result: Vec<F> = mul_by_acc_matrix(&result);
        debug_assert_eq!(result.len(), self.codeword_len());
        result

        // let encoding_matrix: DenseMatrixZ<F> = self.encoding_matrix.map_to_field(field);
        // let alt_result: Vec<F> = encoding_matrix.mat_vec_mul(row);
        // debug_assert_eq!(alt_result.len(), self.codeword_len());
        // assert_eq!(
        //     result, alt_result,
        //     "Encoding matrices should produce the same result"
        // );
        // alt_result
    }
}

/// Matrix content, vector of rows.
#[repr(transparent)]
#[derive(Debug, Clone, Eq, PartialEq)]
struct DenseMatrixZ<L> {
    rows: Vec<Vec<L>>,
}

impl<L: Send + Sync + Clone> DenseMatrixZ<L> {
    fn mat_vec_mul<In, Out>(&self, vec: &[In]) -> Vec<Out>
    where
        In: Send + Sync,
        Out: Zero
            + for<'a> AddAssign<&'a Out>
            + for<'a> From<&'a In>
            + for<'a> From<&'a L>
            + Mul<Output = Out>
            + core::iter::Sum
            + Clone
            + Send
            + Sync,
    {
        self.rows
            .iter()
            .map(|row| {
                let mut sum = Out::zero();
                for (a, b) in row.iter().zip(vec.iter()) {
                    sum += &(Out::from(a) * Out::from(b))
                }
                sum
            })
            .collect()

        // TODO: Enable this parallelized implementation:
        //
        // let mut result = vec![Out::zero(); self.rows.len()];
        // parallelize(&mut result, |(result_row, i)| {
        //     let mut sum = Out::zero();
        //     for (a, b) in self.rows[i].iter().zip(vec.iter()) {
        //         sum += &(Out::from(a) * Out::from(b))
        //     }
        //     result_row[i] = sum;
        // });
        // result
    }

    fn map_to_field<F: Field>(&self, field: F::R) -> DenseMatrixZ<F>
    where
        L: Integer + FieldMap<F, Output = F>,
    {
        DenseMatrixZ {
            rows: self
                .rows
                .iter()
                .map(|row| row.iter().map(|x| x.map_to_field(field)).collect())
                .collect(),
        }
    }

    fn mul_mod_2(&self, rhs: &[Vec<L>]) -> DenseMatrixZ<L>
    where
        L: Integer,
    {
        let mut result = vec![vec![L::ZERO; self.rows[0].len()]; self.rows.len()];
        let modulus = L::ONE + L::ONE; // Modulo 2

        for (i, row) in self.rows.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                for (k, rhs_value) in rhs[j].iter().enumerate() {
                    result[i][k] += &(value.clone() * rhs_value);
                    result[i][k] %= modulus.clone();
                }
            }
        }

        DenseMatrixZ { rows: result }
    }

    fn randomize_sign(&mut self, transcript: &mut impl ZipTranscript<L>)
    where
        L: Integer,
    {
        let mut rng = rand::rngs::StdRng::seed_from_u64(transcript.get_u64());
        for row in &mut self.rows {
            for value in row {
                if !Zero::is_zero(value) && rng.random_bool(0.5) {
                    *value = -value.clone();
                }
            }
        }
    }
}

impl<L: Integer> Mul<Self> for &DenseMatrixZ<L> {
    type Output = DenseMatrixZ<L>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = vec![vec![L::ZERO; self.rows[0].len()]; self.rows.len()];

        for (i, row) in self.rows.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                for (k, rhs_value) in rhs.rows[j].iter().enumerate() {
                    result[i][k] += &(value.clone() * rhs_value);
                }
            }
        }

        DenseMatrixZ { rows: result }
    }
}

impl<L: Integer> Mul<&Self> for DenseMatrixZ<L> {
    type Output = DenseMatrixZ<L>;

    fn mul(self, rhs: &Self) -> Self::Output {
        &self * rhs
    }
}

impl<L: Integer> Mul<&SparseMatrixZ<L>> for &DenseMatrixZ<L> {
    type Output = DenseMatrixZ<L>;

    fn mul(self, rhs: &SparseMatrixZ<L>) -> Self::Output {
        self * &DenseMatrixZ {
            rows: rhs.to_dense(),
        }
    }
}

impl<L: Integer> Mul<&SparseMatrixZ<L>> for DenseMatrixZ<L> {
    type Output = DenseMatrixZ<L>;

    fn mul(self, rhs: &SparseMatrixZ<L>) -> Self::Output {
        &self * rhs
    }
}
