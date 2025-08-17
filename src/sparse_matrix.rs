use crypto_bigint::Int;

use crate::traits::{Field, Ring};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SparseMatrix<F: Clone + Send + Sync> {
    pub n_rows: usize,
    pub n_cols: usize,
    pub coeffs: Vec<Vec<(F, usize)>>,
}

// At the moment only using i128 for the sparse matrix, macro later if needed

impl<F: Ring + for<'a> From<&'a T>, T: Ring> From<&SparseMatrix<T>> for SparseMatrix<F> {
    fn from(value: &SparseMatrix<T>) -> Self {
        let mut matrix = SparseMatrix::<F> {
            n_rows: value.n_rows,
            n_cols: value.n_cols,
            coeffs: Vec::with_capacity(value.coeffs.len()),
        };

        for row in value.coeffs.iter() {
            let mut new_row = Vec::with_capacity(row.len());
            for (value, col) in row.iter() {
                new_row.push((F::from(value), *col));
            }
            matrix.coeffs.push(new_row);
        }

        matrix
    }
}

impl<F: Ring> SparseMatrix<F> {
    pub fn empty() -> Self {
        Self {
            n_rows: 0,
            n_cols: 0,
            coeffs: vec![],
        }
    }

    pub fn to_dense(&self) -> Vec<Vec<F>> {
        let mut r: Vec<Vec<F>> = vec![vec![F::zero(); self.n_cols]; self.n_rows];
        for (row_i, row) in self.coeffs.iter().enumerate() {
            for (value, col_i) in row.iter() {
                r[row_i][*col_i] = *value;
            }
        }
        r
    }

    pub fn nrows(&self) -> usize {
        self.n_rows
    }

    pub fn ncols(&self) -> usize {
        self.n_cols
    }

    pub fn pad_rows(&mut self, new_size: usize) {
        if new_size > self.nrows() {
            self.n_rows = new_size;
        }
    }

    pub fn pad_cols(&mut self, new_size: usize) {
        if new_size > self.ncols() {
            self.n_cols = new_size;
        }
    }
}

impl<F: Ring + for<'a> From<&'a T>, T: Ring> From<Vec<Vec<T>>> for SparseMatrix<F> {
    fn from(value: Vec<Vec<T>>) -> Self {
        let mut matrix = SparseMatrix {
            n_rows: value.len(),
            n_cols: value[0].len(),
            coeffs: Vec::with_capacity(value.len()),
        };

        for m_row in value.iter() {
            let mut row: Vec<(F, usize)> = Vec::new();
            for (col_i, value) in m_row.iter().enumerate() {
                if !value.is_zero() {
                    row.push((F::from(value), col_i));
                }
            }
            matrix.coeffs.push(row);
        }

        matrix
    }
}

impl<const N: usize> SparseMatrix<Int<N>> {
    pub(crate) fn map_to_field<F: Field<LIMBS>, const LIMBS: usize>(&self) -> SparseMatrix<F> {
        let mut matrix = SparseMatrix::<F> {
            n_rows: self.n_rows,
            n_cols: self.n_cols,
            coeffs: Vec::new(),
        };
        for row in self.coeffs.iter() {
            let mut new_row = Vec::new();
            for (value, col) in row.iter() {
                new_row.push((F::from(value.resize()), *col));
            }
            matrix.coeffs.push(new_row);
        }
        matrix
    }
}

#[allow(non_snake_case)]
pub fn compute_eval_table_sparse<F: Ring>(
    M: &SparseMatrix<F>,
    rx: &[F],
    num_rows: usize,
    num_cols: usize,
) -> Vec<F> {
    assert_eq!(rx.len(), num_rows);
    M.coeffs
        .iter()
        .enumerate()
        .fold(vec![F::zero(); num_cols], |mut M_evals, (row, vals)| {
            for (val, col) in vals {
                M_evals[*col] += rx[row] * val;
            }
            M_evals
        })
}

#[allow(non_snake_case)]
pub fn to_Z_matrix<const I: usize>(m: Vec<Vec<usize>>) -> SparseMatrix<Int<I>> {
    let m: Vec<Vec<Int<I>>> = m
        .iter()
        .map(|row| row.iter().map(|v| Int::from(*v as i128)).collect())
        .collect();
    SparseMatrix::from(m)
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{U128, const_monty_params};
    use num_traits::{One, Zero};

    use super::*;
    use crate::field::RandomField;

    const_monty_params!(ModP, U128, "0076F668F4274572E39A3EA8285319B5");
    type Fp = RandomField<ModP, { U128::LIMBS }>;

    fn mat(values: &[[u64; 3]; 2]) -> SparseMatrix<Fp> {
        let dense: Vec<Vec<Fp>> = values
            .iter()
            .map(|row| row.iter().cloned().map(Fp::from).collect())
            .collect();
        SparseMatrix::from(dense)
    }

    #[test]
    fn test_from_vecvec_and_to_dense_roundtrip() {
        let dense = vec![
            vec![Fp::from(0u64), Fp::from(2u64), Fp::from(0u64)],
            vec![Fp::from(5u64), Fp::from(0u64), Fp::from(7u64)],
        ];
        let sm: SparseMatrix<Fp> = SparseMatrix::from(dense.clone());
        assert_eq!(sm.n_rows, 2);
        assert_eq!(sm.n_cols, 3);
        assert_eq!(sm.coeffs[0].len(), 1);
        assert_eq!(sm.coeffs[1].len(), 2);
        assert_eq!(sm.to_dense(), dense);
    }

    #[test]
    fn test_from_sparsematrix_conversion_between_fields() {
        let a = mat(&[[0, 2, 0], [5, 0, 7]]);
        let b: SparseMatrix<Fp> = SparseMatrix::from(&a);
        assert_eq!(a.to_dense(), b.to_dense());
        assert_eq!(a.n_rows, b.n_rows);
        assert_eq!(a.n_cols, b.n_cols);
    }

    #[test]
    fn test_pad_rows_and_cols_only_expand() {
        let mut sm = mat(&[[1, 0, 0], [0, 0, 2]]);
        let orig_rows = sm.n_rows;
        let orig_cols = sm.n_cols;
        sm.pad_rows(orig_rows - 1);
        sm.pad_cols(orig_cols - 1);
        assert_eq!(sm.n_rows, orig_rows);
        assert_eq!(sm.n_cols, orig_cols);
        sm.pad_rows(orig_rows + 3);
        sm.pad_cols(orig_cols + 4);
        assert_eq!(sm.n_rows, orig_rows + 3);
        assert_eq!(sm.n_cols, orig_cols + 4);
        let dense = sm.to_dense();
        assert_eq!(dense.len(), sm.n_rows);
        assert_eq!(dense[0].len(), sm.n_cols);
        assert_eq!(dense[0][0], Fp::one());
        assert_eq!(dense[1][2], Fp::from(2u64));
    }

    #[test]
    fn test_compute_eval_table_sparse_matches_manual() {
        let dense = vec![
            vec![Fp::from(0u64), Fp::from(3u64), Fp::zero(), Fp::from(4u64)],
            vec![Fp::from(5u64), Fp::zero(), Fp::from(6u64), Fp::zero()],
            vec![Fp::zero(), Fp::from(2u64), Fp::from(1u64), Fp::zero()],
        ];
        let sm: SparseMatrix<Fp> = SparseMatrix::from(dense.clone());
        let rx = vec![Fp::from(10u64), Fp::from(20u64), Fp::from(30u64)];
        let evals = compute_eval_table_sparse(&sm, &rx, sm.n_rows, sm.n_cols);
        let mut expected = vec![Fp::zero(); sm.n_cols];
        for i in 0..sm.n_rows {
            for (j, value) in expected.iter_mut().enumerate().take(sm.n_cols) {
                *value += rx[i] * dense[i][j];
            }
        }
        assert_eq!(evals, expected);
    }
}
