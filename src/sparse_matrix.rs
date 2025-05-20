#![allow(non_snake_case)]

use crate::field_config::ConfigRef;
use ark_ff::vec::*;
use ark_ff::Zero;
use ark_std::rand::Rng;
use crypto_bigint::Int;
use crypto_bigint::Random;

use crate::{field::conversion::FieldMap, field::RandomField};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SparseMatrix<R1: Clone + Send + Sync> {
    pub n_rows: usize,
    pub n_cols: usize,
    pub coeffs: Vec<Vec<(R1, usize)>>,
}

impl<R1: Clone + Send + Sync + std::fmt::Display> std::fmt::Display for SparseMatrix<R1> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.n_rows {
            let mut row = vec!["0".to_string(); self.n_cols];
            if i < self.coeffs.len() {
                for (val, col) in &self.coeffs[i] {
                    row[*col] = val.to_string().trim_start_matches('0').to_string();
                    if row[*col].is_empty() {
                        row[*col] = "0".to_string();
                    }
                }
            }
            writeln!(f, "{}", row.join(" "))?;
        }
        Ok(())
    }
}

// At the moment only using i128 for the sparse matrix, macro later if needed
macro_rules! impl_field_map_sparse_matrix {
    ($type:ty) => {
        impl<'cfg, const N: usize> FieldMap<'cfg, N> for SparseMatrix<$type> {
            type Cfg = ConfigRef<'cfg, N>;
            type Output = SparseMatrix<RandomField<'cfg, N>>;
            type Lifetime = ();
            fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
                let mut matrix = SparseMatrix::<RandomField<N>> {
                    n_rows: self.n_rows,
                    n_cols: self.n_cols,
                    coeffs: Vec::new(),
                };
                for row in self.coeffs.iter() {
                    let mut new_row = Vec::new();
                    for (value, col) in row.iter() {
                        new_row.push((value.map_to_field(config), *col));
                    }
                    matrix.coeffs.push(new_row);
                }
                matrix
            }
        }
    };
}
impl<'cfg, const N: usize, const M: usize> FieldMap<'cfg, N> for SparseMatrix<Int<M>> {
    type Cfg = ConfigRef<'cfg, N>;
    type Output = SparseMatrix<RandomField<'cfg, N>>;
    type Lifetime = ();
    fn map_to_field(&self, config: Self::Cfg) -> Self::Output {
        let mut matrix = SparseMatrix::<RandomField<N>> {
            n_rows: self.n_rows,
            n_cols: self.n_cols,
            coeffs: Vec::new(),
        };
        for row in self.coeffs.iter() {
            let mut new_row = Vec::new();
            for (value, col) in row.iter() {
                new_row.push((value.map_to_field(config), *col));
            }
            matrix.coeffs.push(new_row);
        }
        matrix
    }
}
impl_field_map_sparse_matrix!(i64);
impl_field_map_sparse_matrix!(i128);

impl<R: Copy + Send + Sync + Zero + Random> SparseMatrix<R> {
    pub fn empty() -> Self {
        Self {
            n_rows: 0,
            n_cols: 0,
            coeffs: vec![],
        }
    }

    pub fn rand<RND: Rng>(rng: &mut RND, n_rows: usize, n_cols: usize) -> Self {
        const ZERO_VAL_PROBABILITY: f64 = 0.8f64;

        let dense = (0..n_rows)
            .map(|_| {
                (0..n_cols)
                    .map(|_| {
                        if !rng.gen_bool(ZERO_VAL_PROBABILITY) {
                            return R::random(rng);
                        }
                        R::zero()
                    })
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<Vec<R>>>();
        dense_matrix_to_sparse(dense)
    }

    pub fn to_dense(&self) -> Vec<Vec<R>> {
        let mut r: Vec<Vec<R>> = vec![vec![R::zero(); self.n_cols]; self.n_rows];
        for (row_i, row) in self.coeffs.iter().enumerate() {
            for &(value, col_i) in row.iter() {
                r[row_i][col_i] = value;
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

pub fn dense_matrix_to_sparse<R: Copy + Send + Sync + Zero>(m: Vec<Vec<R>>) -> SparseMatrix<R> {
    let mut r = SparseMatrix::<R> {
        n_rows: m.len(),
        n_cols: m[0].len(),
        coeffs: Vec::new(),
    };
    for m_row in m.iter() {
        let mut row: Vec<(R, usize)> = Vec::new();
        for (col_i, value) in m_row.iter().enumerate() {
            if !value.is_zero() {
                row.push((*value, col_i));
            }
        }
        r.coeffs.push(row);
    }
    r
}

pub fn dense_matrix_u64_to_sparse<R>(m: Vec<Vec<u64>>) -> SparseMatrix<R>
where
    R: Copy
        + Send
        + Sync
        + Zero
        + From<u64>
        + ark_serialize::Valid
        + ark_serialize::CanonicalSerialize
        + ark_serialize::CanonicalDeserialize,
{
    let mut r = SparseMatrix::<R> {
        n_rows: m.len(),
        n_cols: m[0].len(),
        coeffs: Vec::new(),
    };
    for m_row in m.iter() {
        let mut row: Vec<(R, usize)> = Vec::new();
        for (col_i, &value) in m_row.iter().enumerate() {
            if value != 0 {
                row.push((R::from(value), col_i));
            }
        }
        r.coeffs.push(row);
    }
    r
}

pub fn compute_eval_table_sparse<'cfg, const N: usize>(
    M: &SparseMatrix<RandomField<'cfg, N>>,
    rx: &[RandomField<'cfg, N>],
    num_rows: usize,
    num_cols: usize,
    config: ConfigRef<'cfg, N>,
) -> Vec<RandomField<'cfg, N>> {
    assert_eq!(rx.len(), num_rows);
    M.coeffs.iter().enumerate().fold(
        vec![0u32.map_to_field(config); num_cols],
        |mut M_evals, (row, vals)| {
            for (val, col) in vals {
                M_evals[*col] += &(rx[row] * val);
            }
            M_evals
        },
    )
}
