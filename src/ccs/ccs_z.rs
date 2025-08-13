//! The arith module provides utility for operating with arithmetic constraint systems.

#![allow(non_snake_case, non_camel_case_types)]

use ark_std::{log2, marker::PhantomData, sync::atomic::AtomicPtr, vec, vec::Vec};

use super::{
    ccs_f::{CCS_F, Statement_F, Witness_F},
    utils::{hadamard, mat_vec_mul, vec_add, vec_scalar_mul},
};
use crate::{
    ccs::error::CSError as Error,
    sparse_matrix::{SparseMatrix, dense_matrix_to_sparse},
    traits::{ConfigReference, Field, FieldMap, Integer},
};

///  * `R: Ring` - the ring algebra over which the constraint system operates
pub trait Arith_Z<I: Integer> {
    /// Checks that the given Arith structure is satisfied by a z vector. Used only for testing.
    fn check_relation(&self, M: &[SparseMatrix<I>], z: &[I]) -> Result<(), Error>;

    /// Returns the bytes that represent the parameters, that is, the matrices sizes, the amount of
    /// public inputs, etc, without the matrices/polynomials values.
    fn params_to_le_bytes(&self) -> Vec<u8>;
}

/// CCS represents the Customizable Constraint Systems structure defined in
/// the [CCS paper](https://eprint.iacr.org/2023/552)
#[derive(Debug, Clone, PartialEq)]
pub struct CCS_Z<I> {
    /// m: number of rows in M_i (such that M_i \in F^{m, n})
    pub m: usize,
    /// n = |z|, number of cols in M_i
    pub n: usize,
    /// l = |io|, size of public input/output
    pub l: usize,
    /// t = |M|, number of matrices
    pub t: usize,
    /// q = |c| = |S|, number of multisets
    pub q: usize,
    /// d: max degree in each variable
    pub d: usize,
    /// s = log(m), dimension of x
    pub s: usize,
    /// s_prime = log(n), dimension of y
    pub s_prime: usize,
    /// vector of multisets
    pub S: Vec<Vec<usize>>,
    /// vector of coefficients
    pub c: Vec<i64>,
    pub _phantom: PhantomData<I>,
}

impl<I: Integer> Arith_Z<I> for CCS_Z<I> {
    /// check that a CCS structure is satisfied by a z vector. Only for testing.
    fn check_relation(&self, M: &[SparseMatrix<I>], z: &[I]) -> Result<(), Error> {
        let mut result = vec![I::ZERO; self.m];
        for m in M.iter() {
            assert_eq!(
                m.n_rows, self.m,
                "Incorrect number of rows, expected {} and got {}.",
                self.m, m.n_rows
            );
            assert_eq!(
                m.n_cols, self.n,
                "Incorrect number of rows, expected {} and got {}.",
                self.n, m.n_cols
            );
        }
        for i in 0..self.q {
            // extract the needed M_j matrices out of S_i
            let vec_M_j: Vec<&SparseMatrix<I>> = self.S[i].iter().map(|j| &M[*j]).collect();

            // complete the hadamard chain
            let mut hadamard_result = vec![I::from_i64(1); self.m];
            for M_j in vec_M_j.into_iter() {
                let mut res = mat_vec_mul(M_j, z)?;
                res.resize(self.m, I::ZERO);
                hadamard_result = hadamard(&hadamard_result, &res)?;
            }

            // multiply by the coefficient of this step
            let c_M_j_z = vec_scalar_mul(&hadamard_result, &(I::from(self.c[i])));

            // add it to the final vector
            result = vec_add(&result, &c_M_j_z)?;
        }

        // make sure the final vector is all zeroes
        result
            .iter()
            .all(|item| *item == I::ZERO)
            .then_some(())
            .ok_or(Error::NotSatisfied)
    }

    fn params_to_le_bytes(&self) -> Vec<u8> {
        [
            self.l.to_le_bytes(),
            self.m.to_le_bytes(),
            self.n.to_le_bytes(),
            self.t.to_le_bytes(),
            self.q.to_le_bytes(),
            self.d.to_le_bytes(),
        ]
        .concat()
    }
}

impl<I: Integer> CCS_Z<I> {
    pub fn pad(&mut self, statement: &mut Statement_Z<I>, size: usize) {
        let size = size.next_power_of_two();
        if size > self.m {
            let log_m = log2(size) as usize;
            self.m = size;
            self.s = log_m;
            self.n = size;
            self.s_prime = log_m;

            // Update matrices
            statement
                .constraints
                .iter_mut()
                .for_each(|mat: &mut SparseMatrix<I>| {
                    mat.pad_cols(size);
                    mat.pad_rows(size);
                });
        }
    }
}

impl<F: Field, I: Integer + FieldMap<F, Output = F>> FieldMap<F> for CCS_Z<I> {
    type Output = CCS_F<F>;

    fn map_to_field(&self, config_ref: F::R) -> Self::Output {
        match config_ref.pointer() {
            Some(config_ptr) => CCS_F {
                m: self.m,
                n: self.n,
                l: self.l,
                t: self.t,
                q: self.q,
                d: self.d,
                s: self.s,
                s_prime: self.s_prime,
                S: self.S.clone(),
                c: self.c.iter().map(|c| c.map_to_field(config_ref)).collect(),
                config: AtomicPtr::new(config_ptr),
            },
            None => panic!("FieldConfig cannot be null"),
        }
    }
}

pub struct Statement_Z<I: Integer> {
    pub constraints: Vec<SparseMatrix<I>>,
    pub public_input: Vec<I>,
}

impl<F: Field, I: Integer> FieldMap<F> for Statement_Z<I>
where
    I: FieldMap<F, Output = F>,
{
    type Output = Statement_F<F>;

    fn map_to_field(&self, config_ref: F::R) -> Self::Output {
        Self::Output {
            constraints: self.constraints.map_to_field(config_ref),
            public_input: self.public_input.map_to_field(config_ref),
        }
    }
}
/// A representation of a CCS witness.
#[derive(Debug, Clone, PartialEq)]
pub struct Witness_Z<I> {
    /// `w_ccs` is the original CCS witness.
    pub w_ccs: Vec<I>,
}

impl<I: Integer> Witness_Z<I> {
    /// Create a [`Witness`] from a ccs witness.
    pub fn new(w_ccs: Vec<I>) -> Self {
        Self { w_ccs }
    }
}

// TODO Refactor, it might have lost performance because it's generalized over Int<M>
impl<F: Field, I: Field> FieldMap<F> for Witness_Z<I>
where
    I: FieldMap<F, Output = F>,
{
    type Output = Witness_F<F>;

    fn map_to_field(&self, config_ref: F::R) -> Self::Output {
        Witness_F {
            w_ccs: self.w_ccs.map_to_field(config_ref),
        }
    }
}

/// A trait for defining the behaviour of a satisfying instance of a constraint system
///
/// # Types
///  - `R: Ring` - the ring in which the constraint system is operating.
///
pub trait Instance_Z<I: Integer> {
    /// Given a witness vector, produce a concatonation of the statement and the witness
    fn get_z_vector(&self, w: &[I]) -> Vec<I>;
}

impl<I: Integer> Instance_Z<I> for Statement_Z<I> {
    fn get_z_vector(&self, w: &[I]) -> Vec<I> {
        let mut z: Vec<_> = Vec::with_capacity(self.public_input.len() + w.len() + 1);

        z.extend_from_slice(&self.public_input);
        z.push(I::from_i64(1));
        z.extend_from_slice(w);

        z
    }
}

#[cfg(test)]
pub(crate) fn get_test_ccs_Z<I: Integer>() -> CCS_Z<I> {
    // R1CS for: x^3 + x + 5 = y (example from article
    // https://www.vitalik.ca/general/2016/12/10/qap.html )

    let m = 4;
    let n = 6;
    CCS_Z {
        m,
        n,
        l: 1,
        t: 3,
        q: 2,
        d: 2,
        s: log2(m) as usize,
        s_prime: log2(n) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![1, -1],
        _phantom: PhantomData,
    }
}

pub fn to_Z_matrix<I: Integer>(M: Vec<Vec<usize>>) -> SparseMatrix<I> {
    dense_matrix_to_sparse(
        M.iter()
            .map(|m| m.iter().map(|c| I::from_i64(*c as i64)).collect())
            .collect(),
    )
}

#[cfg(test)]
pub(crate) fn get_test_ccs_Z_statement<I: Integer>(input: i64) -> Statement_Z<I> {
    let A = to_Z_matrix(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 1, 0, 0],
        vec![1, 0, 0, 0, 1, 0],
        vec![0, 5, 0, 0, 0, 1],
    ]);
    let B = to_Z_matrix(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
    ]);
    let C = to_Z_matrix(vec![
        vec![0, 0, 0, 1, 0, 0],
        vec![0, 0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 0, 1],
        vec![0, 0, 1, 0, 0, 0],
    ]);
    let constraints = vec![A, B, C];
    let public_input = vec![I::from_i64(input)];
    Statement_Z {
        constraints,
        public_input,
    }
}

#[cfg(test)]
pub(crate) fn get_test_z_Z<I: Integer>(input: i64) -> Vec<I> {
    // z = (io, 1, w)
    vec![
        I::from_i64(input), // io
        I::from_i64(1),
        I::from_i64(input * input * input + input + 5), // x^3 + x + 5
        I::from_i64(input * input),                     // x^2
        I::from_i64(input * input * input),             // x^2 * x
        I::from_i64(input * input * input + input),     // x^3 + x
    ]
}

#[cfg(test)]
pub(crate) fn get_test_wit_Z<I: Integer>(input: i64) -> Witness_Z<I> {
    Witness_Z::new(vec![
        I::from_i64(input * input * input + input + 5), // x^3 + x + 5
        I::from_i64(input * input),                     // x^2
        I::from_i64(input * input * input),             // x^2 * x
        I::from_i64(input * input * input + input),     // x^3 + x
    ])
}

#[cfg(test)]
pub(crate) fn get_test_ccs_stuff_Z<I: Integer>(
    input: i64,
) -> (CCS_Z<I>, Statement_Z<I>, Witness_Z<I>, Vec<I>) {
    let mut ccs = get_test_ccs_Z();
    let mut statement = get_test_ccs_Z_statement(input);
    let witness = get_test_wit_Z(input);
    let z = get_test_z_Z(input);
    let len = usize::max(ccs.m.next_power_of_two(), ccs.n.next_power_of_two());
    ccs.pad(&mut statement, len);
    (ccs, statement, witness, z)
}

#[cfg(test)]
mod tests {
    use ark_std::vec::Vec;

    use super::{Arith_Z, get_test_ccs_Z, get_test_ccs_Z_statement, get_test_z_Z};
    use crate::{
        big_int,
        ccs::{
            ccs_f::{Arith, CCS_F, Instance_F},
            ccs_z::CCS_Z,
        },
        field::{ConfigRef, Int, RandomField},
        field_config,
        sparse_matrix::SparseMatrix,
        traits::FieldMap,
    };

    #[test]
    fn test_ccs_z() {
        const N: usize = 3;
        let input = 3;
        let ccs: CCS_Z<Int<N>> = get_test_ccs_Z();
        let statement = get_test_ccs_Z_statement(input);
        let z = get_test_z_Z(input);
        let constraints = statement
            .constraints
            .iter()
            .map(|m_i64| {
                let mut m = SparseMatrix::empty();
                m.n_rows = m_i64.n_rows;
                m.n_cols = m_i64.n_cols;
                m.coeffs = m_i64
                    .coeffs
                    .iter()
                    .map(|vec| vec.iter().map(|(coeff, col)| (*coeff, *col)).collect())
                    .collect();
                m
            })
            .collect::<Vec<_>>();

        let res = ccs.check_relation(&constraints, z.as_slice());
        assert!(res.is_ok())
    }

    #[test]
    fn test_ccs_z_conversion() {
        let input = 3;
        let ccs = get_test_ccs_Z::<Int<N>>();
        let statement = get_test_ccs_Z_statement(input);
        let z = get_test_z_Z(input);
        let wit = z[2..].to_vec();

        let constraints = statement
            .constraints
            .iter()
            .map(|m_i64| {
                let mut m = SparseMatrix::empty();
                m.n_rows = m_i64.n_rows;
                m.n_cols = m_i64.n_cols;
                m.coeffs = m_i64
                    .coeffs
                    .iter()
                    .map(|vec| vec.iter().map(|(coeff, col)| (*coeff, *col)).collect())
                    .collect();
                m
            })
            .collect::<Vec<_>>();

        ccs.check_relation(&constraints, &z)
            .expect("Failed to check relation over Integer Ring");

        const N: usize = 3;
        let config = field_config!(312829638388039969874974628075306023441, N);

        let config_ptr = ConfigRef::from(&config);

        let ccs_f: CCS_F<RandomField<N>> = ccs.map_to_field(config_ptr);
        let statement_f = statement.map_to_field(config_ptr);
        let witness_f = wit.map_to_field(config_ptr);
        let z_f = statement_f.get_z_vector(&witness_f, config_ptr);

        ccs_f
            .check_relation(&statement_f.constraints, &z_f)
            .expect("Failed to check relation over Random Field");
    }
}
