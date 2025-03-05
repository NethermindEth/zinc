//! The arith module provides utility for operating with arithmetic constraint systems.

#![allow(non_snake_case, dead_code, non_camel_case_types)]

use std::sync::atomic::AtomicPtr;

use ark_ff::{One, Zero};

use super::ccs_f::{Statement_F, Witness_F, CCS_F};
use super::utils::{hadamard, mat_vec_mul, vec_add, vec_scalar_mul};
use crate::ccs::error::CSError as Error;
use crate::field_config::FieldConfig;
use crate::sparse_matrix::SparseMatrix;
use crate::field::conversion::FieldMap;
use num_bigint::BigInt as Z;

/// A trait for defining the behaviour of an arithmetic constraint system.
///
/// ## Type Parameters
///
///  * `R: Ring` - the ring algebra over which the constraint system operates
pub trait Arith_Z {
    /// Checks that the given Arith structure is satisfied by a z vector. Used only for testing.
    fn check_relation(&self, M: &[SparseMatrix<Z>], z: &[Z]) -> Result<(), Error>;

    /// Returns the bytes that represent the parameters, that is, the matrices sizes, the amount of
    /// public inputs, etc, without the matrices/polynomials values.
    fn params_to_le_bytes(&self) -> Vec<u8>;
}

/// CCS represents the Customizable Constraint Systems structure defined in
/// the [CCS paper](https://eprint.iacr.org/2023/552)
#[derive(Debug, Clone, PartialEq)]
pub struct CCS_Z {
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
    pub c: Vec<i128>,
}

impl Arith_Z for CCS_Z {
    /// check that a CCS structure is satisfied by a z vector. Only for testing.
    fn check_relation(&self, M: &[SparseMatrix<Z>], z: &[Z]) -> Result<(), Error> {
        let mut result = vec![Z::zero(); self.m];
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
            let vec_M_j: Vec<&SparseMatrix<Z>> = self.S[i].iter().map(|j| &M[*j]).collect();

            // complete the hadamard chain
            let mut hadamard_result = vec![Z::one(); self.m];
            for M_j in vec_M_j.into_iter() {
                let mut res = mat_vec_mul(M_j, z)?;
                res.resize(self.m, Z::zero());
                hadamard_result = hadamard(&hadamard_result, &res)?;
            }

            // multiply by the coefficient of this step
            let c_M_j_z = vec_scalar_mul(&hadamard_result, &self.c[i].into());

            // add it to the final vector
            result = vec_add(&result, &c_M_j_z)?;
        }

        // make sure the final vector is all zeroes
        result
            .iter()
            .all(|item| item.is_zero())
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

impl FieldMap for CCS_Z {
    type Output<const N: usize> = CCS_F<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        CCS_F{
            m: self.m,
            n: self.n,
            l: self.l,
            t: self.t,
            q: self.q,
            d: self.d,
            s: self.s,
            s_prime: self.s_prime,
            S: self.S.clone(),
            c: self.c.iter().map(|c| c.map_to_field(config)).collect(),
            config: AtomicPtr::new(config as *mut FieldConfig<N>),
        }
    }
}

pub struct Statement_Z {
    pub constraints: Vec<SparseMatrix<i128>>,
    pub public_input: Vec<i64>,
}

impl FieldMap for Statement_Z {
    type Output<const N: usize> = Statement_F<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        Statement_F{
            constraints: self.constraints.iter().map(|m| m.map_to_field(config)).collect(),
            public_input: self.public_input.iter().map(|i| i.map_to_field(config)).collect(),
        }
    }
}
/// A representation of a CCS witness.
#[derive(Debug, Clone, PartialEq)]
pub struct Witness_Z {
    /// `w_ccs` is the original CCS witness.
    pub w_ccs: Vec<i64>,
}

impl Witness_Z {
    /// Create a [`Witness`] from a ccs witness.
    pub fn new(w_ccs: Vec<i64>) -> Self {
        Self { w_ccs }
    }
}

impl FieldMap for Witness_Z {
    type Output<const N: usize> = Witness_F<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        Witness_F { w_ccs: self.w_ccs.iter().map(|i| i.map_to_field(config)).collect() }
    }
}

/// A trait for defining the behaviour of a satisfying instance of a constraint system
///
/// # Types
///  - `R: Ring` - the ring in which the constraint system is operating.
///
pub trait Instance_Z {
    /// Given a witness vector, produce a concatonation of the statement and the witness
    fn get_z_vector(&self, x: &[i64], w: &[i64]) -> Vec<Z>;
}
