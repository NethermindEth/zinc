//! The arith module provides utility for operating with arithmetic constraint systems.

#![allow(non_snake_case, dead_code, non_camel_case_types)]

use std::sync::atomic::AtomicPtr;

use ark_std::log2;

use super::ccs_f::{Statement_F, Witness_F, CCS_F};
use super::utils::{hadamard, mat_vec_mul, vec_add, vec_scalar_mul};
use crate::ccs::error::CSError as Error;
use crate::field::conversion::FieldMap;
use crate::field_config::FieldConfig;
use crate::sparse_matrix::{dense_matrix_to_sparse, SparseMatrix};

/// A trait for defining the behaviour of an arithmetic constraint system.
///
/// ## Type Parameters
///
///  * `R: Ring` - the ring algebra over which the constraint system operates
pub trait Arith_Z {
    /// Checks that the given Arith structure is satisfied by a z vector. Used only for testing.
    fn check_relation(&self, M: &[SparseMatrix<i64>], z: &[i64]) -> Result<(), Error>;

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
    pub c: Vec<i64>,
}

impl Arith_Z for CCS_Z {
    /// check that a CCS structure is satisfied by a z vector. Only for testing.
    fn check_relation(&self, M: &[SparseMatrix<i64>], z: &[i64]) -> Result<(), Error> {
        let mut result = vec![0i64; self.m];
        for m in M.iter() {
            assert_eq!(
                m.n_rows, self.n,
                "Incorrect number of rows, expected {} and got {}.",
                self.m, m.n_rows
            );
            assert_eq!(
                m.n_cols, self.m,
                "Incorrect number of rows, expected {} and got {}.",
                self.n, m.n_cols
            );
        }
        for i in 0..self.q {
            // extract the needed M_j matrices out of S_i
            let vec_M_j: Vec<&SparseMatrix<i64>> = self.S[i].iter().map(|j| &M[*j]).collect();

            // complete the hadamard chain
            let mut hadamard_result = vec![1i64; self.m];
            for M_j in vec_M_j.into_iter() {
                let mut res = mat_vec_mul(M_j, z)?;
                res.resize(self.m, 0i64);
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
            .all(|item| *item == 0)
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

impl CCS_Z {
    pub fn pad(&mut self, statement: &mut Statement_Z, size: usize) {
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
                .for_each(|mat: &mut SparseMatrix<i64>| {
                    mat.pad_cols(size);
                    mat.pad_rows(size);
                });
        }
    }
}

impl FieldMap for CCS_Z {
    type Output<const N: usize> = CCS_F<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        CCS_F {
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
    pub constraints: Vec<SparseMatrix<i64>>,
    pub public_input: Vec<i64>,
}

impl FieldMap for Statement_Z {
    type Output<const N: usize> = Statement_F<N>;
    fn map_to_field<const N: usize>(&self, config: *const FieldConfig<N>) -> Self::Output<N> {
        Statement_F {
            constraints: self
                .constraints
                .iter()
                .map(|m| m.map_to_field(config))
                .collect(),
            public_input: self
                .public_input
                .iter()
                .map(|i| i.map_to_field(config))
                .collect(),
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
        Witness_F {
            w_ccs: self.w_ccs.iter().map(|i| i.map_to_field(config)).collect(),
        }
    }
}

/// A trait for defining the behaviour of a satisfying instance of a constraint system
///
/// # Types
///  - `R: Ring` - the ring in which the constraint system is operating.
///
pub trait Instance_Z {
    /// Given a witness vector, produce a concatonation of the statement and the witness
    fn get_z_vector(&self, x: &[i64], w: &[i64]) -> Vec<i128>;
}

#[cfg(test)]
pub(crate) fn get_test_ccs_Z() -> CCS_Z {
    // R1CS for: x^3 + x + 5 = y (example from article
    // https://www.vitalik.ca/general/2016/12/10/qap.html )

    let m = 6;
    let n = 4;
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
    }
}

pub fn to_Z_matrix(M: Vec<Vec<usize>>) -> SparseMatrix<i64> {
    dense_matrix_to_sparse(
        M.iter()
            .map(|m| m.iter().map(|c| *c as i64).collect())
            .collect(),
    )
}

#[cfg(test)]
pub(crate) fn get_test_ccs_Z_statement(input: i64) -> Statement_Z {
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
    let public_input = vec![input];
    Statement_Z {
        constraints,
        public_input,
    }
}

#[cfg(test)]
pub(crate) fn get_test_z_Z(input: i64) -> Vec<i64> {
    // z = (io, 1, w)
    vec![
        input, // io
        1,
        input * input * input + input + 5, // x^3 + x + 5
        input * input,                     // x^2
        input * input * input,             // x^2 * x
        input * input * input + input,     // x^3 + x
    ]
}

#[cfg(test)]
pub(crate) fn get_test_wit_Z(input: i64) -> Witness_Z {
    Witness_Z::new(vec![
        input * input * input + input + 5, // x^3 + x + 5
        input * input,                     // x^2
        input * input * input,             // x^2 * x
        input * input * input + input,     // x^3 + x
    ])
}

#[cfg(test)]
pub(crate) fn get_test_ccs_stuff_Z(input: i64) -> (CCS_Z, Statement_Z, Witness_Z, Vec<i64>) {
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

    use crate::ccs::test_utils::get_dummy_ccs_Z_from_z_length;

    use super::{get_test_ccs_Z, get_test_ccs_Z_statement, get_test_z_Z, Arith_Z};

    #[test]
    fn test_ccs_z() {
        let input = 3;
        let ccs = get_test_ccs_Z();
        let statement = get_test_ccs_Z_statement(input);
        let z = get_test_z_Z(input)
            .into_iter()
            .map(|i| i as i64)
            .collect::<Vec<_>>();

        let res = ccs.check_relation(statement.constraints.as_slice(), z.as_slice());
        assert!(res.is_ok())
    }

    // #[test]
    // fn test_dummy_ccs_z() {
    //     let mut rng = ark_std::test_rng();
    //     let n = 1 << 13;
    //     let (z, ccs, statement, _) = get_dummy_ccs_Z_from_z_length(n, &mut rng);

    //     let res = ccs.check_relation(&statement.constraints, &z);
    //     assert!(res.is_ok())
    // }
}
