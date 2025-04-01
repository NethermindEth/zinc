//! The arith module provides utility for operating with arithmetic constraint systems.

#![allow(non_snake_case, dead_code, non_camel_case_types)]

use std::sync::atomic::AtomicPtr;

use ark_ff::{One, UniformRand, Zero};
use ark_std::{log2, rand};

use crate::ccs::error::CSError as Error;
use crate::field::conversion::FieldMap;
use crate::field_config::FieldConfig;
use crate::poly_f::mle::{DenseMultilinearExtension, SparseMultilinearExtension};
use crate::sparse_matrix::{compute_eval_table_sparse, dense_matrix_to_sparse};
use crate::{field::RandomField, sparse_matrix::SparseMatrix};

use super::utils::{hadamard, mat_vec_mul, vec_add, vec_scalar_mul};

/// A trait for defining the behaviour of an arithmetic constraint system.
///
/// ## Type Parameters
///
///  * `R: Ring` - the ring algebra over which the constraint system operates
pub trait Arith<const N: usize> {
    /// Checks that the given Arith structure is satisfied by a z vector. Used only for testing.
    fn check_relation(
        &self,
        M: &[SparseMatrix<RandomField<N>>],
        z: &[RandomField<N>],
    ) -> Result<(), Error>;

    /// Returns the bytes that represent the parameters, that is, the matrices sizes, the amount of
    /// public inputs, etc, without the matrices/polynomials values.
    fn params_to_le_bytes(&self) -> Vec<u8>;
}

/// CCS represents the Customizable Constraint Systems structure defined in
/// the [CCS paper](https://eprint.iacr.org/2023/552)
#[derive(Debug)]
pub struct CCS_F<const N: usize> {
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
    pub c: Vec<RandomField<N>>,
    /// The field the constraint system operates in
    pub config: AtomicPtr<FieldConfig<N>>,
}

impl<const N: usize> Arith<N> for CCS_F<N> {
    /// check that a CCS structure is satisfied by a z vector. Only for testing.
    fn check_relation(
        &self,
        M: &[SparseMatrix<RandomField<N>>],
        z: &[RandomField<N>],
    ) -> Result<(), Error> {
        let mut result = vec![RandomField::zero(); self.m];
        for m in M.iter() {
            assert_eq!(
                m.n_rows, self.n,
                "Incorrect number of rows, expected {} and got {}.",
                self.m, m.n_rows
            );
            assert_eq!(
                m.n_cols, self.m,
                "Incorrect number of cols, expected {} and got {}.",
                self.n, m.n_cols
            );
        }
        for i in 0..self.q {
            // extract the needed M_j matrices out of S_i
            let vec_M_j: Vec<&SparseMatrix<RandomField<N>>> =
                self.S[i].iter().map(|j| &M[*j]).collect();

            // complete the hadamard chain
            let mut hadamard_result = vec![RandomField::one(); self.m];
            for M_j in vec_M_j.into_iter() {
                let mut res = mat_vec_mul(M_j, z)?;
                res.resize(self.m, RandomField::zero());
                hadamard_result = hadamard(&hadamard_result, &res)?;
            }

            // multiply by the coefficient of this step
            let c_M_j_z = vec_scalar_mul(&hadamard_result, &self.c[i]);

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

impl<const N: usize> CCS_F<N> {
    fn pad(&mut self, statement: &mut Statement_F<N>, size: usize) {
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
                .for_each(|mat: &mut SparseMatrix<RandomField<N>>| {
                    mat.pad_cols(size);
                    mat.pad_rows(size);
                });
        }
    }
}

/// A representation of a CCS statement
#[derive(Debug, Clone, PartialEq)]
pub struct Statement_F<const N: usize> {
    pub constraints: Vec<SparseMatrix<RandomField<N>>>,
    pub public_input: Vec<RandomField<N>>,
}

impl<const N: usize> Statement_F<N> {
    pub fn compute_eval_table_sparse(
        &self,
        num_rows: usize,
        num_cols: usize,
        ccs: &CCS_F<N>,
        evals: &[RandomField<N>],
    ) -> Vec<Vec<RandomField<N>>> {
        assert_eq!(num_rows, ccs.n);
        assert!(num_cols > (ccs.m - ccs.l) - 1);

        self.constraints
            .iter()
            .map(|M| {
                compute_eval_table_sparse(M, evals, num_rows, num_cols, unsafe {
                    *(ccs.config.as_ptr())
                })
            })
            .collect()
    }
}

/// A representation of a linearised CCS statement
#[derive(Debug, Clone, PartialEq)]
pub struct LStatement<const N: usize> {
    constraints: Vec<SparseMultilinearExtension<N>>,
    r: Vec<RandomField<N>>,
}

/// A representation of a CCS witness.
#[derive(Debug, Clone, PartialEq)]
pub struct Witness_F<const N: usize> {
    /// `w_ccs` is the original CCS witness.
    pub w_ccs: Vec<RandomField<N>>,
}

/// A representation of a linearised CCS witness.
#[derive(Debug, Clone, PartialEq)]
pub struct LWitness<const N: usize> {
    /// `w_ccs` is the original CCS witness.
    pub lw_ccs: DenseMultilinearExtension<N>,
}

impl<const N: usize> Witness_F<N> {
    /// Create a [`Witness`] from a ccs witness.
    pub fn new(w_ccs: Vec<RandomField<N>>) -> Self {
        Self { w_ccs }
    }

    /// Generates a random witness by firstly generating a random
    /// vector of arbitrary norm and then computing the rest of the data
    /// needed for a witness.
    ///
    /// # Arguments
    /// * `rng` is a mutable reference to the random number generator.
    /// * `w_ccs_len` is the length of the non-decomposed witness (a.k.a. the CCS witness).
    pub fn rand<Rng: rand::Rng + ?Sized>(rng: &mut Rng, w_ccs_len: usize) -> Self {
        Self::new(
            (0..w_ccs_len)
                .map(|_| RandomField::<N>::rand(rng))
                .collect(),
        )
    }
}

/// A trait for defining the behaviour of a satisfying instance of a constraint system
///
/// # Types
///  - `R: Ring` - the ring in which the constraint system is operating.
///
pub trait Instance_F<const N: usize> {
    /// Given a witness vector, produce a concatonation of the statement and the witness
    fn get_z_vector(
        &self,
        w: &[RandomField<N>],
        config: *const FieldConfig<N>,
    ) -> Vec<RandomField<N>>;
}

impl<const N: usize> Instance_F<N> for Statement_F<N> {
    fn get_z_vector(
        &self,
        w: &[RandomField<N>],
        config: *const FieldConfig<N>,
    ) -> Vec<RandomField<N>> {
        let mut z: Vec<RandomField<N>> = Vec::with_capacity(self.public_input.len() + w.len() + 1);

        z.extend_from_slice(&self.public_input);
        z.push(1u32.map_to_field(config));
        z.extend_from_slice(w);

        z
    }
}

/// Returns a sparse matrix of field elements given a matrix of unsigned ints
pub fn to_F_matrix<const N: usize>(
    config: *const FieldConfig<N>,
    M: Vec<Vec<usize>>,
) -> SparseMatrix<RandomField<N>> {
    dense_matrix_to_sparse(
        M.iter()
            .map(|m| m.iter().map(|c| (*c as u64).map_to_field(config)).collect())
            .collect(),
    )
}

/// Returns a dense matrix of field elements given a matrix of unsigned ints
pub fn to_F_dense_matrix<const N: usize>(
    config: *const FieldConfig<N>,
    M: Vec<Vec<usize>>,
) -> Vec<Vec<RandomField<N>>> {
    M.iter()
        .map(|m| m.iter().map(|c| (*c as u64).map_to_field(config)).collect())
        .collect()
}

/// Returns a vector of field elements given a vector of unsigned ints
pub fn to_F_vec<const N: usize>(z: Vec<u64>, config: *const FieldConfig<N>) -> Vec<RandomField<N>> {
    z.iter().map(|c| (*c).map_to_field(config)).collect()
}

#[cfg(test)]
pub(crate) fn get_test_ccs_F<const N: usize>(config: *const FieldConfig<N>) -> CCS_F<N> {
    use std::ops::Neg;

    use crate::field::conversion::FieldMap;
    // R1CS for: x^3 + x + 5 = y (example from article
    // https://www.vitalik.ca/general/2016/12/10/qap.html )

    let m = 6;
    let n = 4;
    CCS_F {
        m,
        n,
        l: 1,
        t: 3,
        q: 2,
        d: 2,
        s: log2(m) as usize,
        s_prime: log2(n) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![1u32.map_to_field(config), (1u32.map_to_field(config)).neg()],
        config: AtomicPtr::new(config as *mut FieldConfig<N>),
    }
}

#[cfg(test)]
pub(crate) fn get_test_ccs_F_statement<const N: usize>(
    input: u64,
    config: *const FieldConfig<N>,
) -> Statement_F<N> {
    let A = to_F_matrix(
        config,
        vec![
            vec![1, 0, 0, 0, 0, 0],
            vec![0, 0, 0, 1, 0, 0],
            vec![1, 0, 0, 0, 1, 0],
            vec![0, 5, 0, 0, 0, 1],
        ],
    );
    let B = to_F_matrix(
        config,
        vec![
            vec![1, 0, 0, 0, 0, 0],
            vec![1, 0, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0, 0],
        ],
    );
    let C = to_F_matrix(
        config,
        vec![
            vec![0, 0, 0, 1, 0, 0],
            vec![0, 0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 0, 1],
            vec![0, 0, 1, 0, 0, 0],
        ],
    );
    let constraints = vec![A, B, C];
    let public_input = to_F_vec(vec![input], config);
    Statement_F {
        constraints,
        public_input,
    }
}

#[cfg(test)]
pub(crate) fn get_test_z_F<const N: usize>(
    input: u64,
    config: *const FieldConfig<N>,
) -> Vec<RandomField<N>> {
    // z = (io, 1, w)
    to_F_vec(
        vec![
            input, // io
            1,
            input * input * input + input + 5, // x^3 + x + 5
            input * input,                     // x^2
            input * input * input,             // x^2 * x
            input * input * input + input,     // x^3 + x
        ],
        config,
    )
}

#[cfg(test)]
pub(crate) fn get_test_wit_F<const N: usize>(
    input: u64,
    config: *const FieldConfig<N>,
) -> Witness_F<N> {
    Witness_F::new(to_F_vec(
        vec![
            input * input * input + input + 5, // x^3 + x + 5
            input * input,                     // x^2
            input * input * input,             // x^2 * x
            input * input * input + input,     // x^3 + x
        ],
        config,
    ))
}

#[cfg(test)]
pub(crate) fn get_test_ccs_stuff_F<const N: usize>(
    input: u64,
    config: *const FieldConfig<N>,
) -> (CCS_F<N>, Statement_F<N>, Witness_F<N>, Vec<RandomField<N>>) {
    let mut ccs = get_test_ccs_F::<N>(config);
    let mut statement = get_test_ccs_F_statement::<N>(input, config);
    let witness = get_test_wit_F(input, config);
    let z = get_test_z_F(input, config);
    let len = usize::max(ccs.m.next_power_of_two(), ccs.n.next_power_of_two());
    ccs.pad(&mut statement, len);
    (ccs, statement, witness, z)
}

#[cfg(test)]
mod tests {
    use crate::{
        biginteger::BigInt, ccs::test_utils::get_dummy_ccs_F_from_z_length,
        field_config::FieldConfig,
    };

    use super::{get_test_ccs_F, get_test_ccs_F_statement, get_test_z_F, Arith};
    #[test]
    fn test_ccs_f() {
        use std::str::FromStr;

        const N: usize = 2;
        let config: *const FieldConfig<N> =
            &FieldConfig::new(BigInt::from_str("75671012754143952277701807739").unwrap());
        let input = 3;
        let ccs = get_test_ccs_F::<N>(config);
        let statement = get_test_ccs_F_statement::<N>(input, config);
        let z = get_test_z_F::<N>(input, config);

        let res = ccs.check_relation(&statement.constraints, &z);
        assert!(res.is_ok())
    }

    #[test]
    fn test_dummy_ccs_f() {
        use std::str::FromStr;

        const N: usize = 2;
        let config: *const FieldConfig<N> =
            &FieldConfig::new(BigInt::from_str("75671012754143952277701807739").unwrap());
        let mut rng = ark_std::test_rng();
        let n = 1 << 13;
        let (z, ccs, statement, _) = get_dummy_ccs_F_from_z_length::<N>(n, &mut rng, config);

        let res = ccs.check_relation(&statement.constraints, &z);
        assert!(res.is_ok())
    }
}
