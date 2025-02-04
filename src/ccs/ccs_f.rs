//! The arith module provides utility for operating with arithmetic constraint systems.

#![allow(non_snake_case, dead_code, non_camel_case_types)]

use ark_ff::{One, UniformRand, Zero};
use ark_std::{log2, rand};

use crate::ccs::error::CSError as Error;
use crate::field_config::FieldConfig;
use crate::{biginteger::BigInt, field::RandomField, sparse_matrix::SparseMatrix};

use super::ccs_z::CCS_Z;
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
#[derive(Debug, Clone, PartialEq)]
pub struct CCS_RF<const N: usize> {
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
    pub config: *const FieldConfig<N>,
}

impl<const N: usize> Arith<N> for CCS_RF<N> {
    /// check that a CCS structure is satisfied by a z vector. Only for testing.
    fn check_relation(
        &self,
        M: &[SparseMatrix<RandomField<N>>],
        z: &[RandomField<N>],
    ) -> Result<(), Error> {
        let mut result = vec![RandomField::zero(); self.m];
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

impl<const N: usize> CCS_RF<N> {
    fn pad_rows_to(&mut self, M: &mut [SparseMatrix<RandomField<N>>], size: usize) {
        let size = size.next_power_of_two();
        if size > self.m {
            self.m = size;
            self.s = log2(size) as usize;

            // Update matrices
            M.iter_mut().for_each(|mat| mat.pad_rows(size));
        }
    }
}

/// A representation of a CCS witness.
#[derive(Debug, Clone, PartialEq)]
pub struct Witness<const N: usize> {
    /// `w_ccs` is the original CCS witness.
    pub w_ccs: Vec<RandomField<N>>,
}

impl<const N: usize> Witness<N> {
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
        x: &[SparseMatrix<RandomField<N>>],
        w: &[RandomField<N>],
    ) -> Vec<RandomField<N>>;
}

pub fn from_ccs_z<const N: usize>(
    ccs_z: &CCS_Z,
    config: *const FieldConfig<N>,
) -> Result<CCS_RF<N>, ()> {
    for c in ccs_z.c.iter() {
        let bigint: Result<BigInt<N>, _> = c.magnitude().clone().try_into();
        if bigint.is_err() || bigint.unwrap() >= unsafe { *config }.modulus {
            return Err(());
        }
    }
    // now we can safely convert all the integers into field elements
    let c: Vec<RandomField<N>> = ccs_z
        .c
        .iter()
        .map(|c| match c.sign() {
            num_bigint::Sign::Minus => {
                -RandomField::from_bigint(config, BigInt::try_from(c.magnitude().clone()).unwrap())
                    .unwrap()
            }
            num_bigint::Sign::NoSign => {
                RandomField::from_bigint(config, BigInt::try_from(c.magnitude().clone()).unwrap())
                    .unwrap()
            }
            num_bigint::Sign::Plus => {
                RandomField::from_bigint(config, BigInt::try_from(c.magnitude().clone()).unwrap())
                    .unwrap()
            }
        })
        .collect();
    Ok(CCS_RF {
        m: ccs_z.m,
        n: ccs_z.n,
        l: ccs_z.l,
        t: ccs_z.t,
        q: ccs_z.q,
        d: ccs_z.d,
        s: ccs_z.s,
        s_prime: ccs_z.s_prime,
        S: ccs_z.S.clone(),
        c,
        config,
    })
}
