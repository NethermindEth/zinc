use std::{
    marker::PhantomData,
    ops::{Add, Mul},
};

use ark_std::log2;
use crypto_bigint::{Int, Random};
use num_traits::One;
use rand_core::RngCore;
use thiserror::Error;

use crate::{
    poly::{dense::DenseMultilinearExtension, sparse::SparseMultilinearExtension},
    sparse_matrix::{SparseMatrix, compute_eval_table_sparse},
    traits::{Field, Ring},
};

#[derive(Debug, Error)]
pub enum CSError {
    /// The provided witness does not satisfy the constraint system
    #[error("constraint system is not satisfied")]
    NotSatisfied,

    /// The constraint system is not of length $2^k$ for any $k \in \mathbb{N}$.
    ///
    /// More to the point, the witness length will not be a power of 2,
    /// so we cannot use it as a MLE.
    #[error("constraint system matrices rows length (m) not a power of 2: {0}")]
    MatricesRowsLengthNotPowerOf2(usize),

    /// This error occurs if the CCS instance is not correctly padded
    ///
    /// See [definition 4.3](https://eprint.iacr.org/2024/257.pdf#page=40) of LatticeFold paper.
    #[error("constraint system has invalid size bounds: m = {0}, n = {1}, L = {2}")]
    InvalidSizeBounds(usize, usize, usize),

    /// This error occurs when performing operations on vectors of differing lengths.
    #[error("vectors {0} and {1} have different lengths: {0} and {1}")]
    LengthsNotEqual(String, String, usize, usize),
}

pub trait Arith<F: Ring> {
    /// Checks that the given Arith structure is satisfied by a z vector. Used only for testing.
    #[allow(non_snake_case)]
    fn check_relation(&self, M: &[SparseMatrix<F>], z: &[F]) -> Result<(), CSError>;

    /// Returns the bytes that represent the parameters, that is, the matrices sizes, the amount of
    /// public inputs, etc, without the matrices/polynomials values.
    fn params_to_le_bytes(&self) -> Vec<u8>;
}

// CCS represents the Customizable Constraint Systems structure defined in
/// the [CCS paper](https://eprint.iacr.org/2023/552)
#[allow(non_snake_case)]
#[derive(Debug, Clone, PartialEq)]
pub struct CcsF<F: Ring> {
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
    pub c: Vec<F>,
}

/// A representation of a CCS statement
#[derive(Debug, Clone, PartialEq)]
pub struct Statement<F: Clone + Send + Sync> {
    pub constraints: Vec<SparseMatrix<F>>,
    pub public_input: Vec<F>,
}

/// A representation of a linearised CCS statement
#[derive(Debug, Clone, PartialEq)]
pub struct LStatement<F> {
    constraints: Vec<SparseMultilinearExtension<F>>,
    r: Vec<F>,
}

/// A representation of a CCS witness.
#[derive(Debug, Clone, PartialEq)]
pub struct Witness<F> {
    /// `w_ccs` is the original CCS witness.
    pub w_ccs: Vec<F>,
}

/// A representation of a linearised CCS witness.
#[derive(Debug, Clone, PartialEq)]
pub struct LWitness<F> {
    /// `w_ccs` is the original CCS witness.
    pub lw_ccs: DenseMultilinearExtension<F>,
}

impl<F: Ring> Witness<F> {
    /// Create a [`Witness`] from a ccs witness.
    pub fn new(w_ccs: Vec<F>) -> Self {
        Self { w_ccs }
    }

    /// Generates a random witness by firstly generating a random
    /// vector of arbitrary norm and then computing the rest of the data
    /// needed for a witness.
    ///
    /// # Arguments
    /// * `rng` is a mutable reference to the random number generator.
    /// * `w_ccs_len` is the length of the non-decomposed witness (a.k.a. the CCS witness).
    pub fn random<Rng: RngCore + ?Sized>(rng: &mut Rng, w_ccs_len: usize) -> Self {
        Self::new((0..w_ccs_len).map(|_| F::random(rng)).collect())
    }
}

/// A trait for defining the behaviour of a satisfying instance of a constraint system
///
/// # Types
///  - `R: Ring` - the ring in which the constraint system is operating.
///
pub trait Instance<F> {
    /// Given a witness vector, produce a concatonation of the statement and the witness
    fn get_z_vector(&self, w: &[F]) -> Vec<F>;
}

impl<F: Ring> Instance<F> for Statement<F> {
    fn get_z_vector(&self, w: &[F]) -> Vec<F> {
        let mut z: Vec<F> = Vec::with_capacity(self.public_input.len() + w.len() + 1);

        z.extend_from_slice(&self.public_input);
        z.push(F::one());
        z.extend_from_slice(w);

        z
    }
}

impl<const N: usize> Statement<Int<N>> {
    pub(crate) fn map_to_field<F: Field<LIMBS>, const LIMBS: usize>(&self) -> Statement<F> {
        Statement {
            constraints: self.constraints.iter().map(|m| m.map_to_field()).collect(),
            public_input: F::map_iterable(&self.public_input),
        }
    }
}

impl<F: Ring> Statement<F> {
    #[allow(non_snake_case)]
    pub fn compute_eval_table_sparse(
        &self,
        num_rows: usize,
        num_cols: usize,
        ccs: &CcsF<F>,
        evals: &[F],
    ) -> Vec<Vec<F>> {
        assert_eq!(num_rows, ccs.n);
        assert!(num_cols > (ccs.m - ccs.l) - 1);

        self.constraints
            .iter()
            .map(|M| compute_eval_table_sparse(M, evals, num_rows, num_cols))
            .collect()
    }
}

#[allow(non_snake_case)]
pub fn into_dense_matrix<F: Ring + From<u64>>(M: Vec<Vec<usize>>) -> Vec<Vec<F>> {
    M.iter()
        .map(|m| m.iter().map(|v| F::from(*v as u64)).collect())
        .collect()
}

pub fn into_vec<F: Ring + From<u64>>(z: Vec<u64>) -> Vec<F> {
    z.iter().map(|v| F::from(*v)).collect()
}

#[allow(non_snake_case)]
#[derive(Debug, Clone, PartialEq)]
pub struct CcsZ<F> {
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
    pub _phantom: PhantomData<F>,
}

impl<F: Ring> CcsZ<F> {
    pub fn pad(&mut self, statement: &mut Statement<F>, size: usize) {
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
                .for_each(|mat: &mut SparseMatrix<F>| {
                    mat.pad_cols(size);
                    mat.pad_rows(size);
                });
        }
    }
}

impl<const N: usize> CcsZ<Int<N>> {
    pub fn map_to_field<F: Field<LIMBS>, const LIMBS: usize>(&self) -> CcsF<F> {
        CcsF {
            m: self.m,
            n: self.n,
            l: self.l,
            t: self.t,
            q: self.q,
            d: self.d,
            s: self.s,
            s_prime: self.s_prime,
            S: self.S.clone(),
            c: self.c.iter().map(|v| F::from(*v)).collect(),
        }
    }
}

impl<F: Ring + From<i64>> Arith<F> for CcsZ<F> {
    /// check that a CCS structure is satisfied by a z vector. Only for testing.
    #[allow(non_snake_case)]
    fn check_relation(&self, M: &[SparseMatrix<F>], z: &[F]) -> Result<(), CSError> {
        let mut result = vec![F::ZERO; self.m];
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
            let vec_M_j: Vec<&SparseMatrix<F>> = self.S[i].iter().map(|j| &M[*j]).collect();

            // complete the hadamard chain
            let mut hadamard_result = vec![F::one(); self.m];
            for M_j in vec_M_j.into_iter() {
                let mut res = mat_vec_mul(M_j, z)?;
                res.resize(self.m, F::ZERO);
                hadamard_result = hadamard(&hadamard_result, &res)?;
            }

            // multiply by the coefficient of this step
            let c_M_j_z = vec_scalar_mul(&hadamard_result, &(F::from(self.c[i])));

            // add it to the final vector
            result = vec_add(&result, &c_M_j_z)?;
        }

        // make sure the final vector is all zeroes
        result
            .iter()
            .all(|item| *item == F::ZERO)
            .then_some(())
            .ok_or(CSError::NotSatisfied)
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

// Adds two ring vectors
pub(crate) fn vec_add<R: Clone + Add<R, Output = R>>(a: &[R], b: &[R]) -> Result<Vec<R>, CSError> {
    if a.len() != b.len() {
        return Err(CSError::LengthsNotEqual(
            String::from("a"),
            String::from("b"),
            a.len(),
            b.len(),
        ));
    }
    Ok(a.iter()
        .zip(b.iter())
        .map(|(x, y)| x.clone() + y.clone())
        .collect())
}

pub(crate) fn vec_scalar_mul<R: Clone + Mul<R, Output = R>>(vec: &[R], c: &R) -> Vec<R> {
    vec.iter().map(|a| a.clone() * c.clone()).collect()
}

pub(crate) fn hadamard<R: Clone + Mul<R, Output = R>>(a: &[R], b: &[R]) -> Result<Vec<R>, CSError> {
    if a.len() != b.len() {
        return Err(CSError::LengthsNotEqual(
            "a".into(),
            "b".into(),
            a.len(),
            b.len(),
        ));
    }
    Ok(a.iter()
        .zip(b)
        .map(|(a, b)| a.clone() * b.clone())
        .collect())
}

#[allow(non_snake_case)]
pub(crate) fn mat_vec_mul<R>(M: &SparseMatrix<R>, z: &[R]) -> Result<Vec<R>, CSError>
where
    R: Clone + Send + Sync + Mul<R, Output = R> + Add<Output = R> + Default,
    for<'a> R: Mul<&'a R, Output = R>,
{
    if M.n_cols != z.len() {
        return Err(CSError::LengthsNotEqual(
            "M".into(),
            "z".into(),
            M.n_cols,
            z.len(),
        ));
    }

    let mut result = Vec::with_capacity(M.coeffs.len());

    for row in &M.coeffs {
        let mut acc = R::default(); // Assuming Default gives the additive identity (e.g., 0)
        for (value, col_i) in row {
            acc = acc + (z[*col_i].clone() * value);
        }
        result.push(acc);
    }

    Ok(result)
}

#[allow(non_snake_case)]
pub fn get_dummy_ccs_Z_from_z_length<const I: usize>(
    n: usize,
    rng: &mut impl RngCore,
) -> (
    Vec<Int<I>>,
    CcsZ<Int<I>>,
    Statement<Int<I>>,
    Witness<Int<I>>,
) {
    let mut z: Vec<_> = (0..n).map(|_| Int::random(rng)).collect();
    let pub_io_len = 1;
    z[pub_io_len] = Int::one();
    let (ccs, statement, wit) = get_dummy_ccs_Z_from_z(&z, pub_io_len);

    (z, ccs, statement, wit)
}

#[allow(non_snake_case)]
fn get_dummy_ccs_Z_from_z<const I: usize>(
    z: &[Int<I>],
    pub_io_len: usize,
) -> (CcsZ<Int<I>>, Statement<Int<I>>, Witness<Int<I>>) {
    let ccs = CcsZ {
        m: z.len(),
        n: z.len(),
        l: pub_io_len,
        t: 3,
        q: 2,
        d: 2,
        s: log2(z.len()) as usize,
        s_prime: log2(z.len()) as usize,
        S: vec![vec![0, 1], vec![2]],
        c: vec![1, -1],
        _phantom: PhantomData,
    };

    let a = create_dummy_identity_sparse_matrix_Z(z.len(), z.len());
    let b = a.clone();
    let c = create_dummy_squaring_sparse_matrix_Z(z.len(), z.len(), z);

    let statement = Statement {
        constraints: vec![a, b, c],
        public_input: z[..pub_io_len].to_vec(),
    };

    let wit = Witness {
        w_ccs: z[pub_io_len + 1..].to_vec(),
    };

    (ccs, statement, wit)
}

#[allow(non_snake_case)]
pub(crate) fn create_dummy_identity_sparse_matrix_Z<const I: usize>(
    rows: usize,
    columns: usize,
) -> SparseMatrix<Int<I>> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((Int::one(), i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
#[allow(non_snake_case)]
pub(crate) fn create_dummy_squaring_sparse_matrix_Z<const I: usize>(
    rows: usize,
    columns: usize,
    witness: &[Int<I>],
) -> SparseMatrix<Int<I>> {
    assert_eq!(
        rows,
        witness.len(),
        "Length of witness vector must be equal to ccs width"
    );
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((witness[i], i));
    }
    matrix
}

#[cfg(test)]
mod tests {
    use std::marker::PhantomData;

    use crypto_bigint::{U128, const_monty_params};

    use super::*;
    use crate::{field::RandomField, sparse_matrix::SparseMatrix};

    const_monty_params!(ModP, U128, "0076F668F4274572E39A3EA8285319B5");
    type Fp = RandomField<ModP, { U128::LIMBS }>;

    #[test]
    fn test_statement_get_z_vector() {
        let constraints: Vec<SparseMatrix<Fp>> = vec![];
        let public_input = vec![Fp::from(2u64), Fp::from(3u64)];
        let stmt = Statement {
            constraints,
            public_input,
        };
        let w = vec![Fp::from(5u64), Fp::from(7u64)];
        let z = <Statement<Fp> as Instance<Fp>>::get_z_vector(&stmt, &w);
        assert_eq!(
            z,
            vec![
                Fp::from(2u64),
                Fp::from(3u64),
                Fp::from(1u64),
                Fp::from(5u64),
                Fp::from(7u64)
            ]
        );
    }

    #[allow(non_snake_case)]
    #[test]
    fn test_into_dense_matrix_and_vec() {
        let M_usize: Vec<Vec<usize>> = vec![vec![0, 2, 0], vec![5, 0, 7]];
        let M_fp: Vec<Vec<Fp>> = into_dense_matrix::<Fp>(M_usize.clone());
        let expected: Vec<Vec<Fp>> = vec![
            vec![Fp::from(0u64), Fp::from(2u64), Fp::from(0u64)],
            vec![Fp::from(5u64), Fp::from(0u64), Fp::from(7u64)],
        ];
        assert_eq!(M_fp, expected);

        let z_u64 = vec![0u64, 1u64, 10u64, 42u64];
        let z_fp = into_vec::<Fp>(z_u64.clone());
        let expected_z: Vec<Fp> = z_u64.iter().copied().map(Fp::from).collect();
        assert_eq!(z_fp, expected_z);
    }

    #[test]
    fn test_witness_new_and_random_len() {
        // new
        let w = vec![Fp::from(1u64), Fp::from(2u64), Fp::from(3u64)];
        let wit = Witness::new(w.clone());
        assert_eq!(wit.w_ccs, w);

        // random (length check) with a dummy RNG to satisfy trait bounds
        let mut rng = DummyRng::new(12345);
        let len = 8usize;
        let wr = Witness::<Fp>::random(&mut rng, len);
        assert_eq!(wr.w_ccs.len(), len);
    }

    struct DummyRng(u64);
    impl DummyRng {
        fn new(seed: u64) -> Self {
            Self(seed)
        }
    }
    impl RngCore for DummyRng {
        fn next_u32(&mut self) -> u32 {
            // simple xorshift32 from state
            let mut x = (self.0 as u32).wrapping_add(0x9E3779B9);
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            self.0 = (self.0 << 32) | x as u64;
            x
        }
        fn next_u64(&mut self) -> u64 {
            let hi = self.next_u32() as u64;
            let lo = self.next_u32() as u64;
            (hi << 32) | lo
        }
        fn fill_bytes(&mut self, dest: &mut [u8]) {
            let mut i = 0;
            while i < dest.len() {
                let v = self.next_u64().to_le_bytes();
                let take = core::cmp::min(8, dest.len() - i);
                dest[i..i + take].copy_from_slice(&v[..take]);
                i += take;
            }
        }
    }

    #[test]
    fn test_ccsz_pad_expands_and_updates_matrices() {
        // Start with m=3 so padding to size=5 should move to 8
        let mut ccsz: CcsZ<Fp> = CcsZ {
            m: 3,
            n: 3,
            l: 0,
            t: 1,
            q: 0,
            d: 0,
            s: 1,       // arbitrary; will be overwritten
            s_prime: 1, // arbitrary; will be overwritten
            S: vec![],
            c: vec![],
            _phantom: PhantomData,
        };

        // One constraint 3x3 with some entries
        let dense: Vec<Vec<Fp>> = vec![
            vec![Fp::from(1u64), Fp::from(0u64), Fp::from(2u64)],
            vec![Fp::from(0u64), Fp::from(0u64), Fp::from(0u64)],
            vec![Fp::from(3u64), Fp::from(0u64), Fp::from(4u64)],
        ];
        let mut statement = Statement::<Fp> {
            constraints: vec![SparseMatrix::from(dense)],
            public_input: vec![],
        };

        ccsz.pad(&mut statement, 5);

        assert_eq!(ccsz.m, 8);
        assert_eq!(ccsz.n, 8);
        assert_eq!(ccsz.s, 3);
        assert_eq!(ccsz.s_prime, 3);
        assert_eq!(statement.constraints.len(), 1);
        let sm = &statement.constraints[0];
        assert_eq!(sm.nrows(), 8);
        assert_eq!(sm.ncols(), 8);
    }

    #[test]
    fn test_ccsz_pad_noop_when_size_not_greater() {
        let mut ccsz: CcsZ<Fp> = CcsZ {
            m: 16,
            n: 16,
            l: 0,
            t: 0,
            q: 0,
            d: 0,
            s: 4,
            s_prime: 4,
            S: vec![],
            c: vec![],
            _phantom: PhantomData,
        };
        let mut statement = Statement::<Fp> {
            constraints: vec![SparseMatrix {
                n_rows: 16,
                n_cols: 16,
                coeffs: vec![],
            }],
            public_input: vec![],
        };
        ccsz.pad(&mut statement, 8);
        assert_eq!(ccsz.m, 16);
        assert_eq!(ccsz.n, 16);
        assert_eq!(ccsz.s, 4);
        assert_eq!(ccsz.s_prime, 4);
        let sm = &statement.constraints[0];
        assert_eq!(sm.nrows(), 16);
        assert_eq!(sm.ncols(), 16);
    }

    #[test]
    fn test_cserror_formatting_includes_values() {
        let err = CSError::MatricesRowsLengthNotPowerOf2(7);
        let msg = format!("{err}");
        assert!(msg.contains("7"));
    }
}
