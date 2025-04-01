use std::{marker::PhantomData, slice};

use ark_std::rand::RngCore;
use sha3::{digest::Output, Keccak256};

use crate::{
    brakedown::code::{Brakedown, BrakedownSpec, LinearCodes},
    field::RandomField as F,
    field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension,
};

#[derive(Debug)]
pub struct MultilinearBrakedown<const N: usize, S: BrakedownSpec>(PhantomData<S>);

impl<const N: usize, S: BrakedownSpec> Clone for MultilinearBrakedown<N, S> {
    fn clone(&self) -> Self {
        Self(PhantomData)
    }
}

#[derive(Clone, Debug)]
pub struct MultilinearBrakedownParams<const N: usize> {
    num_vars: usize,
    num_rows: usize,
    brakedown: Brakedown<N>,
}

impl<const N: usize> MultilinearBrakedownParams<N> {
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    pub fn brakedown(&self) -> &Brakedown<N> {
        &self.brakedown
    }
}

/// Representantation of a brakedown commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearBrakedownCommitment<const N: usize> {
    /// The encoded rows of the polynomial matrix representation
    rows: Vec<F<N>>,
    /// Hashes of the merkle tree with the encoded columns as leaves
    intermediate_hashes: Vec<Output<Keccak256>>,
    /// Root of the merkle tree with the encoded columns as leaves
    root: Output<Keccak256>,
}

impl<const N: usize> MultilinearBrakedownCommitment<N> {
    pub fn new(
        rows: Vec<F<N>>,
        intermediate_hashes: Vec<Output<Keccak256>>,
        root: Output<Keccak256>,
    ) -> MultilinearBrakedownCommitment<N> {
        MultilinearBrakedownCommitment {
            rows,
            intermediate_hashes,
            root,
        }
    }
    pub fn from_root(root: Output<Keccak256>) -> Self {
        Self {
            root,
            ..Default::default()
        }
    }

    pub fn rows(&self) -> &[F<N>] {
        &self.rows
    }

    pub fn intermediate_hashes(&self) -> &[Output<Keccak256>] {
        &self.intermediate_hashes
    }

    pub fn root(&self) -> &Output<Keccak256> {
        &self.root
    }
}

impl<const N: usize> AsRef<[Output<Keccak256>]> for MultilinearBrakedownCommitment<N> {
    fn as_ref(&self) -> &[Output<Keccak256>] {
        slice::from_ref(&self.root)
    }
}

impl<const N: usize, S> MultilinearBrakedown<N, S>
where
    S: BrakedownSpec,
{
    pub type Param = MultilinearBrakedownParams<N>;
    pub type ProverParam = MultilinearBrakedownParams<N>;
    pub type VerifierParam = MultilinearBrakedownParams<N>;
    pub type Polynomial = DenseMultilinearExtension<N>;
    pub type Commitment = MultilinearBrakedownCommitment<N>;
    pub type CommitmentChunk = Output<Keccak256>;

    pub fn setup(
        poly_size: usize,
        rng: impl RngCore,
        config: *const FieldConfig<N>,
    ) -> Self::Param {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        let brakedown =
            Brakedown::new_multilinear::<S>(num_vars, 20.min((1 << num_vars) - 1), rng, config);
        MultilinearBrakedownParams {
            num_vars,
            num_rows: (1 << num_vars) / brakedown.row_len(),
            brakedown,
        }
    }
}
