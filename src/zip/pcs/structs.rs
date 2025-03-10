use std::{marker::PhantomData, slice};

use ark_std::rand::RngCore;
use i256::I256;
use sha3::{digest::Output, Keccak256};

use crate::{
    poly_z::mle::DenseMultilinearExtension,
    zip::code::{LinearCodes, Zip, ZipSpec},
};

#[derive(Debug)]
pub struct MultilinearZip<const N: usize, S: ZipSpec>(PhantomData<S>);

impl<const N: usize, S: ZipSpec> Clone for MultilinearZip<N, S> {
    fn clone(&self) -> Self {
        Self(PhantomData)
    }
}

#[derive(Clone, Debug)]
pub struct MultilinearZipParams<const N: usize> {
    num_vars: usize,
    num_rows: usize,
    zip: Zip<N>,
}

impl<const N: usize> MultilinearZipParams<N> {
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    pub fn zip(&self) -> &Zip<N> {
        &self.zip
    }
}

/// Representantation of a zip commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearZipCommitment<const N: usize> {
    /// The encoded rows of the polynomial matrix representation
    rows: Vec<I256>,
    /// Hashes of the merkle tree with the encoded columns as leaves
    intermediate_hashes: Vec<Output<Keccak256>>,
    /// Root of the merkle tree with the encoded columns as leaves
    root: Output<Keccak256>,
}

impl<const N: usize> MultilinearZipCommitment<N> {
    pub fn new(
        rows: Vec<I256>,
        intermediate_hashes: Vec<Output<Keccak256>>,
        root: Output<Keccak256>,
    ) -> MultilinearZipCommitment<N> {
        MultilinearZipCommitment {
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

    pub fn rows(&self) -> &[I256] {
        &self.rows
    }

    pub fn intermediate_hashes(&self) -> &[Output<Keccak256>] {
        &self.intermediate_hashes
    }

    pub fn root(&self) -> &Output<Keccak256> {
        &self.root
    }
}

impl<const N: usize> AsRef<[Output<Keccak256>]> for MultilinearZipCommitment<N> {
    fn as_ref(&self) -> &[Output<Keccak256>] {
        slice::from_ref(&self.root)
    }
}

impl<const N: usize, S> MultilinearZip<N, S>
where
    S: ZipSpec,
{
    pub type Param = MultilinearZipParams<N>;
    pub type ProverParam = MultilinearZipParams<N>;
    pub type VerifierParam = MultilinearZipParams<N>;
    pub type Polynomial = DenseMultilinearExtension;
    pub type Commitment = MultilinearZipCommitment<N>;
    pub type CommitmentChunk = Output<Keccak256>;

    pub fn setup(poly_size: usize, rng: impl RngCore) -> Self::Param {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        let zip = Zip::new_multilinear::<S>(num_vars, 20.min((1 << num_vars) - 1), rng);

        MultilinearZipParams {
            num_vars,
            num_rows: ((1 << num_vars) / zip.row_len()).next_power_of_two(),
            zip,
        }
    }
}
