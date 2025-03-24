use std::{collections::BTreeSet, marker::PhantomData, slice};

use i256::I512;
use sha3::{digest::Output, Keccak256};

use crate::{
    poly_z::mle::DenseMultilinearExtension,
    zip::code::{LinearCodes, Zip, ZipSpec},
};

#[derive(Debug)]
pub struct MultilinearZip<const N: usize, S: ZipSpec, T: ZipTranscript>(
    PhantomData<S>,
    PhantomData<T>,
);

impl<const N: usize, S: ZipSpec, T: ZipTranscript> Clone for MultilinearZip<N, S, T> {
    fn clone(&self) -> Self {
        Self(PhantomData, PhantomData)
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
    rows: Vec<I512>,
    /// Hashes of the merkle tree with the encoded columns as leaves
    intermediate_hashes: Vec<Output<Keccak256>>,
    /// Root of the merkle tree with the encoded columns as leaves
    root: Output<Keccak256>,
}

impl<const N: usize> MultilinearZipCommitment<N> {
    pub fn new(
        rows: Vec<I512>,
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

    pub fn rows(&self) -> &[I512] {
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

pub trait ZipTranscript {
    fn get_encoding_element(&mut self) -> i128;
    fn sample_unique_columns(
        &mut self,
        range: std::ops::Range<usize>,
        columns: &mut BTreeSet<usize>,
        count: usize,
    ) -> usize;
}
impl<const N: usize, S, T> MultilinearZip<N, S, T>
where
    S: ZipSpec,
    T: ZipTranscript,
{
    pub type Param = MultilinearZipParams<N>;
    pub type ProverParam = MultilinearZipParams<N>;
    pub type VerifierParam = MultilinearZipParams<N>;
    pub type Polynomial = DenseMultilinearExtension;
    pub type Commitment = MultilinearZipCommitment<N>;
    pub type CommitmentChunk = Output<Keccak256>;

    pub fn setup(poly_size: usize, transcript: &mut T) -> Self::Param {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        let zip = Zip::new_multilinear::<S, T>(num_vars, 20.min((1 << num_vars) - 1), transcript);

        MultilinearZipParams {
            num_vars,
            num_rows: ((1 << num_vars) / zip.row_len()).next_power_of_two(),
            zip,
        }
    }
}
