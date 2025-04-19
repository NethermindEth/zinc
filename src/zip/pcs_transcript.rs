#![allow(non_snake_case)]
use std::io::{Cursor, Read, Write};

use i256::{I256, I512};
use sha3::digest::Output;
use sha3::Keccak256;

use crate::biginteger::BigInt;
use crate::field::RandomField as F;
use crate::field_config::FieldConfig;
use crate::traits::FromBytes;
use crate::transcript::KeccakTranscript;

use super::pcs::utils::MerkleProof;
use super::Error;

#[derive(Default, Clone)]
pub struct PcsTranscript<const N: usize> {
    pub fs_transcript: KeccakTranscript,
    pub stream: Cursor<Vec<u8>>,
}

impl<const N: usize> PcsTranscript<N> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn into_proof(self) -> Vec<u8> {
        self.stream.into_inner()
    }

    pub fn from_proof(proof: &[u8]) -> Self {
        Self {
            fs_transcript: KeccakTranscript::default(),
            stream: Cursor::new(proof.to_vec()),
        }
    }

    pub fn common_field_element(&mut self, fe: &F<N>) {
        self.fs_transcript.absorb_random_field(fe);
    }

    pub fn read_commitment(&mut self) -> Result<Output<Keccak256>, Error> {
        let mut buf = Output::<Keccak256>::default();
        self.stream
            .read_exact(&mut buf)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        Ok(*Output::<Keccak256>::from_slice(&buf))
    }

    pub fn write_commitment(&mut self, comm: &Output<Keccak256>) -> Result<(), Error> {
        self.stream
            .write_all(comm)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        Ok(())
    }

    // TODO if we change this to an iterator we may be able to save some memory
    pub fn write_field_elements(&mut self, elems: &[F<N>]) -> Result<(), Error> {
        for elem in elems {
            self.write_field_element(elem)?;
        }

        Ok(())
    }

    pub fn read_field_elements(
        &mut self,
        n: usize,
        config: *const FieldConfig<N>,
    ) -> Result<Vec<F<N>>, Error> {
        (0..n)
            .map(|_| self.read_field_element(config))
            .collect::<Result<Vec<_>, _>>()
    }

    pub fn read_field_element(&mut self, config: *const FieldConfig<N>) -> Result<F<N>, Error> {
        let mut bytes: Vec<u8> = vec![0; N * 8];

        self.stream
            .read_exact(&mut bytes)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;

        let fe = F::new_unchecked(config, BigInt::from_bytes_be(&bytes).unwrap());

        self.common_field_element(&fe);
        Ok(fe)
    }

    pub fn write_field_element(&mut self, fe: &F<N>) -> Result<(), Error> {
        self.common_field_element(fe);
        let repr = fe.value().to_bytes_be();
        self.stream
            .write_all(repr.as_ref())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))
    }

    pub fn write_integer(&mut self, int: &i64) -> Result<(), Error> {
        self.stream
            .write_all(int.to_be_bytes().as_ref())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))
    }

    pub fn write_integers(&mut self, ints: &[i64]) -> Result<(), Error> {
        for int in ints {
            self.write_integer(int)?;
        }
        Ok(())
    }
    pub fn write_I256(&mut self, int: &I256) -> Result<(), Error> {
        self.stream
            .write_all(int.to_be_bytes().as_ref())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))
    }

    pub fn write_I256_vec(&mut self, ints: &[I256]) -> Result<(), Error> {
        for int in ints {
            self.write_I256(int)?;
        }
        Ok(())
    }

    pub fn write_I512(&mut self, int: &I512) -> Result<(), Error> {
        self.stream
            .write_all(int.to_be_bytes().as_ref())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))
    }

    pub fn write_I512_vec(&mut self, ints: &[I512]) -> Result<(), Error> {
        for int in ints {
            self.write_I512(int)?;
        }
        Ok(())
    }
    pub fn read_integer(&mut self) -> Result<i64, Error> {
        let mut bytes = [0; 8];

        self.stream
            .read_exact(&mut bytes)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        Ok(i64::from_be_bytes(bytes))
    }

    pub fn read_integers(&mut self, n: usize) -> Result<Vec<i64>, Error> {
        (0..n)
            .map(|_| self.read_integer())
            .collect::<Result<Vec<_>, _>>()
    }

    pub fn read_I256(&mut self) -> Result<I256, Error> {
        let mut bytes = [0; 32];

        self.stream
            .read_exact(&mut bytes)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        Ok(I256::from_be_bytes(bytes))
    }

    pub fn read_I256_vec(&mut self, n: usize) -> Result<Vec<I256>, Error> {
        (0..n)
            .map(|_| self.read_I256())
            .collect::<Result<Vec<_>, _>>()
    }

    pub fn read_I512(&mut self) -> Result<I512, Error> {
        let mut bytes = [0; 64];

        self.stream
            .read_exact(&mut bytes)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        Ok(I512::from_be_bytes(bytes))
    }

    pub fn read_I512_vec(&mut self, n: usize) -> Result<Vec<I512>, Error> {
        (0..n)
            .map(|_| self.read_I512())
            .collect::<Result<Vec<_>, _>>()
    }

    pub fn read_commitments(&mut self, n: usize) -> Result<Vec<Output<Keccak256>>, Error> {
        (0..n).map(|_| self.read_commitment()).collect()
    }

    pub fn write_commitments<'a>(
        &mut self,
        comms: impl IntoIterator<Item = &'a Output<Keccak256>>,
    ) -> Result<(), Error> {
        for comm in comms.into_iter() {
            self.write_commitment(comm)?;
        }
        Ok(())
    }

    pub fn squeeze_challenge_idx(&mut self, config: *const FieldConfig<N>, cap: usize) -> usize {
        let challenge = self.fs_transcript.get_challenge(config);
        let bytes = challenge.value().to_bytes_le();
        let num = u32::from_le_bytes(bytes[..4].try_into().unwrap()) as usize;
        num % cap
    }

    pub fn read_merkle_proof(&mut self) -> Result<MerkleProof, Error> {
        // Read the length of the merkle_path first
        let mut length_bytes = [0u8; 8];
        self.stream
            .read_exact(&mut length_bytes)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        let path_length = u64::from_be_bytes(length_bytes);

        // Read each element of the merkle_path
        let mut merkle_path = Vec::with_capacity(path_length as usize);
        for _ in 0..path_length {
            merkle_path.push(self.read_commitment()?);
        }

        Ok(MerkleProof { merkle_path })
    }

    pub fn write_merkle_proof(&mut self, proof: &MerkleProof) -> Result<(), Error> {
        // Write the length of the merkle_path first
        let path_length = proof.merkle_path.len() as u64;
        self.stream
            .write_all(&path_length.to_be_bytes())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;

        // Write each element of the merkle_path
        for hash in &proof.merkle_path {
            self.write_commitment(hash)?;
        }

        Ok(())
    }
}

macro_rules! test_read_write {
    ($write_fn:ident, $read_fn:ident, $original_value:expr, $assert_msg:expr) => {{
        let mut transcript = PcsTranscript::<N>::new();
        transcript
            .$write_fn(&$original_value)
            .expect(&format!("Failed to write {}", $assert_msg));
        let proof = transcript.into_proof();
        let mut transcript = PcsTranscript::<N>::from_proof(&proof);
        let read_value = transcript
            .$read_fn()
            .expect(&format!("Failed to read {}", $assert_msg));
        assert_eq!(
            $original_value, read_value,
            "{} read does not match original",
            $assert_msg
        );
    }};
}

macro_rules! test_read_write_vec {
    ($write_fn:ident, $read_fn:ident, $original_values:expr, $assert_msg:expr) => {{
        let mut transcript = PcsTranscript::<N>::new();
        transcript
            .$write_fn(&$original_values)
            .expect(&format!("Failed to write {}", $assert_msg));
        let proof = transcript.into_proof();
        let mut transcript = PcsTranscript::<N>::from_proof(&proof);
        let read_values = transcript
            .$read_fn($original_values.len())
            .expect(&format!("Failed to read {}", $assert_msg));
        assert_eq!(
            $original_values, read_values,
            "{} read does not match original",
            $assert_msg
        );
    }};
}

const N: usize = 4;

#[test]
fn test_pcs_transcript_read_write() {
    const N: usize = 4;

    // Test integer
    let original_int: i64 = 42;
    test_read_write!(write_integer, read_integer, original_int, "integer");

    // Test I256
    let original_i256 = I256::from(340282366920938463463374607431768211455u128);
    test_read_write!(write_I256, read_I256, original_i256, "I256");

    // Test I512
    let bytes_array = [42u8; 64];
    let original_i512 = I512::from_be_bytes(bytes_array);
    test_read_write!(write_I512, read_I512, original_i512, "I512");

    // Test commitment
    let original_commitment = Output::<Keccak256>::default();
    test_read_write!(
        write_commitment,
        read_commitment,
        original_commitment,
        "commitment"
    );

    // Test vector of integers
    let original_ints = vec![1i64; 1024];
    test_read_write_vec!(write_integers, read_integers, original_ints, "integers");

    // Test vector of I256
    let original_i256_vec = vec![I256::from(1); 1024];
    test_read_write_vec!(
        write_I256_vec,
        read_I256_vec,
        original_i256_vec,
        "I256 vector"
    );

    // Test vector of I512
    let original_i512_vec = vec![I512::from(1); 1024];
    test_read_write_vec!(
        write_I512_vec,
        read_I512_vec,
        original_i512_vec,
        "I512 vector"
    );

    // Test vector of commitments
    let original_commitments = vec![Output::<Keccak256>::default(); 1024];
    test_read_write_vec!(
        write_commitments,
        read_commitments,
        original_commitments,
        "commitments vector"
    );
}
