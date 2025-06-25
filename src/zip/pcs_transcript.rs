#![allow(non_snake_case)]
use ark_std::{
    io::{Cursor, Read, Write},
    vec,
    vec::Vec,
};
use crypto_bigint::Int;
use sha3::{digest::Output, Keccak256};

use super::{pcs::utils::MerkleProof, Error};
use crate::{
    biginteger::BigInt,
    field::RandomField,
    field_config::ConfigRef,
    poly::alloc::string::ToString,
    traits::{Field, FromBytes},
    transcript::KeccakTranscript,
};

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

    pub fn common_field_element(&mut self, fe: &RandomField<N>) {
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
    pub fn write_field_elements(&mut self, elems: &[RandomField<N>]) -> Result<(), Error> {
        for elem in elems {
            self.write_field_element(elem)?;
        }

        Ok(())
    }

    pub fn read_field_elements<'cfg>(
        &mut self,
        n: usize,
        config: ConfigRef<'cfg, N>,
    ) -> Result<Vec<RandomField<'cfg, N>>, Error> {
        (0..n)
            .map(|_| self.read_field_element(config))
            .collect::<Result<Vec<_>, _>>()
    }

    pub fn read_field_element<'cfg>(
        &mut self,
        config: ConfigRef<'cfg, N>,
    ) -> Result<RandomField<'cfg, N>, Error> {
        let mut bytes: Vec<u8> = vec![0; N * 8];

        self.stream
            .read_exact(&mut bytes)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;

        let fe = RandomField::new_unchecked(config, BigInt::from_bytes_be(&bytes).unwrap());

        self.common_field_element(&fe);
        Ok(fe)
    }

    pub fn write_field_element(&mut self, fe: &RandomField<N>) -> Result<(), Error> {
        self.common_field_element(fe);
        let repr = fe.value().to_bytes_be();
        self.stream
            .write_all(repr.as_ref())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))
    }

    pub fn write_integer<const M: usize>(&mut self, int: &Int<M>) -> Result<(), Error> {
        for &word in int.as_words().iter() {
            let bytes = word.to_le_bytes();
            self.stream
                .write_all(&bytes)
                .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        }
        Ok(())
    }
    pub fn write_integers<const M: usize>(&mut self, ints: &[Int<M>]) -> Result<(), Error> {
        for int in ints {
            self.write_integer(int)?;
        }
        Ok(())
    }

    pub fn read_integer<const M: usize>(&mut self) -> Result<Int<M>, Error> {
        let mut words = [0u64; M];

        for word in &mut words {
            let mut buf = [0u8; 8];
            self.stream
                .read_exact(&mut buf)
                .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;

            *word = u64::from_le_bytes(buf);
        }
        Ok(Int::<M>::from_words(words))
    }

    pub fn read_integers<const M: usize>(&mut self, n: usize) -> Result<Vec<Int<M>>, Error> {
        (0..n)
            .map(|_| self.read_integer())
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

    pub fn squeeze_challenge_idx(&mut self, config: ConfigRef<N>, cap: usize) -> usize {
        let challenge: RandomField<N> = self.fs_transcript.get_challenge(config);
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

#[allow(unused_macros)]
macro_rules! test_read_write {
    ($write_fn:ident, $read_fn:ident, $original_value:expr, $assert_msg:expr) => {{
        use ark_std::format;
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

#[allow(unused_macros)]
macro_rules! test_read_write_vec {
    ($write_fn:ident, $read_fn:ident, $original_values:expr, $assert_msg:expr) => {{
        use ark_std::format;
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

#[test]
fn test_pcs_transcript_read_write() {
    const N: usize = 4;

    // Test commitment
    let original_commitment = Output::<Keccak256>::default();
    test_read_write!(
        write_commitment,
        read_commitment,
        original_commitment,
        "commitment"
    );
    //TODO put the tests back in for Int<N> types
    // Test vector of commitments
    let original_commitments = vec![Output::<Keccak256>::default(); 1024];
    test_read_write_vec!(
        write_commitments,
        read_commitments,
        original_commitments,
        "commitments vector"
    );
}
