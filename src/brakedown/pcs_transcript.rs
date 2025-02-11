use std::io::{Cursor, Read, Write};

use sha3::digest::Output;
use sha3::Keccak256;

use crate::biginteger::BigInt;
use crate::field::RandomField as F;
use crate::field_config::FieldConfig;
use crate::traits::FromBytes;
use crate::transcript::KeccakTranscript;

use super::Error;

#[derive(Default)]
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
        let mut bytes: [u8; N] = [0; N];

        self.stream
            .read_exact(&mut bytes)
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;

        let fe =
            F::from_bigint(config, BigInt::from_bytes_be(&bytes).unwrap()).ok_or_else(|| {
                Error::Transcript(
                    std::io::ErrorKind::Other,
                    "Invalid field element encoding in proof".to_string(),
                )
            })?;
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
}
