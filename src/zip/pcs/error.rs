use thiserror::Error;

#[derive(Error, Debug)]
pub enum MerkleError {
    #[error("Invalid Merkle proof: {0}")]
    InvalidMerkleProof(String),

    #[error("Failed to read merkle proof")]
    FailedMerkleProofReading,

    #[error("Failed to write merkle proof")]
    FailedMerkleProofWriting,
}
