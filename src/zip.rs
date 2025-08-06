pub mod code;
pub mod code_raa;
pub mod pcs;
pub mod pcs_transcript;
pub mod utils;

#[cfg(test)]
mod tests;

use ark_std::string::String;
use thiserror::Error;
#[derive(Clone, Debug, PartialEq, Error)]
pub enum Error {
    #[error("Invalid PCS param: {0}")]
    InvalidPcsParam(String),
    #[error("Invalid commitment opening: {0}")]
    InvalidPcsOpen(String),
    #[error("Bad Snark: {0}")]
    InvalidSnark(String),
    #[error("Serialization Error: {0}")]
    Serialization(String),
    #[error("Transcript failure: {1}")]
    Transcript(ark_std::io::ErrorKind, String),
}
