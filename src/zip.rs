//! Defines behaviour for Zip, a PCS over the integers
pub mod code;
pub mod pcs;
pub mod pcs_transcript;
mod utils;

#[cfg(test)]
mod tests;
use thiserror::Error;
#[derive(Clone, Debug, PartialEq, Error)]
/// Errors that can occur in the Zip PCS scheme
pub enum ZipError {
    /// Parameters do not match
    #[error("Invalid PCS param: {0}")]
    InvalidPcsParam(String),
    /// The proof of an evaluation was not valid
    #[error("Invalid commitment opening: {0}")]
    InvalidPcsOpen(String),
    /// Not able to read or write from the shared transcript
    #[error("Transcript failure: {1}")]
    Transcript(std::io::ErrorKind, String),
}
