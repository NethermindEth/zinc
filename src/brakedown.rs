pub mod code;
pub mod pcs;
pub mod pcs_transcript;
pub mod utils;

#[cfg(test)]
mod tests;
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    InvalidSumcheck(String),
    InvalidPcsParam(String),
    InvalidPcsOpen(String),
    InvalidSnark(String),
    Serialization(String),
    Transcript(std::io::ErrorKind, String),
}
