//! The core SNARK over the integers
#![allow(non_snake_case)]

mod errors;

pub mod prover;
pub mod structs;
#[cfg(test)]
mod tests;
// TODO: private utils
pub mod utils;
pub mod verifier;

pub mod prelude {
    pub use crate::{
        ccs::ccs_z::*,
        transcript::KeccakTranscript,
        zinc::utils::draw_random_field,
        zinc::{
            prover::Prover,
            structs::{ZincProver, ZincVerifier},
            verifier::Verifier,
        },
        zip::code::ZipSpec1,
    };
}
