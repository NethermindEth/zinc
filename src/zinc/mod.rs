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
        zinc::{
            prover::Prover,
            structs::{ZincProver, ZincVerifier},
            utils::draw_random_field,
            verifier::Verifier,
        },
        zip::code::ZipLinearCodeSpec1,
    };
}
