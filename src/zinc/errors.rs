#![allow(clippy::enum_variant_names)]
use thiserror::Error;

use crate::{
    ccs::error::CSError, poly_f::polynomials::ArithErrors, sumcheck::SumCheckError, zip::ZipError,
};

#[derive(Debug, Error)]
pub enum ZincError<const N: usize> {
    #[error("lookup error: {0}")]
    LookupError(#[from] LookupError),
    #[error("spartan error: {0}")]
    SpartanError(#[from] SpartanError<N>),
}

#[derive(Debug, Error)]
pub enum SpartanError<const N: usize> {
    #[error("sum check failed at linearization step: {0}")]
    SumCheckError(#[from] SumCheckError<N>),
    #[error("parameters error: {0}")]
    ParametersError(String),
    #[error("constraint system related error: {0}")]
    ConstraintSystemError(#[from] CSError),
    #[error("Arithmetic error: {0}")]
    ArithmeticError(#[from] ArithErrors),
    #[error("mle evaluation failed: {0}")]
    EvaluationError(#[from] MleEvaluationError),
    #[error("Verification Failed {0}")]
    VerificationError(String),
    #[error("")]
    ZipError(#[from] ZipError),
}

#[derive(Debug, Error)]
pub enum LookupError {
    #[error("lookup failed")]
    LookupFailed,
}

#[derive(Debug, Error)]
pub enum MleEvaluationError {
    #[error("lengths of evaluation point and evaluations are not consistent: 1 << {0} != {1}")]
    IncorrectLength(usize, usize),
}
