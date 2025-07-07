#![allow(clippy::enum_variant_names)]
use ark_std::string::String;
use thiserror::Error;

use crate::poly::ArithErrors;
use crate::{ccs::error::CSError, sumcheck::SumCheckError, zip::Error as ZipError};

#[derive(Debug, Error)]
pub enum ZincError {
    #[error("spartan error: {0}")]
    SpartanError(#[from] SpartanError),
    #[error("field config error")]
    FieldConfigError,
}

#[derive(Debug, Error)]
pub enum SpartanError {
    #[error("sum check failed at linearization step: {0}")]
    SumCheckError(#[from] SumCheckError),
    #[error("parameters error: {0}")]
    ParametersError(String),
    #[error("constraint system related error: {0}")]
    ConstraintSystemError(#[from] CSError),
    #[error("Arithmetic error: {0}")]
    ArithmeticError(#[from] ArithErrors),
    #[error("mle evaluation failed: {0}")]
    EvaluationError(#[from] MleEvaluationError),
    #[error("Verification Failed {0}")]
    PCSVerificationError(String),
    #[error("")]
    ZipError(#[from] ZipError),
}

#[derive(Debug, Error)]
pub enum MleEvaluationError {
    #[error("lengths of evaluation point and evaluations are not consistent: 1 << {0} != {1}")]
    IncorrectLength(usize, usize),
}
