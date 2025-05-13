//! Zinc crate implements the Zinc Protocol[https://eprint.iacr.org/2025/316]
//! A hash-based succinct argument system for integer arithmetic.
//!
//! This crate implements a practically efficient, hash-based succinct proof protocol
//! designed to bypass the significant arithmetization overhead of many existing schemes.
//! It enables proving statements directly over integers (ℤ) or rationals (ℚ), supporting
//! modular operations over arbitrary moduli (including non-prime and multiple moduli).
//!
//! The protocol combines:
//! 1. A reduction framework that maps rational statements modulo a random prime;
//! 2. A novel integer-proximity IOP-based polynomial commitment scheme (Brakedown-style)
//!    ensuring polynomial coefficients remain close to integral.
//!
//! This results in a pure code-and-hash-based system without hidden-order groups,
//! achieving efficient succinct proofs for integer arithmetic with minimal overhead.

#![allow(incomplete_features)]
#![feature(inherent_associated_types, generic_const_exprs)]
#![warn(missing_docs)]
pub mod biginteger;
pub mod ccs;
mod const_helpers;
pub mod field;
pub mod field_config;
pub mod poly_f;
pub mod poly_z;
pub mod prime_gen;
pub(crate) mod sparse_matrix;
pub mod sumcheck;
pub(crate) mod traits;
pub mod transcript;
pub mod zinc;
pub mod zip;
