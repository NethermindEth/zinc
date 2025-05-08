mod errors;
mod multilinear_polynomial;

//TODO add exports
pub use errors::ArithErrors;
pub(crate) use multilinear_polynomial::random_mle_list;

extern crate alloc;

// ark-std v0.5 should re-export alloc/std::sync.
// While already released on crates.io, related versioning changes were reverted on the GitHub repo.
// Let's wait for this issue to stabilize.
#[cfg(target_has_atomic = "ptr")]
pub use alloc::sync::Arc as RefCounter;
#[cfg(not(target_has_atomic = "ptr"))]
pub use ark_std::rc::Rc as RefCounter;
