//! Utility types for the Zinc protocol
use std::marker::PhantomData;

use crate::{
    field::RandomField,
    sumcheck,
    zip::{code::ZipSpec, pcs::structs::MultilinearZipCommitment},
};

/// Non-interactive proof of the Spartan Subprotocol
#[derive(Debug, Clone)]
pub struct SpartanProof<const N: usize> {
    /// A proof of the sumcheck in step 5
    /// https://eprint.iacr.org/2019/550.pdf#page=20  
    pub first_sumcheck: sumcheck::SumcheckProof<N>,
    /// A proof of the sumcheck in step 11
    /// https://eprint.iacr.org/2019/550.pdf#page=20  
    pub second_sumcheck: sumcheck::SumcheckProof<N>,
    /// Evaluations of the matrix extensions at the first
    /// sumcheck challenge point.
    pub V_s: Vec<RandomField<N>>,
}

/// A proof of a PCS opening
pub struct ZipProof<const N: usize> {
    pub(crate) z_comm: MultilinearZipCommitment<N>,
    pub(crate) v: RandomField<N>,
    pub(crate) pcs_proof: Vec<u8>,
}

/// Proof for the lookup performed at the end
/// of the Zinc Protocol
pub struct LookupProof<const N: usize> {}

/// Proof for the lookup performed at the end
/// of the Zinc Protocol
pub struct ZincProof<const N: usize> {
    /// The proof of the spartan subprotocol
    pub spartan_proof: SpartanProof<N>,
    /// The PCS opening proof
    pub zip_proof: ZipProof<N>,
    /// The lookup proof that all witness elements
    /// are integers
    pub lookup_proof: LookupProof<N>,
}

/// The implementation of the `LinearizationProver` trait is defined in the main linearization file.
pub struct ZincProver<const N: usize, S>
where
    S: ZipSpec,
{
    /// holds the parameters of the underlying PCS
    pub data: PhantomData<S>,
}

/// The implementation of the `LinearizationVerifier` trait is defined in the main linearization file.
pub struct ZincVerifier<const N: usize, S>
where
    S: ZipSpec,
{
    /// holds the parameters of the underlying PCS
    pub data: PhantomData<S>,
}
