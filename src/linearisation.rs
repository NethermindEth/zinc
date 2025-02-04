#![allow(dead_code)]
use errors::LinearizationError;
use structs::LinearisationProof;

use crate::{
    ccs::ccs_f::{LStatement, LWitness, Statement, Witness, CCS_F},
    transcript::KeccakTranscript,
};

mod errors;
mod structs;
mod utils;
/// Prover for the Linearization subprotocol
pub trait LinearisationProver<const N: usize> {
    /// Generates a proof for the linearization subprotocol
    ///
    /// # Arguments
    ///
    /// * `cm_i` - A reference to a committed CCS statement to be linearized, i.e. a CCCS<C, NTT>.
    /// * `wit` - A reference to a CCS witness for the statement cm_i.
    /// * `transcript` - A mutable reference to a sponge for generating NI challenges.
    /// * `ccs` - A reference to a Customizable Constraint System circuit representation.
    ///
    /// # Returns
    ///
    /// On success, returns a tuple `(LCCCS<C, NTT>, LinearizationProof<NTT>)` where:
    ///   * `LCCCS<C, NTT>` is a linearized version of the CCS witness commitment.
    ///   * `LinearizationProof<NTT>` is a proof that the linearization subprotocol was executed correctly.
    ///
    /// # Errors
    ///
    /// Returns an error if asked to evaluate MLEs with incorrect number of variables
    ///
    fn prove(
        statement: &Statement<N>,
        wit: &Witness<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<(LStatement<N>, LWitness<N>, LinearisationProof<N>), LinearizationError<N>>;
}

/// Verifier for the Linearization subprotocol.
pub trait LinearisationVerifier<const N: usize> {
    /// Verifies a proof for the linearization subprotocol.
    ///
    /// # Arguments
    ///
    /// * `cm_i` - A reference to a `CCCS<C, NTT>`, which represents a CCS statement and a commitment to a witness.
    /// * `proof` - A reference to a `LinearizationProof<NTT>` containing the linearization proof.
    /// * `transcript` - A mutable reference to a sponge for generating NI challenges.
    /// * `ccs` - A reference to a Customizable Constraint System instance used in the protocol.
    ///
    /// # Returns
    ///
    /// * `Ok(LCCCS<C, NTT>)` - On success, returns a linearized version of the CCS witness commitment.
    /// * `Err(LinearizationError<NTT>)` - If verification fails, returns a `LinearizationError<NTT>`.
    ///
    fn verify(
        cm_i: &Statement<N>,
        proof: &LinearisationProof<N>,
        transcript: &mut KeccakTranscript,
        ccs: &CCS_F<N>,
    ) -> Result<LStatement<N>, LinearizationError<N>>;
}
