use crate::{field::RandomField as F, transcript::KeccakTranscript};
pub struct MultisetHashes<const N: usize> {
    /// Multiset hash of "read" tuples
    pub read_hashes: Vec<F<N>>,
    /// Multiset hash of "write" tuples
    pub write_hashes: Vec<F<N>>,
    /// Multiset hash of "init" tuples
    pub init_hashes: Vec<F<N>>,
    /// Multiset hash of "final" tuples
    pub final_hashes: Vec<F<N>>,
}
impl<const N: usize> MultisetHashes<N> {
    pub fn append_to_transcript(&self, transcript: &mut KeccakTranscript) {
        transcript.absorb_slice(&self.read_hashes);
        transcript.absorb_slice(&self.write_hashes);
        transcript.absorb_slice(&self.init_hashes);
        transcript.absorb_slice(&self.final_hashes);
    }
}
