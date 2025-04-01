use crate::{
    brakedown::pcs::structs::MultilinearBrakedownCommitment, field::RandomField as F,
    poly_f::polynomials::dense_interleaved_polynomial::DenseInterleavedPolynomial,
    sumcheck::SumcheckProof,
};

pub struct QuarkGrandProductProof<const N: usize> {
    sumcheck_proof: SumcheckProof<N>,
    g_commitment: MultilinearBrakedownCommitment<N>,
    g_r_sumcheck: F<N>,
    g_r_prime: (F<N>, F<N>),
    v_r_prime: (F<N>, F<N>),
    pub num_vars: usize,
}
pub struct QuarkGrandProduct<const N: usize> {
    batch_size: usize,
    quark_poly: Option<Vec<F<N>>>,
    base_layers: Vec<DenseInterleavedPolynomial<N>>,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct QuarkGrandProductConfig {
    pub hybrid_layer_depth: QuarkHybridLayerDepth,
}
#[derive(Clone, Copy, Debug, Default)]
pub enum QuarkHybridLayerDepth {
    #[default]
    Default,
    Min,
    Max,
    Custom(usize),
}

impl QuarkHybridLayerDepth {
    // The depth in the product tree of the grand product at which the
    // hybrid implementation will switch to using quarks grand product proofs
    pub fn get_crossover_depth(&self) -> usize {
        match self {
            QuarkHybridLayerDepth::Min => 0, // Always use quarks
            QuarkHybridLayerDepth::Default => 4,
            QuarkHybridLayerDepth::Max => usize::MAX, // Never use quarks
            QuarkHybridLayerDepth::Custom(depth) => *depth,
        }
    }
}
