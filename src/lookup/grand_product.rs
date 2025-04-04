use crate::{brakedown::code::BrakedownSpec, field::RandomField, sumcheck::SumcheckProof};

use super::grand_product_quarks::QuarkGrandProductProof;

pub struct GrandProductLayerProof<const N: usize> {
    pub sumcheck_proof: SumcheckProof<N>,
    pub left_claim: RandomField<N>,
    pub right_claim: RandomField<N>,
}

pub struct GrandProductProof<const N: usize, S: BrakedownSpec> {
    pub gkr_layers: Vec<GrandProductLayerProof<N>>,
    pub quark_proof: Option<QuarkGrandProductProof<N, S>>,
}
