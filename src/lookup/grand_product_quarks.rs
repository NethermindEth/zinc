use crate::{
    brakedown::pcs::structs::MultilinearBrakedownCommitment, field::RandomField as F, lookup::Math,
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

impl<const N: usize> QuarkGrandProduct<N> {
    /// The bottom/input layer of the grand products
    // (leaf values, batch size)
    type Leaves = (Vec<F<N>>, usize);
    type Config = QuarkGrandProductConfig;

    fn construct_with_config(leaves: Self::Leaves, config: Self::Config) -> Self {
        let (leaves, batch_size) = leaves;
        assert!(leaves.len() % batch_size == 0);
        assert!((leaves.len() / batch_size).is_power_of_two());

        let tree_depth = (leaves.len() / batch_size).log_2();
        let crossover = config.hybrid_layer_depth.get_crossover_depth();
        let num_layers = if tree_depth <= crossover {
            tree_depth - 1
        } else {
            crossover
        };

        // Taken 1 to 1 from the code in the BatchedDenseGrandProduct implementation
        let mut layers = Vec::<DenseInterleavedPolynomial<N>>::new();
        layers.push(DenseInterleavedPolynomial::new(leaves));

        for i in 0..num_layers {
            let previous_layer = &layers[i];
            let new_layer = previous_layer.layer_output();
            layers.push(new_layer);
        }

        // If the tree depth is too small we just do the GKR grand product
        if tree_depth <= num_layers {
            return Self {
                batch_size,
                quark_poly: None,
                base_layers: layers,
            };
        }

        // Take the top layer and then turn it into a quark poly
        // Note - We always push the base layer so the unwrap will work even with depth = 0
        let quark_poly = layers.pop().unwrap().coeffs;
        Self {
            batch_size,
            quark_poly: Some(quark_poly),
            base_layers: layers,
        }
    }

    fn num_layers(&self) -> usize {
        self.base_layers.len()
    }
}
