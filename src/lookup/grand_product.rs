use crate::{
    brakedown::code::BrakedownSpec,
    field::RandomField as F,
    lookup::Math,
    poly_f::polynomials::{
        dense_interleaved_polynomial::DenseInterleavedPolynomial, split_eq_poly::SplitEqPolynomial,
    },
    sumcheck::SumcheckProof,
    transcript::KeccakTranscript,
    zip::code::ZipSpec,
};

use super::grand_product_quarks::QuarkGrandProductProof;

pub struct GrandProductLayerProof<const N: usize> {
    pub sumcheck_proof: SumcheckProof<N>,
    pub left_claim: F<N>,
    pub right_claim: F<N>,
}

pub struct BatchedGrandProductProof<const N: usize, S: BrakedownSpec> {
    pub gkr_layers: Vec<GrandProductLayerProof<N>>,
    pub quark_proof: Option<QuarkGrandProductProof<N, S>>,
}

pub struct BatchedDenseGrandProduct<const N: usize> {
    layers: Vec<DenseInterleavedPolynomial<N>>,
}

pub trait BatchedGrandProduct<const N: usize, S: BrakedownSpec>: Sized {
    /// The bottom/input layer of the grand products
    type Leaves;
    type Config: Default + Clone + Copy;

    /// Constructs the grand product circuit(s) from `leaves` with the default configuration
    fn construct(leaves: Self::Leaves) -> Self {
        Self::construct_with_config(leaves, Self::Config::default())
    }
    /// Constructs the grand product circuit(s) from `leaves` with a config
    fn construct_with_config(leaves: Self::Leaves, config: Self::Config) -> Self;
    /// The number of layers in the grand product.
    fn num_layers(&self) -> usize;
    /// The claimed outputs of the grand products.
    fn claimed_outputs(&self) -> Vec<F<N>>;
    /// Returns an iterator over the layers of this batched grand product circuit.
    /// Each layer is mutable so that its polynomials can be bound over the course
    /// of proving.
    fn layers(&'_ mut self) -> impl Iterator<Item = &'_ mut dyn BatchedGrandProductLayer<N>>;

    /// Computes a batched grand product proof, layer by layer.

    fn prove_grand_product() -> (BatchedGrandProductProof<N, S>, Vec<F<N>>) {
        todo!()
    }
}

pub trait BatchedGrandProductLayer<const N: usize> {}

impl<const N: usize> BatchedDenseGrandProduct<N> {
    /// The bottom/input layer of the grand products
    // (leaf values, batch size)
    type Leaves = (Vec<F<N>>, usize);
    type Config = ();

    fn construct(leaves: Self::Leaves) -> Self {
        let (leaves, batch_size) = leaves;
        assert!(leaves.len() % batch_size == 0);
        assert!((leaves.len() / batch_size).is_power_of_two());

        let num_layers = (leaves.len() / batch_size).log_2();
        let mut layers: Vec<DenseInterleavedPolynomial<N>> = Vec::with_capacity(num_layers);
        layers.push(DenseInterleavedPolynomial::new(leaves));

        for i in 0..num_layers - 1 {
            let previous_layer = &layers[i];
            let new_layer = previous_layer.layer_output();
            layers.push(new_layer);
        }

        Self { layers }
    }
}
