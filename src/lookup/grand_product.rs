use ark_ff::{One, Zero};
use ark_std::log2;
use itertools::Itertools;

use crate::{
    brakedown::code::BrakedownSpec,
    field::RandomField as F,
    field_config::FieldConfig,
    lookup::Math,
    poly_f::{
        mle::DenseMultilinearExtension,
        polynomials::{
            dense_interleaved_polynomial::DenseInterleavedPolynomial,
            split_eq_poly::SplitEqPolynomial,
        },
    },
    sumcheck::SumcheckProof,
    transcript::KeccakTranscript,
    zip::code::ZipSpec,
};

use super::grand_product_quarks::QuarkGrandProductProof;

pub struct BatchedGrandProductLayerProof<const N: usize> {
    pub sumcheck_proof: SumcheckProof<N>,
    pub left_claim: F<N>,
    pub right_claim: F<N>,
}

impl<const N: usize> BatchedGrandProductLayerProof<N> {
    pub fn verify(
        &self,
        claim: F<N>,
        num_rounds: usize,
        degree_bound: usize,
        transcript: &mut KeccakTranscript,
    ) -> (F<N>, Vec<F<N>>) {
        todo!()
    }
}

pub struct BatchedGrandProductProof<const N: usize, S: BrakedownSpec> {
    pub gkr_layers: Vec<BatchedGrandProductLayerProof<N>>,
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

    fn prove_grand_product(
        &mut self,
        transcript: &mut KeccakTranscript,
        config: *const FieldConfig<N>,
    ) -> (BatchedGrandProductProof<N, S>, Vec<F<N>>) {
        let mut proof_layers = Vec::with_capacity(self.num_layers());

        // Evaluate the MLE of the output layer at a random point to reduce the outputs to
        // a single claim.
        let mut outputs = self.claimed_outputs();
        transcript.absorb_slice(&outputs);
        outputs.resize(outputs.len().next_power_of_two(), F::zero());
        let output_mle = DenseMultilinearExtension::from_evaluations_slice(
            log2(outputs.len()) as usize,
            &outputs,
            config,
        );
        let mut r: Vec<F<N>> = transcript.get_challenges(output_mle.num_vars, config);
        let mut claim = output_mle
            .evaluate(&r, config)
            .expect("Evaluation of output has not worked!");

        for layer in self.layers() {
            proof_layers.push(layer.prove_layer(&mut claim, &mut r, transcript));
        }

        (
            BatchedGrandProductProof {
                gkr_layers: proof_layers,
                quark_proof: None,
            },
            r,
        )
    }
    fn verify_sumcheck_claim(
        layer_proofs: &[BatchedGrandProductLayerProof<N>],
        layer_index: usize,
        sumcheck_claim: F<N>,
        eq_eval: F<N>,
        grand_product_claim: &mut F<N>,
        r_grand_product: &mut Vec<F<N>>,
        transcript: &mut KeccakTranscript,
        config: *const FieldConfig<N>,
    ) {
        let layer_proof = &layer_proofs[layer_index];
        let expected_sumcheck_claim: F<N> =
            layer_proof.left_claim * layer_proof.right_claim * eq_eval;
        assert_eq!(expected_sumcheck_claim, sumcheck_claim);

        // produce a random challenge to condense two claims into a single claim
        let r_layer = transcript.get_challenge(config);
        *grand_product_claim =
            layer_proof.left_claim + r_layer * (layer_proof.right_claim - layer_proof.left_claim);

        r_grand_product.push(r_layer);
    }

    /// Function used for layer sumchecks in the generic batch verifier as well as the quark layered sumcheck hybrid
    fn verify_layers(
        proof_layers: &[BatchedGrandProductLayerProof<N>],
        mut claim: F<N>,
        transcript: &mut KeccakTranscript,
        r_start: Vec<F<N>>,
        config: *const FieldConfig<N>,
    ) -> (F<N>, Vec<F<N>>) {
        let mut r_grand_product = r_start.clone();
        let fixed_at_start = r_start.len();

        for (layer_index, layer_proof) in proof_layers.iter().enumerate() {
            let (sumcheck_claim, r_sumcheck) =
                layer_proof.verify(claim, layer_index + fixed_at_start, 3, transcript);

            transcript.absorb_random_field(&layer_proof.left_claim);
            transcript.absorb_random_field(&layer_proof.right_claim);

            let eq_eval: F<N> = r_grand_product
                .iter()
                .zip_eq(r_sumcheck.iter().rev())
                .map(|(&r_gp, &r_sc)| r_gp * r_sc + (F::one() - r_gp) * (F::one() - r_sc))
                .product();

            r_grand_product = r_sumcheck.into_iter().rev().collect();

            Self::verify_sumcheck_claim(
                proof_layers,
                layer_index,
                sumcheck_claim,
                eq_eval,
                &mut claim,
                &mut r_grand_product,
                transcript,
                config,
            );
        }

        (claim, r_grand_product)
    }
    /// Verifies the given grand product proof.
    fn verify_grand_product(
        proof: &BatchedGrandProductProof<N, S>,
        claimed_outputs: &[F<N>],
        transcript: &mut KeccakTranscript,
        config: *const FieldConfig<N>,
    ) -> (F<N>, Vec<F<N>>) {
        // Evaluate the MLE of the output layer at a random point to reduce the outputs to
        // a single claim.
        transcript.absorb_slice(claimed_outputs);
        let r: Vec<F<N>> =
            transcript.get_challenges(claimed_outputs.len().next_power_of_two().log_2(), config);
        let claim = DenseMultilinearExtension::from_evaluations_vec(
            r.len(),
            claimed_outputs.to_vec(),
            config,
        )
        .evaluate(&r, config)
        .expect("Evalution has not worked!");

        Self::verify_layers(&proof.gkr_layers, claim, transcript, r, config)
    }

    fn quark_poly(&self) -> Option<&[F<N>]> {
        None
    }
}

pub trait BatchedGrandProductLayer<const N: usize> {
    fn prove_layer(
        &mut self,
        claim: &mut F<N>,
        r_grand_product: &mut Vec<F<N>>,
        transcript: &mut KeccakTranscript,
    ) -> BatchedGrandProductLayerProof<N> {
        todo!()
    }
}

impl<const N: usize, S: BrakedownSpec> BatchedGrandProduct<N, S> for BatchedDenseGrandProduct<N> {
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

    fn construct_with_config(leaves: Self::Leaves, config: Self::Config) -> Self {
        todo!()
    }

    fn num_layers(&self) -> usize {
        todo!()
    }

    fn claimed_outputs(&self) -> Vec<F<N>> {
        todo!()
    }

    fn layers(&'_ mut self) -> impl Iterator<Item = &'_ mut dyn BatchedGrandProductLayer<N>> {
        std::iter::empty()
    }
}
