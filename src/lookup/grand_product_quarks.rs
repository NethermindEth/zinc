use std::marker::PhantomData;

use ark_ff::{One, Zero};

use crate::{
    brakedown::{
        code::BrakedownSpec,
        pcs::structs::{
            MultilinearBrakedown, MultilinearBrakedownCommitment, MultilinearBrakedownParams,
        },
    },
    field::RandomField  as F,
    field_config::FieldConfig,
    lookup::Math,
    poly_f::{
        mle::{DenseMultilinearExtension, MultilinearExtension},
        polynomials::dense_interleaved_polynomial::DenseInterleavedPolynomial,
    },
    sumcheck::{utils::build_eq_x_r, MLSumcheck, SumcheckProof},
    transcript::KeccakTranscript,
};

use super::grand_product::GrandProductProof;

pub struct QuarkGrandProductProof<const N: usize, S: BrakedownSpec> {
    sumcheck_proof: SumcheckProof<N>,
    g_commitment: MultilinearBrakedownCommitment<N>,
    g_r_sumcheck: F<N>,
    g_r_prime: (F<N>, F<N>),
    v_r_prime: (F<N>, F<N>),
    pub num_vars: usize,
    _brakedown_marker: PhantomData<S>,
}
pub struct QuarkGrandProduct<const N: usize, S: BrakedownSpec> {
    batch_size: usize,
    quark_poly: Option<Vec<F<N>>>,
    base_layers: Vec<DenseInterleavedPolynomial<N>>,
    _brakedown_marker: PhantomData<S>,
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

impl<const N: usize, S: BrakedownSpec> QuarkGrandProduct<N, S> {
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
                _brakedown_marker: PhantomData,
            };
        }

        // Take the top layer and then turn it into a quark poly
        // Note - We always push the base layer so the unwrap will work even with depth = 0
        let quark_poly = layers.pop().unwrap().coeffs;
        Self {
            batch_size,
            quark_poly: Some(quark_poly),
            base_layers: layers,
            _brakedown_marker: PhantomData,
        }
    }

    fn num_layers(&self) -> usize {
        self.base_layers.len()
    }

    /// The claimed outputs of the grand products.
    fn claimed_outputs(&self) -> Vec<F<N>> {
        if let Some(quark_poly) = &self.quark_poly {
            let chunk_size = quark_poly.len() / self.batch_size;
            quark_poly
                .chunks(chunk_size)
                .map(|chunk| chunk.iter().product())
                .collect()
        } else {
            let top_layer = &self.base_layers[self.base_layers.len() - 1];
            top_layer
                .chunks(2)
                .map(|chunk| chunk[0] * chunk[1])
                .collect()
        }
    }

    /// Returns an iterator over the layers of this batched grand product circuit.
    /// Each layer is mutable so that its polynomials can be bound over the course
    /// of proving.
    fn layers(&'_ mut self) -> impl Iterator<Item = &'_ mut DenseInterleavedPolynomial<N>> {
        self.base_layers.iter_mut().map(|layer| layer).rev()
    }

    fn quark_poly(&self) -> Option<&[F<N>]> {
        self.quark_poly.as_deref()
    }

    fn prove_grand_product(&mut self) -> (GrandProductProof<N, S>, Vec<F<N>>) {
        todo!()
    }

    fn verify_grand_product(
        &self,
        proof: &GrandProductProof<N, S>,
        claimed_outputs: &[F<N>],
    ) -> (F<N>, Vec<F<N>>) {
        todo!()
    }
}

pub struct QuarkGrandProductBase<const N: usize, S: BrakedownSpec> {
    _brakedown_marker: PhantomData<S>,
}

impl<const N: usize, S: BrakedownSpec> QuarkGrandProductBase<N, S> {
    pub fn prove_quark_grand_product(
        grand_product: &QuarkGrandProduct<N, S>,
    ) -> (GrandProductProof<N, S>, Vec<F<N>>) {
        todo!()
    }

    pub fn verify_quark_grand_product(
        proof: &GrandProductProof<N, S>,
        claimed_outputs: &[F<N>],
    ) -> (GrandProductProof<N, S>, Vec<F<N>>) {
        todo!()
    }
}

impl<const N: usize, S: BrakedownSpec> QuarkGrandProductProof<N, S> {
    /// Computes a grand product proof using Section 5 technique from Quarks Paper
    /// First - Extends the evals of v to create an f poly, then commits to it and evals
    /// Then - Constructs a g poly and preforms sumcheck proof that sum == 0
    /// Finally - computes opening proofs for a random sampled during sumcheck proof and returns
    /// Returns a random point and evaluation to be verified by the caller (which our hybrid prover does with GKR)
    fn prove(
        v: &[F<N>],
        r_outputs: &[F<N>],
        claim: F<N>,
        transcript: &mut KeccakTranscript,
        params: &MultilinearBrakedownParams<N>,
        config: *const FieldConfig<N>,
    ) -> (Self, Vec<F<N>>, F<N>) {
        let v_length = v.len();
        let v_variables = v_length.log_2();

        let (g_polynomial, f_x_0, f_x_1) = v_into_f_s(v, config);
        let mut sumcheck_polys = vec![g_polynomial.clone(), f_x_0, f_x_1];

        // Commit to g(x) = f(1, x)
        let g_commitment = MultilinearBrakedown::<N, S>::commit(params, &g_polynomial).unwrap();
        transcript.absorb(g_commitment.root()); // TODO: add intermediate hashes?

        let tau = transcript.get_challenges(v_variables, config);
        let eq_tau = build_eq_x_r(&tau, config).unwrap();

        sumcheck_polys.push(eq_tau);

        // This is where things start to deviate from the protocol described in
        // Quarks Section 5.
        //
        // We batch our grand products by laying out the circuits side-by-side, and
        // proving them together as one big circuit with k outputs, where k is the batch size
        // In `prove_grand_product`, we evaluate the MLE of these outputs at random point,
        //   claim := /tilde{outputs}(r_outputs)
        //
        // Quarks Section 5 assumes there's only one output, P = f(1, ..., 1, 0).
        // But claim != f(1, ..., 1, 0), so we have to use a different sumcheck expression
        //
        // If you closely examin `v_into_f_s` and work it out, you'll find that  our k grand product
        // outputs are contained in f(1, x) at x = (1, ..., 1, 0, b), where b /in {0, 1}^{log2(k)}.
        // So we have:
        //   claim = \tilde{outputs}(r_outputs)
        //         = \sum_b EQ(r_outputs, b) * outputs(b)
        //         = \sum_x EQ(1, ..., 1, 0, r_outputs, x) * f(1, x)        where r_outputs ∈ 𝔽^{log2(k)}, x ∈ {0, 1}^{log2(kn)}
        //
        // Modifying the sumcheck instance described in Section  of the Quarks aper, we will
        // be proving:
        //   claim = \sum_x (EQ(\tau, x) * (f(1, x) - f(x, 0) * f(x, 1)) + EQ(1, ..., 1, 0, r_outputs, x) * f(1, x))
        //
        // Note that the first half of the summand EQ(\tau, x) * (f(1, x) - f(x, 0) * f(x, 1))
        // should equal 0 for all x ∈ {0, 1}^{log2(kn)}, ensuring that every output value f(1, x) is equal to the
        // product of its input values f(x, 0) and f(x, 1).

        // First we compute EQ(1, ..., 1, 0, r_outputs, x)
        let mut one_padded_r_outputs = vec![F::<N>::one(); v_variables];
        let slice_index = one_padded_r_outputs.len() - r_outputs.len();
        one_padded_r_outputs[slice_index..].copy_from_slice(r_outputs);
        one_padded_r_outputs[slice_index - 1] = F::<N>::zero();
        let eq_output = build_eq_x_r(&one_padded_r_outputs, config).unwrap(); // TODO: Check this

        // We add eq_output as the last polynomial in the sumcheck
        sumcheck_polys.push(eq_output);

        // This is the sumcheck polynomial
        // EQ(\tau, x) * (f(1, x) - f(x, 0) * f(x, 1)) + EQ(1, ..., 1, 0, r_outputs, x) * f(1, x)
        let output_check_fn = |vals: &[F<N>]| -> F<N> {
            let f_1x = vals[0];
            let f_x0 = vals[1];
            let f_x1 = vals[2];
            let eq_tau = vals[3];
            let eq_output = vals[4];

            eq_tau * (f_1x - f_x0 * f_x1) + eq_output * f_1x
        };

        // Now run the sumcheck
        let (sum_check_proof, prover_state) = MLSumcheck::prove_as_subprotocol(
            transcript,
            sumcheck_polys,
            v_variables,
            3,
            output_check_fn,
            config,
        );
        let r_sumcheck = prover_state.randomness;

        todo!()
    }
}

/// Computes the polyynomials f(1, x) f(x, 0), and f(x, 1) from the v polynomial,
/// as described in Lemma 5.1 of the Quarks paper.
fn v_into_f_s<const N: usize>(
    v_evals: &[F<N>],
    config: *const FieldConfig<N>,
) -> (
    DenseMultilinearExtension<N>,
    DenseMultilinearExtension<N>,
    DenseMultilinearExtension<N>,
) {
    let v_length = v_evals.len();
    assert!(v_length.is_power_of_two(), "Is this already assured?"); // TODO: Check this
    let v_variables = v_length.log_2();
    let mut f_evals = vec![F::<N>::zero(); 2 * v_length];
    f_evals[..v_length].copy_from_slice(v_evals);

    for i in v_length..2 * v_length {
        let i_shift_mod = (i << 1) % (2 * v_length);
        // The transform follows the logic of the paper and to accumulate
        // the partial sums into the correct indices.
        f_evals[i] = f_evals[i_shift_mod] * f_evals[i_shift_mod + 1];
    }

    // We pull out the coefficient which instantiate the lower d polys for the sumcheck
    let f_1_x = f_evals[v_length..].to_vec();

    let mut f_x_0: Vec<F<N>> = Vec::new();
    let mut f_x_1: Vec<F<N>> = Vec::new();
    for (i, &x) in f_evals.iter().enumerate() {
        if i % 2 == 0 {
            f_x_0.push(x);
        } else {
            f_x_1.push(x);
        }
    }

    (
        DenseMultilinearExtension::from_evaluations_vec(v_variables, f_x_0, config),
        DenseMultilinearExtension::from_evaluations_vec(v_variables, f_x_1, config),
        DenseMultilinearExtension::from_evaluations_vec(v_variables, f_1_x, config),
    )
}
