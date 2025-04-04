use ark_ff::Zero;

use crate::{field::RandomField as F, poly_f::polynomials::util::unsafe_allocate_zero_vec};
use ark_std::slice::Chunks;
/// Represents a single layer of a grand product circuit.
///
/// A layer is assumed to be arranged in "interleaved" order, i.e. the natural
/// order in the visual representation of the circuit:
///      Λ        Λ        Λ        Λ
///     / \      / \      / \      /  \
///   L0   R0  L1   R1  L2   R2  L3   R3   <- This is layer would be represented as [L0, R0, L1, R1, L2, R2, L3, R3]
///                                           (as opposed to e.g. [L0, L1, L2, L3, R0, R1, R2, R3])
#[derive(Default, Debug, Clone)]
pub struct DenseInterleavedPolynomial<const N: usize> {
    /// The coefficients for the "left" and "right" polynomials comprising a
    /// dense grand product layer.
    /// The coefficients are in interleaved order:
    /// [L0, R0, L1, R1, L2, R2, L3, R3, ...]
    pub(crate) coeffs: Vec<F<N>>,
    /// The effective length of `coeffs`. When binding, we update this length
    /// instead of truncating `coeffs`, which incurs the cost of dropping the
    /// truncated values.
    len: usize,
    /// A reused buffer where bound values are written to during `bind`.
    /// With every bind, `coeffs` and `binding_scratch_space` are swapped.
    binding_scratch_space: Vec<F<N>>,
}
impl<const N: usize> DenseInterleavedPolynomial<N> {
    pub fn new(coeffs: Vec<F<N>>) -> Self {
        assert!(coeffs.len() % 2 == 0);
        let len = coeffs.len();
        Self {
            coeffs,
            len,
            binding_scratch_space: unsafe_allocate_zero_vec(len.next_multiple_of(4) / 2),
        }
    }
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn iter(&self) -> impl Iterator<Item = &F<N>> {
        self.coeffs[..self.len].iter()
    }

    pub fn chunks(&self, chunk_size: usize) -> Chunks<'_, F<N>> {
        self.coeffs[..self.len].chunks(chunk_size)
    }

    #[cfg(test)]
    pub fn interleave(left: &Vec<F<N>>, right: &Vec<F<N>>) -> Self {
        assert_eq!(left.len(), right.len());
        let mut interleaved = vec![];
        for i in 0..left.len() {
            interleaved.push(left[i]);
            interleaved.push(right[i]);
        }
        Self::new(interleaved)
    }

    pub fn uninterleave(&self) -> (Vec<F<N>>, Vec<F<N>>) {
        let left: Vec<F<N>> = self.coeffs[..self.len].iter().copied().step_by(2).collect();
        let mut right: Vec<F<N>> = self.coeffs[..self.len]
            .iter()
            .copied()
            .skip(1)
            .step_by(2)
            .collect();
        if right.len() < left.len() {
            right.resize(left.len(), F::zero());
        }
        (left, right)
    }

    pub fn layer_output(&self) -> Self {
        let output = self.chunks(2).map(|chunk| chunk[0] * chunk[1]).collect();
        Self::new(output)
    }
}
