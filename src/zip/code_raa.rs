use ark_std::{fmt::Debug, marker::PhantomData, ops::AddAssign, vec::Vec};
use num_traits::Zero;

use crate::{
    traits::{Field, FieldMap, Integer, Words, ZipTypes},
    zip::{
        code::{LinearCode, LinearCodeSpec},
        pcs::structs::ZipTranscript,
        utils::shuffle_seeded,
    },
};

/// Implementation of a repeat-accumulate-accumulate (RAA) codes over the binary field,
/// as defined by the Blaze paper (https://eprint.iacr.org/2024/1609)
#[derive(Debug, Clone)]
pub struct RaaCode<ZT: ZipTypes> {
    row_len: usize,

    repetition_factor: usize,

    num_column_opening: usize,

    num_proximity_testing: usize,

    /// Randomness seed for the first permutation
    perm_1_seed: u64,

    /// Randomness seed for the second permutation
    perm_2_seed: u64,

    phantom: PhantomData<ZT>,
}

impl<ZT: ZipTypes> RaaCode<ZT> {
    pub fn new<S: LinearCodeSpec, T: ZipTranscript<ZT::L>>(
        spec: &S,
        poly_size: usize,
        transcript: &mut T,
    ) -> Self {
        // Taken from original Zip codes

        let num_vars = poly_size.ilog2() as usize;
        let row_len = ((1 << num_vars) as u64).isqrt().next_power_of_two() as usize;
        let repetition_factor = spec.repetition_factor();

        let num_column_opening = spec.num_column_opening();
        let log2_q = <ZT::N as Integer>::W::num_words();
        let n_0 = 20.min((1 << num_vars) - 1);
        let num_proximity_testing = spec.num_proximity_testing(log2_q, row_len, n_0);

        // Width of each entry in codeword vector, in bits.
        // For RAA it's initial_bits + 2*log(repetition_factor) + num_variables
        let codeword_width_bits = {
            let initial_bits = ZT::N::num_bits();
            let rep_factor_log: usize = repetition_factor
                .checked_next_power_of_two()
                .expect("Repetition factor is too large")
                .ilog2()
                .try_into()
                .expect("Repetition factor logarithm is too large");
            let num_vars_even = if num_vars.is_multiple_of(2) {
                num_vars
            } else {
                num_vars + 1
            };
            initial_bits + num_vars_even + (2 * rep_factor_log)
        };
        assert!(
            ZT::K::num_bits() >= codeword_width_bits,
            "Cannot fit {codeword_width_bits}-bit wide codeword entries in {} bits integers",
            ZT::K::num_bits()
        );

        let perm_1_seed = transcript.get_u64();
        let perm_2_seed = transcript.get_u64();

        Self {
            row_len,
            repetition_factor,
            num_column_opening,
            num_proximity_testing,
            perm_1_seed,
            perm_2_seed,
            phantom: PhantomData,
        }
    }

    /// Do the actual encoding, as per RAA spec
    fn encode_inner<In, Out>(&self, row: &[In]) -> Vec<Out>
    where
        Out: Zero + AddAssign<Out> + for<'a> From<&'a In> + Clone,
    {
        debug_assert_eq!(
            row.len(),
            self.row_len,
            "Row length must match the code's row length"
        );
        let mut result: Vec<Out> = repeat(row, self.repetition_factor);
        shuffle_seeded(&mut result, self.perm_1_seed);
        accumulate(&mut result);
        shuffle_seeded(&mut result, self.perm_2_seed);
        accumulate(&mut result);
        debug_assert_eq!(result.len(), self.codeword_len());
        result
    }
}

impl<ZT: ZipTypes> LinearCode<ZT> for RaaCode<ZT> {
    fn row_len(&self) -> usize {
        self.row_len
    }

    fn codeword_len(&self) -> usize {
        self.row_len * self.repetition_factor
    }

    fn num_column_opening(&self) -> usize {
        self.num_column_opening
    }

    fn num_proximity_testing(&self) -> usize {
        self.num_proximity_testing
    }

    fn encode_wide<In, Out>(&self, row: &[In]) -> Vec<Out>
    where
        In: Integer,
        Out: Integer + for<'a> From<&'a In> + for<'a> From<&'a ZT::L>,
    {
        self.encode_inner(row)
    }

    fn encode_f<F: Field>(&self, row: &[F], _field: F::R) -> Vec<F>
    where
        ZT::L: FieldMap<F, Output = F>,
    {
        self.encode_inner(row)
    }
}

/// Repeat the given slice N times, e.g `[1,2,3] => [1,2,3,1,2,3]`
fn repeat<In, Out: for<'a> From<&'a In> + Clone>(
    input: &[In],
    repetition_factor: usize,
) -> Vec<Out> {
    input
        .iter()
        .map(|i| Out::from(i))
        .cycle()
        .take(input.len() * repetition_factor)
        .collect()
}

/// Perform an operation equivalent to multiplying the slice in-place by the accumulation matrix
/// from the RAA code - a lower triangular matrix of the appropriate size, i.e. a matrix looking
/// like this:
///
/// ```text
/// 1 0 0 0
/// 1 1 0 0
/// 1 1 1 0
/// 1 1 1 1
/// ```
fn accumulate<I>(input: &mut [I])
where
    I: Zero + AddAssign<I> + Clone,
{
    for i in 1..input.len() {
        input[i] += input[i - 1].clone();
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{vec, vec::Vec};
    use num_traits::Zero;

    use crate::{
        define_random_field_zip_types,
        field::Int,
        implement_random_field_zip_types,
        traits::ZipTypes,
        zip::{
            code::{DefaultLinearCodeSpec, LinearCode},
            code_raa::{accumulate, repeat, RaaCode},
            pcs::tests::MockTranscript,
            utils::shuffle_seeded,
        },
    };

    // Define common types for testing
    const INT_LIMBS: usize = 1;
    define_random_field_zip_types!();
    implement_random_field_zip_types!(INT_LIMBS);
    type ZT = RandomFieldZipTypes<INT_LIMBS>;
    type I = Int<INT_LIMBS>;

    #[test]
    fn repeat_function_duplicates_row_correctly() {
        let input = vec![Int::<INT_LIMBS>::from(10), Int::<INT_LIMBS>::from(20)];

        let repetition_factor = 3;

        let repeated_output = repeat::<_, I>(&input, repetition_factor);

        let expected_output: Vec<_> = [10, 20, 10, 20, 10, 20]
            .into_iter()
            .map(Int::<INT_LIMBS>::from)
            .collect();
        assert_eq!(
            repeated_output, expected_output,
            "Failed on repetition factor > 1"
        );

        let empty_input: Vec<I> = vec![];
        let repeated_empty = repeat::<_, I>(&empty_input, 5);
        assert!(repeated_empty.is_empty(), "Failed on empty input vector");

        let repeated_once = repeat::<_, I>(&input, 1);
        assert_eq!(repeated_once, input, "Failed on repetition factor of 1");
    }

    #[test]
    fn accumulate_function_computes_cumulative_sum() {
        let mut input1: Vec<I> = [1, 2, 3, 4].into_iter().map(I::from).collect();
        let expected1: Vec<I> = [1, 3, 6, 10].into_iter().map(I::from).collect();
        accumulate(&mut input1);
        assert_eq!(input1, expected1, "Failed on positive integers");

        let mut input2: Vec<I> = [5, 0, 2, 0].into_iter().map(I::from).collect();
        let expected2: Vec<I> = [5, 5, 7, 7].into_iter().map(I::from).collect();
        accumulate(&mut input2);
        assert_eq!(input2, expected2, "Failed on vector with zeros");

        let mut input3: Vec<I> = [-1, 5, -10, 2].into_iter().map(I::from).collect();
        let expected3: Vec<I> = [-1, 4, -6, -4].into_iter().map(I::from).collect();
        accumulate(&mut input3);
        assert_eq!(input3, expected3, "Failed on vector with negative numbers");

        let mut empty_input: Vec<I> = vec![];
        let expected_empty: Vec<I> = vec![];
        accumulate(&mut empty_input);
        assert_eq!(empty_input, expected_empty, "Failed on empty vector");
    }

    #[test]
    fn shuffle_is_deterministic_for_a_given_seed() {
        let original: Vec<I> = (1..=10).map(I::from).collect();
        let mut vec1 = original.clone();
        let mut vec2 = original.clone();
        let mut vec3 = original.clone();

        let seed1 = 12345;
        let seed2 = 54321;

        shuffle_seeded(&mut vec1, seed1);
        shuffle_seeded(&mut vec2, seed1);
        shuffle_seeded(&mut vec3, seed2);

        assert_eq!(
            vec1, vec2,
            "Shuffling with the same seed should produce the same result"
        );
        assert_ne!(
            vec1, vec3,
            "Shuffling with different seeds should produce different results"
        );
        assert_ne!(
            vec1, original,
            "Shuffled vector should not be the same as the original"
        );
        assert_ne!(
            vec3, original,
            "Shuffled vector should not be the same as the original"
        );
    }

    #[test]
    fn encoding_preserves_linearity() {
        let mut transcript = MockTranscript::default();
        let code = RaaCode::<ZT>::new(&DefaultLinearCodeSpec, 16, &mut transcript);

        let a: Vec<I> = (1..=4).map(I::from).collect();
        let b: Vec<I> = (5..=8).map(I::from).collect();
        let sum_ab: Vec<I> = a.iter().zip(b.iter()).map(|(x, y)| *x + y).collect();

        let encode_a: Vec<<ZT as ZipTypes>::M> = code.encode(&a);
        let encode_b: Vec<<ZT as ZipTypes>::M> = code.encode(&b);
        let encode_sum_ab: Vec<<ZT as ZipTypes>::M> = code.encode(&sum_ab);

        let sum_encode_ab: Vec<<ZT as ZipTypes>::M> = encode_a
            .iter()
            .zip(encode_b.iter())
            .map(|(x, y)| *x + y)
            .collect();

        assert_eq!(encode_sum_ab, sum_encode_ab);
    }

    #[test]
    fn encoding_zero_vector_results_in_zero_codeword() {
        let mut transcript = MockTranscript::default();
        let code = RaaCode::<ZT>::new(&DefaultLinearCodeSpec, 16, &mut transcript);

        let zero_vector: Vec<I> = vec![I::zero(); code.row_len()];
        let encoded_vector: Vec<<ZT as ZipTypes>::M> = code.encode(&zero_vector);

        let expected_codeword: Vec<<ZT as ZipTypes>::M> =
            vec![<ZT as ZipTypes>::M::zero(); code.codeword_len()];

        assert_eq!(
            encoded_vector, expected_codeword,
            "Encoding a zero vector should result in a zero codeword"
        );
    }

    #[test]
    #[should_panic(expected = "Cannot fit 96-bit wide codeword entries in 64 bits integers")]
    fn constructor_panics_on_insufficient_codeword_width() {
        #[derive(Debug, Clone)]
        struct MismatchedZipTypes;
        impl ZipTypes for MismatchedZipTypes {
            type N = Int<1>; // 64 bits
            type L = Int<2>; // 128 bits
            type K = Int<1>; // 64 bits - INSUFFICIENT
            type M = Int<4>; // 256 bits
        }

        let mut transcript = MockTranscript::default();
        let _code =
            RaaCode::<MismatchedZipTypes>::new(&DefaultLinearCodeSpec, 1 << 30, &mut transcript);
    }

    #[test]
    #[should_panic(expected = "Row length must match the code's row length")]
    #[cfg(debug_assertions)]
    fn encode_panics_on_mismatched_row_length() {
        let mut transcript = MockTranscript::default();
        let code = RaaCode::<ZT>::new(&DefaultLinearCodeSpec, 16, &mut transcript);
        let incorrect_row = vec![I::from(1), I::from(2), I::from(3)];
        let _: Vec<<ZT as ZipTypes>::M> = code.encode(&incorrect_row);
    }
}
