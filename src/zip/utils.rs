use ark_std::{
    cfg_iter_mut,
    ops::{Add, Mul},
    vec,
    vec::Vec,
};
use num_integer::Integer as NumInteger;
use rand::{SeedableRng, rngs::StdRng, seq::SliceRandom};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::traits::{Integer, Words};

pub(crate) fn inner_product<'a, 'b, T, L, R>(lhs: L, rhs: R) -> T
where
    T: Clone + Mul<Output = T> + Add<Output = T> + Default + 'a + 'b,
    L: IntoIterator<Item = &'a T>,
    R: IntoIterator<Item = &'b T>,
{
    lhs.into_iter()
        .zip(rhs)
        .map(|(lhs, rhs)| lhs.clone() * rhs.clone())
        .reduce(|acc, product| acc + product)
        .unwrap_or_default()
}

pub(crate) fn div_ceil(dividend: usize, divisor: usize) -> usize {
    NumInteger::div_ceil(&dividend, &divisor)
}

pub(crate) fn num_threads() -> usize {
    #[cfg(feature = "parallel")]
    return rayon::current_num_threads();

    #[cfg(not(feature = "parallel"))]
    return 1;
}

/// Computes a linear combination of multiple evaluation rows into a single combined row.
///
/// Given a flat `evaluations` vector, interpreted as a matrix with `row_len` columns,
/// this function treats each consecutive `row_len` values as one row. The output is a single row,
/// computed by multiplying each input row by the corresponding coefficient from `coeffs`,
/// and summing these scaled rows column-wise.
///
/// This is equivalent to performing a matrix-vector multiplication where the matrix is formed by
/// the evaluations and the vector is formed by the coefficients.
///
/// # Arguments
///
/// - `coeffs`: A slice of coefficients to be applied to each row.
/// - `evaluations`: A flattened slice of evaluations, arranged row-wise.
/// - `row_len`: The number of columns in each evaluation row.
///
/// # Returns
///
/// A vector of length `row_len` representing the combined row.
pub(super) fn combine_rows<F>(coeffs: &[F], evaluations: &[F], row_len: usize) -> Vec<F>
where
    F: Clone
        + Default
        + Send
        + Sync
        + for<'b> ark_std::ops::AddAssign<&'b F>
        + for<'b> ark_std::ops::Mul<&'b F, Output = F>,
{
    let mut combined_row = vec![F::default(); row_len];

    cfg_iter_mut!(combined_row)
        .enumerate()
        .for_each(|(column_index, combined_val)| {
            // Get an iterator over all values in the current column.
            let column_evals = evaluations.iter().skip(column_index).step_by(row_len);

            // Compute the dot product of the coefficients and the column values.
            let sum = coeffs
                .iter()
                .zip(column_evals)
                .map(|(coeff, eval)| coeff.clone() * eval)
                .fold(F::default(), |mut acc, val| {
                    acc += &val;
                    acc
                });

            *combined_val = sum;
        });

    combined_row
}

/// Computes a linear combination of multiple evaluation rows into a single combined row.
///
/// Given a flat `evaluations` vector, interpreted as a matrix with `row_len` columns,
/// this function treats each consecutive `row_len` values as one row. The output is a single row,
/// computed by multiplying each input row by the corresponding coefficient from `coeffs`,
/// and summing these scaled rows column-wise.
///
/// This is equivalent to performing a matrix-vector multiplication where the matrix is formed by
/// the evaluations and the vector is formed by the coefficients.
///
/// # Arguments
///
/// - `coeffs`: Coefficients applied to each row.
/// - `evaluations`: Flattened evaluations arranged row-wise.
/// - `row_len`: Number of columns per evaluation row.
///
/// # Returns
///
/// A vector of length `row_len` representing the combined row.
pub(super) fn combine_rows_as_iter<'a, F, C, E>(coeffs: C, evaluations: E, row_len: usize) -> Vec<F>
where
    F: Clone
        + Default
        + Send
        + Sync
        + for<'b> ark_std::ops::AddAssign<&'b F>
        + for<'b> ark_std::ops::Mul<&'b F, Output = F>,
    C: IntoIterator<Item = F> + Sync,
    E: IntoIterator<Item = F> + Sync,
    C::IntoIter: Clone + Send + Sync,
    E::IntoIter: Clone + Send + Sync,
{
    let coeffs_iter = coeffs.into_iter();
    let evaluations_iter = evaluations.into_iter();

    let mut combined_row = vec![F::default(); row_len];
    cfg_iter_mut!(combined_row)
        .enumerate()
        .for_each(|(offset, combined)| {
            *combined = F::default();
            coeffs_iter
                .clone()
                .zip(evaluations_iter.clone().skip(offset).step_by(row_len))
                .for_each(|(coeff, eval)| {
                    *combined += &(coeff * &eval);
                });
        });

    combined_row
}

pub(super) fn expand<N: Integer, M: Integer + for<'a> From<&'a N>>(narrow_int: &N) -> M {
    assert!(
        N::W::num_words() <= M::W::num_words(),
        "Cannot squeeze a wide integer into a narrow integer."
    );

    M::from(narrow_int)
}

/// Reorder the elements in slice using the given randomness seed
pub(super) fn shuffle_seeded<T>(slice: &mut [T], seed: u64) {
    let mut rng = StdRng::seed_from_u64(seed);
    slice.shuffle(&mut rng);
}

#[cfg(test)]
mod test {

    use crypto_bigint::Random;
    use num_traits::{ConstOne, ConstZero};

    use crate::{
        field::Int,
        traits::Integer,
        zip::utils::{expand, inner_product},
    };

    #[test]
    fn test_inner_product_basic() {
        let lhs = [1, 2, 3];
        let rhs = [4, 5, 6];
        assert_eq!(inner_product(lhs.iter(), rhs.iter()), 4 + 2 * 5 + 3 * 6);
    }

    #[test]
    fn test_expand_normal() {
        let input_words = [1u64, 2u64];
        let input = Int::<2>::from(input_words);
        let expanded = expand::<Int<2>, Int<4>>(&input);

        let expected_words = [1u64, 2u64, 0u64, 0u64];
        assert_eq!(expanded.as_words(), expected_words);
    }

    #[test]
    fn test_expand_identity() {
        let input_words = [42u64, 99u64];
        let input = Int::<2>::from(input_words);
        let expanded = expand::<Int<2>, Int<2>>(&input);

        let expected_words = [42u64, 99u64];
        assert_eq!(expanded.as_words(), expected_words);
    }

    #[test]
    #[should_panic(expected = "Cannot squeeze a wide integer into a narrow integer.")]
    fn test_expand_invalid() {
        let input = Int::<4>::from([1, 2, 3, 4]);
        // N = 4, M = 2 → should panic
        let _ = expand::<Int<4>, Int<2>>(&input);
    }

    #[test]
    fn test_expand_zero_padding() {
        let input = Int::<1>::from([123]);
        let expanded = expand::<Int<1>, Int<3>>(&input);

        let expected_words = [123u64, 0u64, 0u64];
        assert_eq!(expanded.as_words(), expected_words);
    }

    #[test]
    fn test_expand_all_zeros() {
        let input = Int::<2>::from([0u64, 0u64]);
        let expanded = expand::<Int<2>, Int<4>>(&input);

        let expected_words = [0u64, 0u64, 0u64, 0u64];
        assert_eq!(expanded.as_words(), expected_words);
    }
    #[test]
    fn test_expand_negative_number_identity() {
        // Example negative number in two's complement for 2 words
        let negative_val = Int::<2>::from([!0u64, !0u64]); // -1
        let expanded = expand::<Int<2>, Int<2>>(&negative_val);

        assert_eq!(expanded, Int::<2>::ZERO - &Int::<2>::ONE);
    }

    #[test]
    fn test_expand_negative_number_wider() {
        let mut rg = ark_std::test_rng();

        let mut positive_val = Int::<2>::random(&mut rg);
        if positive_val < Int::ZERO {
            positive_val = Int::<2>::ZERO - &positive_val;
        }

        let expanded_positive = expand::<Int<2>, Int<4>>(&positive_val);

        let negative_val = Int::<2>::ZERO - &positive_val;
        let expanded_negative = expand::<Int<2>, Int<4>>(&negative_val);

        let expected_negative = Int::<4>::ZERO - &expanded_positive;

        assert_eq!(expanded_negative, expected_negative);
    }
}
