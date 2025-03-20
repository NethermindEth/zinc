use std::ops::{Add, Mul};

use num_integer::Integer;

pub(crate) fn inner_product<'a, 'b, T, L, R>(lhs: L, rhs: R) -> T
where
    T: Copy + Mul<Output = T> + Add<Output = T> + Default + 'a + 'b,
    L: IntoIterator<Item = &'a T>,
    R: IntoIterator<Item = &'b T>,
{
    lhs.into_iter()
        .zip(rhs)
        .map(|(lhs, rhs)| *lhs * *rhs)
        .reduce(|acc, product| acc + product)
        .unwrap_or_default()
}

pub(crate) fn div_ceil(dividend: usize, divisor: usize) -> usize {
    Integer::div_ceil(&dividend, &divisor)
}

pub(crate) fn num_threads() -> usize {
    #[cfg(feature = "parallel")]
    return rayon::current_num_threads();

    #[cfg(not(feature = "parallel"))]
    return 1;
}

pub(crate) fn parallelize_iter<I, T, F>(iter: I, f: F)
where
    I: Send + Iterator<Item = T>,
    T: Send,
    F: Fn(T) + Send + Sync + Clone,
{
    #[cfg(feature = "parallel")]
    rayon::scope(|scope| {
        iter.for_each(|item| {
            let f = &f;
            scope.spawn(move |_| f(item))
        })
    });

    #[cfg(not(feature = "parallel"))]
    iter.for_each(f);
}

pub(crate) fn parallelize<T, F>(v: &mut [T], f: F)
where
    T: Send,
    F: Fn((&mut [T], usize)) + Send + Sync + Clone,
{
    #[cfg(feature = "parallel")]
    {
        use crate::zip::utils::div_ceil;
        let num_threads = num_threads();
        let chunk_size = div_ceil(v.len(), num_threads);
        if chunk_size < num_threads {
            f((v, 0));
        } else {
            parallelize_iter(v.chunks_mut(chunk_size).zip((0..).step_by(chunk_size)), f);
        }
    }

    #[cfg(not(feature = "parallel"))]
    f((v, 0));
}

// Define function that performs a row operation on the evaluation matrix
// [t_0]^T * M]
pub(super) fn combine_rows<'a, F, C, E>(coeffs: C, evaluations: E, row_len: usize) -> Vec<F>
where
    F: Copy
        + Default
        + Send
        + Sync
        + for<'b> std::ops::AddAssign<&'b F>
        + for<'b> std::ops::Mul<&'b F, Output = F>,
    C: IntoIterator<Item = F> + Sync,
    E: IntoIterator<Item = F> + Sync,
    C::IntoIter: Clone + Send + Sync,
    E::IntoIter: Clone + Send + Sync,
{
    let coeffs_iter = coeffs.into_iter();
    let evaluations_iter = evaluations.into_iter();

    let mut combined_row = vec![F::default(); row_len];
    parallelize(&mut combined_row, |(combined_row, offset)| {
        combined_row
            .iter_mut()
            .zip(offset..)
            .for_each(|(mut combined, column)| {
                *combined = F::default();
                coeffs_iter
                    .clone()
                    .zip(evaluations_iter.clone().skip(column).step_by(row_len))
                    .for_each(|(coeff, eval)| {
                        *combined += &(coeff * &eval);
                    });
            })
    });

    combined_row
}

#[cfg(test)]
mod test {
    use crate::zip::utils::inner_product;

    #[test]
    fn test_inner_product_basic() {
        let lhs = [1, 2, 3];
        let rhs = [4, 5, 6];
        assert_eq!(inner_product(lhs.iter(), rhs.iter()), 4 + 2 * 5 + 3 * 6);
    }
}
