use crate::field::RandomField as F;
use ark_ff::Zero;
use itertools::Itertools;
use num_integer::Integer;

pub fn horner_new<const N: usize>(coeffs: &[F<N>], x: &F<N>) -> F<N> {
    let coeff_vec: Vec<&F<N>> = coeffs.iter().rev().collect();
    let mut acc = F::zero();
    for c in coeff_vec {
        acc = acc * x + *c;
    }
    acc
    //2
    //.fold(F::ZERO, |acc, coeff| acc * x + coeff)
}

pub fn horner<const N: usize>(coeffs: &[F<N>], x: &F<N>) -> F<N> {
    let coeff_vec: Vec<&F<N>> = coeffs.iter().rev().collect();
    let mut acc = F::zero();
    for c in coeff_vec {
        acc = acc * x + *c;
    }
    acc
    //2
    //.fold(F::ZERO, |acc, coeff| acc * x + coeff)
}

pub fn horner_orig<const N: usize>(coeffs: &[F<N>], x: &F<N>) -> F<N> {
    coeffs
        .iter()
        .rev()
        .fold(F::zero(), |acc, coeff| acc * x + *coeff)
}

pub fn inner_product<'a, 'b, const N: usize>(
    lhs: impl IntoIterator<Item = &'a F<N>>,
    rhs: impl IntoIterator<Item = &'b F<N>>,
) -> F<N> {
    lhs.into_iter()
        .zip_eq(rhs)
        .map(|(lhs, rhs)| *lhs * rhs)
        .reduce(|acc, product| acc + product)
        .unwrap_or_default()
}

pub fn div_rem(dividend: usize, divisor: usize) -> (usize, usize) {
    Integer::div_rem(&dividend, &divisor)
}

pub fn div_ceil(dividend: usize, divisor: usize) -> usize {
    Integer::div_ceil(&dividend, &divisor)
}

pub fn num_threads() -> usize {
    #[cfg(feature = "parallel")]
    {
        let nt = rayon::current_num_threads();
        return nt;
    }

    #[cfg(not(feature = "parallel"))]
    return 1;
}

pub fn parallelize_iter<I, T, F>(iter: I, f: F)
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

pub fn parallelize<T, F>(v: &mut [T], f: F)
where
    T: Send,
    F: Fn((&mut [T], usize)) + Send + Sync + Clone,
{
    #[cfg(feature = "parallel")]
    {
        use crate::util::arithmetic::div_ceil;
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

pub fn par_sort_unstable<T>(v: &mut [T])
where
    T: Ord + Send,
{
    #[cfg(feature = "parallel")]
    {
        use rayon::slice::ParallelSliceMut;
        v.par_sort_unstable();
    }

    #[cfg(not(feature = "parallel"))]
    v.sort_unstable();
}

#[cfg(feature = "parallel")]
pub fn par_map_collect<T, R, C>(
    v: impl rayon::prelude::IntoParallelIterator<Item = T>,
    f: impl Fn(T) -> R + Send + Sync,
) -> C
where
    T: Send + Sync,
    R: Send,
    C: rayon::prelude::FromParallelIterator<R>,
{
    use rayon::prelude::ParallelIterator;
    v.into_par_iter().map(f).collect()
}

#[cfg(not(feature = "parallel"))]
pub fn par_map_collect<T, R, C>(v: impl IntoIterator<Item = T>, f: impl Fn(T) -> R) -> C
where
    C: FromIterator<R>,
{
    v.into_iter().map(f).collect()
}
