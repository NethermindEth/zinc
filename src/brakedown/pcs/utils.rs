use ark_ff::Zero;
use ark_std::iterable::Iterable;

use crate::{
    brakedown::Error, field::RandomField as F, field_config::FieldConfig,
    poly::mle::DenseMultilinearExtension, sumcheck::utils::build_eq_x_r,
};

fn err_too_many_variates(function: &str, upto: usize, got: usize) -> Error {
    Error::InvalidPcsParam(
        format!(
            "Too many variates of poly to {function} (param supports variates up to {upto} but got {got})"
        )
    )
}

// Ensures that polynomials and evaluation points are of appropriate size
pub(super) fn validate_input<'a, const N: usize>(
    function: &str,
    param_num_vars: usize,
    polys: impl Iterable<Item = &'a DenseMultilinearExtension<N>>,
    points: impl Iterable<Item = &'a Vec<F<N>>>,
) -> Result<(), Error> {
    // Ensure all the number of variables in the polynomials don't exceed the limit
    for poly in polys.iter() {
        if param_num_vars < poly.num_vars {
            return Err(err_too_many_variates(
                function,
                param_num_vars,
                poly.num_vars,
            ));
        }
    }

    // Ensure all the points are of correct length
    let input_num_vars = polys
        .iter()
        .map(|poly| poly.num_vars)
        .chain(points.iter().map(|point| point.len()))
        .next()
        .expect("To have at least 1 poly or point");

    for point in points.iter() {
        if point.len() != input_num_vars {
            return Err(Error::InvalidPcsParam(format!(
                "Invalid point (expect point to have {input_num_vars} variates but got {})",
                point.len()
            )));
        }
    }
    Ok(())
}
pub(super) fn point_to_tensor<const N: usize>(
    num_rows: usize,
    point: &[F<N>],
    config: *const FieldConfig<N>,
) -> Result<(Vec<F<N>>, Vec<F<N>>), Error> {
    assert!(num_rows.is_power_of_two());
    let (hi, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let t_0 = if !lo.is_empty() {
        build_eq_x_r(lo, config).unwrap()
    } else {
        DenseMultilinearExtension::<N>::zero()
    };

    let t_1 = if !hi.is_empty() {
        build_eq_x_r(hi, config).unwrap()
    } else {
        DenseMultilinearExtension::<N>::zero()
    };

    Ok((t_0.evaluations, t_1.evaluations))
}
