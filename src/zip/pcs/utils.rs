use ark_ff::Zero;
use ark_std::iterable::Iterable;
use i256::I256;

use crate::{
    field::RandomField as F,
    field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension as MLE_F,
    poly_z::mle::{build_eq_x_r as build_eq_x_r_z, DenseMultilinearExtension as MLE_Z},
    sumcheck::utils::build_eq_x_r as build_eq_x_r_f,
    zip::{utils::parallelize, Error},
};

fn err_too_many_variates(function: &str, upto: usize, got: usize) -> Error {
    Error::InvalidPcsParam(
        format!(
            "Too many variates of poly to {function} (param supports variates up to {upto} but got {got})"
        )
    )
}

// Ensures that polynomials and evaluation points are of appropriate size
pub(super) fn validate_input<'a>(
    function: &str,
    param_num_vars: usize,
    polys: impl Iterable<Item = &'a MLE_Z>,
    points: impl Iterable<Item = &'a Vec<i64>>,
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

pub(super) fn point_to_tensor_z(
    num_rows: usize,
    point: &[i64],
) -> Result<(Vec<i64>, Vec<i64>), Error> {
    assert!(num_rows.is_power_of_two());
    let (hi, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let t_0 = if !lo.is_empty() {
        build_eq_x_r_z(lo).unwrap()
    } else {
        MLE_Z::zero()
    };

    let t_1 = if !hi.is_empty() {
        build_eq_x_r_z(hi).unwrap()
    } else {
        MLE_Z::zero()
    };

    Ok((t_0.evaluations, t_1.evaluations))
}

pub(super) fn point_to_tensor_f<const N: usize>(
    num_rows: usize,
    point: &[F<N>],
    config: *const FieldConfig<N>,
) -> Result<(Vec<F<N>>, Vec<F<N>>), Error> {
    assert!(num_rows.is_power_of_two());
    let (hi, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let t_0 = if !lo.is_empty() {
        build_eq_x_r_f(lo, config).unwrap()
    } else {
        MLE_F::<N>::zero()
    };

    let t_1 = if !hi.is_empty() {
        build_eq_x_r_f(hi, config).unwrap()
    } else {
        MLE_F::<N>::zero()
    };

    Ok((t_0.evaluations, t_1.evaluations))
}

// Define function that performs a row operation on the evaluation matrix
// [t_0]^T * M]
pub(super) fn combine_rows_f<const N: usize>(
    coeffs: &[F<N>],
    evaluations: &[F<N>],
    row_len: usize,
) -> Vec<F<N>> {
    let mut combined_row = vec![F::zero(); row_len];
    parallelize(&mut combined_row, |(combined_row, offset)| {
        combined_row
            .iter_mut()
            .zip(offset..)
            .for_each(|(combined, column)| {
                *combined = F::zero();
                coeffs
                    .iter()
                    .zip(evaluations.iter().skip(column).step_by(row_len))
                    .for_each(|(coeff, eval)| {
                        *combined += &(*coeff * eval);
                    });
            })
    });

    combined_row
}
// Define function that performs a row operation on the evaluation matrix
// [t_0]^T * M]
pub(super) fn combine_rows_z(coeffs: &[i64], evaluations: &[i64], row_len: usize) -> Vec<I256> {
    let mut combined_row = vec![I256::from(0); row_len];
    parallelize(&mut combined_row, |(combined_row, offset)| {
        combined_row
            .iter_mut()
            .zip(offset..)
            .for_each(|(combined, column)| {
                *combined = I256::from(0);
                coeffs
                    .iter()
                    .zip(evaluations.iter().skip(column).step_by(row_len))
                    .for_each(|(coeff, eval)| {
                        *combined += &(I256::from(*coeff) * I256::from(*eval));
                    });
            })
    });

    combined_row
}
