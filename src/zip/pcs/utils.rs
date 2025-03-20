use ark_ff::Zero;
use ark_std::iterable::Iterable;

use crate::{
    field::RandomField as F,
    field_config::FieldConfig,
    poly_f::mle::DenseMultilinearExtension as MLE_F,
    poly_z::mle::{build_eq_x_r as build_eq_x_r_z, DenseMultilinearExtension as MLE_Z},
    sumcheck::utils::build_eq_x_r as build_eq_x_r_f,
    zip::Error,
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

#[cfg(test)]
mod tests {
    use crate::zip::utils::combine_rows;

    #[test]
    fn test_basic_combination() {
        let coeffs = vec![1, 2];
        let evaluations = vec![3, 4, 5, 6];
        let row_len = 2;

        let result = combine_rows(coeffs, evaluations, row_len);

        assert_eq!(result, vec![(3 + 2 * 5), (4 + 2 * 6)]);
    }

    #[test]
    fn test_second_combination() {
        let coeffs = vec![3, 4];
        let evaluations = vec![2, 4, 6, 8];
        let row_len = 2;

        let result = combine_rows(coeffs, evaluations, row_len);

        assert_eq!(result, vec![(3 * 2 + 4 * 6), (3 * 4 + 4 * 8)]);
    }
    #[test]
    fn test_large_values() {
        let coeffs = vec![1000, -500];
        let evaluations = vec![2000, -3000, 4000, -5000];
        let row_len = 2;

        let result = combine_rows(coeffs, evaluations, row_len);

        assert_eq!(
            result,
            vec![
                (1000 * 2000 + (-500) * 4000),
                (1000 * -3000 + (-500) * -5000)
            ]
        );
    }
}
