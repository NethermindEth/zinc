use std::borrow::Cow;
use std::{marker::PhantomData, slice};

use ark_ff::Zero;
use ark_std::iterable::Iterable;
use ark_std::rand::RngCore;
use itertools::Itertools;
use sha3::Digest;
use sha3::{digest::Output, Keccak256};

use crate::brakedown::code::LinearCodes;
use crate::brakedown::pcs_transcript::PcsTranscript;
use crate::brakedown::utils::inner_product;
use crate::field_config::FieldConfig;
use crate::sumcheck::utils::build_eq_x_r;
use crate::{field::RandomField as F, poly::mle::DenseMultilinearExtension};

use super::{
    code::{Brakedown, BrakedownSpec},
    utils::{div_ceil, num_threads, parallelize, parallelize_iter},
    Error,
};

#[derive(Debug)]
pub struct MultilinearBrakedown<const N: usize, S: BrakedownSpec>(PhantomData<S>);

impl<const N: usize, S: BrakedownSpec> Clone for MultilinearBrakedown<N, S> {
    fn clone(&self) -> Self {
        Self(PhantomData)
    }
}

#[derive(Clone, Debug)]
pub struct MultilinearBrakedownParams<const N: usize> {
    num_vars: usize,
    num_rows: usize,
    brakedown: Brakedown<N>,
}

impl<const N: usize> MultilinearBrakedownParams<N> {
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    pub fn brakedown(&self) -> &Brakedown<N> {
        &self.brakedown
    }
}

/// Representantation of a brakedown commitment to a multilinear polynomial
#[derive(Clone, Debug, Default)]
pub struct MultilinearBrakedownCommitment<const N: usize> {
    /// The encoded rows of the polynomial matrix representation
    rows: Vec<F<N>>,
    /// Hashes of the merkle tree with the encoded columns as leaves
    intermediate_hashes: Vec<Output<Keccak256>>,
    /// Root of the merkle tree with the encoded columns as leaves
    root: Output<Keccak256>,
}

impl<const N: usize> MultilinearBrakedownCommitment<N> {
    pub fn from_root(root: Output<Keccak256>) -> Self {
        Self {
            root,
            ..Default::default()
        }
    }

    pub fn rows(&self) -> &[F<N>] {
        &self.rows
    }

    pub fn intermediate_hashes(&self) -> &[Output<Keccak256>] {
        &self.intermediate_hashes
    }

    pub fn root(&self) -> &Output<Keccak256> {
        &self.root
    }
}

impl<const N: usize> AsRef<[Output<Keccak256>]> for MultilinearBrakedownCommitment<N> {
    fn as_ref(&self) -> &[Output<Keccak256>] {
        slice::from_ref(&self.root)
    }
}

fn err_too_many_variates(function: &str, upto: usize, got: usize) -> Error {
    Error::InvalidPcsParam(
        format!(
            "Too many variates of poly to {function} (param supports variates up to {upto} but got {got})"
        )
    )
}

// Ensures that polynomials and evaluation points are of appropriate size
fn validate_input<'a, const N: usize>(
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

impl<const N: usize, S> MultilinearBrakedown<N, S>
where
    S: BrakedownSpec,
{
    pub type Param = MultilinearBrakedownParams<N>;
    pub type ProverParam = MultilinearBrakedownParams<N>;
    pub type VerifierParam = MultilinearBrakedownParams<N>;
    pub type Polynomial = DenseMultilinearExtension<N>;
    pub type Commitment = MultilinearBrakedownCommitment<N>;
    pub type CommitmentChunk = Output<Keccak256>;

    pub fn setup(poly_size: usize, rng: impl RngCore) -> Self::Param {
        assert!(poly_size.is_power_of_two());
        let num_vars = poly_size.ilog2() as usize;
        let brakedown = Brakedown::new_multilinear::<S>(num_vars, 20.min((1 << num_vars) - 1), rng);
        MultilinearBrakedownParams {
            num_vars,
            num_rows: (1 << num_vars) / brakedown.row_len(),
            brakedown,
        }
    }

    pub fn trim(
        param: &Self::Param,
        poly_size: usize,
        _: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), Error> {
        assert!(poly_size.is_power_of_two());
        if poly_size == 1 << param.num_vars {
            Ok((param.clone(), param.clone()))
        } else {
            Err(Error::InvalidPcsParam(
                "Can't trim MultilinearBrakedownParams into different poly_size".to_string(),
            ))
        }
    }

    pub fn commit(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, Error> {
        validate_input("commit", pp.num_vars(), [poly], None)?;

        let row_len = pp.brakedown.row_len();
        let codeword_len = pp.brakedown.codeword_len();
        let merkle_depth = codeword_len.next_power_of_two().ilog2() as usize;

        let mut rows = vec![F::zero(); pp.num_rows * codeword_len];
        let mut hashes = vec![Output::<Keccak256>::default(); (2 << merkle_depth) - 1];

        Self::encode_rows(pp, codeword_len, row_len, &mut rows, poly);
        Self::compute_column_hashes(&mut hashes, codeword_len, &rows);
        Self::merklize_column_hashes(merkle_depth, &mut hashes);

        let (intermediate_hashes, root) = {
            let mut intermediate_hashes = hashes;
            let root = intermediate_hashes.pop().unwrap();
            (intermediate_hashes, root)
        };

        Ok(MultilinearBrakedownCommitment {
            rows,
            intermediate_hashes,
            root,
        })
    }

    pub fn batch_commit<'a>(
        pp: &Self::ProverParam,
        polys: impl IntoIterator<Item = &'a DenseMultilinearExtension<N>>,
    ) -> Result<Vec<Self::Commitment>, Error> {
        let polys_vec: Vec<&Self::Polynomial> = polys.into_iter().collect();
        polys_vec
            .iter()
            .map(|poly| Self::commit(pp, poly))
            .collect()
    }

    pub fn open(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        comm: &Self::Commitment,
        point: &Vec<F<N>>,
        eval: &F<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<Output<Keccak256>>, Error> {
        validate_input("open", pp.num_vars(), [poly], [point])?;

        let row_len = pp.brakedown.row_len();

        let codeword_len = pp.brakedown.codeword_len();

        // prove proximity
        let (t_0, _t_1) = point_to_tensor(pp.num_rows, point, eval.config_ptr()).unwrap();
        let t_0_combined_row = if pp.num_rows > 1 {
            let combine = |combined_row: &mut [F<N>], coeffs: &[F<N>]| {
                parallelize(combined_row, |(combined_row, offset)| {
                    combined_row
                        .iter_mut()
                        .zip(offset..)
                        .for_each(|(combined, column)| {
                            *combined = F::zero();
                            coeffs
                                .iter()
                                .zip(poly.evaluations.iter().skip(column).step_by(row_len))
                                .for_each(|(coeff, eval)| {
                                    *combined += &(*coeff * eval);
                                });
                        })
                });
            };
            let mut combined_row = vec![F::zero(); row_len];

            for _ in 0..pp.brakedown.num_proximity_testing() {
                let coeffs = transcript
                    .fs_transcript
                    .get_challenges(eval.config_ptr(), pp.num_rows);
                combine(&mut combined_row, &coeffs);
                transcript.write_field_elements(&combined_row)?;
            }
            combine(&mut combined_row, &t_0);
            Cow::<Vec<F<N>>>::Owned(combined_row)
        } else {
            Cow::Borrowed(&poly.evaluations)
        };
        transcript.write_field_elements(&t_0_combined_row)?;

        // open merkle tree
        let depth = codeword_len.next_power_of_two().ilog2() as usize;
        let mut proof: Vec<Output<Keccak256>> = vec![];
        for _ in 0..pp.brakedown.num_column_opening() {
            let column = squeeze_challenge_idx(transcript, eval.config_ptr(), codeword_len);

            transcript.write_field_elements(
                &comm
                    .rows
                    .iter()
                    .copied()
                    .skip(column)
                    .step_by(codeword_len)
                    .collect::<Vec<_>>(),
            )?;

            let mut offset = 0;
            for (idx, width) in (1..=depth).rev().map(|depth| 1 << depth).enumerate() {
                let neighbor_idx = (column >> idx) ^ 1;
                transcript.write_commitment(&comm.intermediate_hashes[offset + neighbor_idx])?;
                proof.push(comm.intermediate_hashes[offset + neighbor_idx]);
                offset += width;
            }
        }

        Ok(proof)
    }

    // TODO: Apply 2022/1355
    pub fn batch_open<'a>(
        pp: &Self::ProverParam,
        polys: impl IntoIterator<Item = &'a DenseMultilinearExtension<N>>,
        comms: impl IntoIterator<Item = &'a MultilinearBrakedownCommitment<N>>,
        points: &[Vec<F<N>>],
        evals: &[F<N>],
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<Vec<Output<Keccak256>>>, Error> {
        //	use std::env;
        //	let key = "RAYON_NUM_THREADS";
        //	env::set_var(key, "8");
        let polys = polys.into_iter().collect_vec();
        let comms = comms.into_iter().collect_vec();
        let mut proofs = vec![];
        for (i, eval) in evals.iter().enumerate() {
            proofs.push(Self::open(
                pp, polys[i], comms[i], &points[i], eval, transcript,
            )?);
        }
        Ok(proofs)
    }

    pub fn read_commitments(
        _: &Self::VerifierParam,
        num_polys: usize,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<Vec<Self::Commitment>, Error> {
        transcript.read_commitments(num_polys).map(|roots| {
            roots
                .into_iter()
                .map(MultilinearBrakedownCommitment::from_root)
                .collect_vec()
        })
    }

    pub fn verify(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &Vec<F<N>>,
        eval: &F<N>,
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        validate_input("verify", vp.num_vars(), [], [point])?;

        let row_len = vp.brakedown.row_len();
        let codeword_len = vp.brakedown.codeword_len();

        let (t_0, t_1) = point_to_tensor(vp.num_rows, point, eval.config_ptr())?;
        let mut combined_rows = Vec::with_capacity(vp.brakedown.num_proximity_testing() + 1);
        if vp.num_rows > 1 {
            let coeffs = transcript
                .fs_transcript
                .get_challenges(eval.config_ptr(), vp.num_rows);
            let mut combined_row = transcript.read_field_elements(row_len, eval.config_ptr())?;
            combined_row.resize(codeword_len, F::zero());
            vp.brakedown.encode(&mut combined_row);
            combined_rows.push((coeffs, combined_row));
        }
        combined_rows.push({
            let mut combined_row = transcript.read_field_elements(row_len, eval.config_ptr())?;
            combined_row.resize(codeword_len, F::zero());
            vp.brakedown.encode(&mut combined_row);
            (t_0, combined_row)
        });

        let depth = codeword_len.next_power_of_two().ilog2() as usize;

        for _ in 0..vp.brakedown.num_column_opening() {
            let column = squeeze_challenge_idx(transcript, eval.config_ptr(), codeword_len);
            let items = transcript.read_field_elements(vp.num_rows, eval.config_ptr())?;
            let path = transcript.read_commitments(depth)?;

            // verify proximity
            for (coeff, encoded) in combined_rows.iter() {
                let item = if vp.num_rows > 1 {
                    inner_product(coeff, &items)
                } else {
                    items[0]
                };
                if item != encoded[column] {
                    return Err(Error::InvalidPcsOpen("Proximity failure".to_string()));
                }
            }

            // verify merkle tree opening
            let mut hasher = Keccak256::default();
            let mut output = {
                for item in items.iter() {
                    <Keccak256 as sha3::digest::Update>::update(
                        &mut hasher,
                        &item.value().to_bytes_be(),
                    );
                }

                hasher.clone().finalize()
            };
            for (idx, neighbor) in path.iter().enumerate() {
                if (column >> idx) & 1 == 0 {
                    <Keccak256 as sha3::digest::Update>::update(&mut hasher, &output);
                    <Keccak256 as sha3::digest::Update>::update(&mut hasher, neighbor);
                } else {
                    <Keccak256 as sha3::digest::Update>::update(&mut hasher, neighbor);
                    <Keccak256 as sha3::digest::Update>::update(&mut hasher, &output);
                }
                output = hasher.clone().finalize();
            }
            if &output != comm.root() {
                return Err(Error::InvalidPcsOpen(
                    "Invalid merkle tree opening".to_string(),
                ));
            }
        }

        // verify consistency
        let t_0_combined_row = combined_rows
            .last()
            .map(|(_, combined_row)| &combined_row[..row_len])
            .unwrap();
        if inner_product(t_0_combined_row, &t_1) != *eval {
            return Err(Error::InvalidPcsOpen("Consistency failure".to_string()));
        }

        Ok(())
    }

    pub fn batch_verify<'a>(
        vp: &Self::VerifierParam,
        comms: impl IntoIterator<Item = &'a MultilinearBrakedownCommitment<N>>,
        points: &[Vec<F<N>>],
        evals: &[F<N>],
        transcript: &mut PcsTranscript<N>,
    ) -> Result<(), Error> {
        let comms = comms.into_iter().collect_vec();
        for (i, eval) in evals.iter().enumerate() {
            Self::verify(vp, comms[i], &points[i], eval, transcript)?;
        }
        Ok(())
    }

    fn encode_rows(
        pp: &Self::ProverParam,
        codeword_len: usize,
        row_len: usize,
        rows: &mut [F<N>],
        poly: &Self::Polynomial,
    ) {
        let chunk_size = div_ceil(pp.num_rows, num_threads());
        parallelize_iter(
            rows.chunks_exact_mut(chunk_size * codeword_len)
                .zip(poly.evaluations.chunks_exact(chunk_size * row_len)),
            |(rows, evals)| {
                for (row, evals) in rows
                    .chunks_exact_mut(codeword_len)
                    .zip(evals.chunks_exact(row_len))
                {
                    row[..evals.len()].copy_from_slice(evals);
                    pp.brakedown.encode(row);
                }
            },
        );
    }

    fn compute_column_hashes(hashes: &mut [Output<Keccak256>], codeword_len: usize, rows: &[F<N>]) {
        parallelize(&mut hashes[..codeword_len], |(hashes, start)| {
            let mut hasher = Keccak256::new();
            for (hash, column) in hashes.iter_mut().zip(start..) {
                rows.iter()
                    .skip(column)
                    .step_by(codeword_len)
                    .for_each(|item| {
                        <Keccak256 as sha3::digest::Update>::update(
                            &mut hasher,
                            &item.value().to_bytes_be(),
                        )
                    });
                hasher.finalize_into_reset(hash);
            }
        });
    }

    fn merklize_column_hashes(depth: usize, hashes: &mut [Output<Keccak256>]) {
        let mut offset = 0;
        for width in (1..=depth).rev().map(|depth| 1 << depth) {
            let (input, output) = hashes[offset..].split_at_mut(width);
            //	    let num_threads = env::var("RAYON_NUM_THREADS").unwrap();
            let chunk_size = div_ceil(output.len(), num_threads());
            parallelize_iter(
                input
                    .chunks(2 * chunk_size)
                    .zip(output.chunks_mut(chunk_size)),
                |(input, output)| {
                    let mut hasher = Keccak256::new();

                    for (input, output) in input.chunks_exact(2).zip(output.iter_mut()) {
                        <Keccak256 as sha3::digest::Update>::update(&mut hasher, &input[0]);
                        <Keccak256 as sha3::digest::Update>::update(&mut hasher, &input[1]);
                        hasher.finalize_into_reset(output);
                    }
                },
            );
            offset += width;
        }
    }
}

fn point_to_tensor<const N: usize>(
    num_rows: usize,
    point: &[F<N>],
    config: *const FieldConfig<N>,
) -> Result<(Vec<F<N>>, Vec<F<N>>), Error> {
    assert!(num_rows.is_power_of_two());
    let (hi, lo) = point.split_at(point.len() - num_rows.ilog2() as usize);
    // TODO: get rid of these unwraps.
    let t_0 = build_eq_x_r(lo, config).unwrap();
    let t_1 = build_eq_x_r(hi, config).unwrap();
    Ok((t_0.evaluations, t_1.evaluations))
}

fn squeeze_challenge_idx<const N: usize>(
    transcript: &mut PcsTranscript<N>,
    config: *const FieldConfig<N>,
    cap: usize,
) -> usize {
    let challenge = transcript.fs_transcript.get_challenge(config);
    let mut bytes = [0; size_of::<u32>()];
    bytes.copy_from_slice(&challenge.value().to_bytes_be()[..size_of::<u32>()]);
    u32::from_le_bytes(bytes) as usize % cap
}
