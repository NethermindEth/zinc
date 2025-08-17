use ark_std::{collections::BTreeSet, ops::Range, vec::Vec};
use crypto_bigint::{Int, Word};

use crate::{
    poly::dense::DenseMultilinearExtension,
    traits::Transcript,
    zip::{
        code::{DefaultLinearCodeSpec, ZipLinearCode},
        pcs::structs::MultilinearZipParams,
    },
};

const INT_LIMBS: usize = 1;

const N: usize = INT_LIMBS;
const L: usize = INT_LIMBS * 2;
const K: usize = INT_LIMBS * 4;
const M: usize = INT_LIMBS * 8;

#[derive(Default)]
pub struct MockTranscript {
    pub counter: i64,
}

impl<const L: usize> Transcript<L> for MockTranscript {
    fn get_encoding_element(&mut self) -> Int<L> {
        self.counter += 1;
        Int::from(self.counter)
    }
    fn get_word(&mut self) -> Word {
        self.counter += 1;
        self.counter as Word
    }
    fn sample_unique_columns(
        &mut self,
        range: Range<usize>,
        columns: &mut BTreeSet<usize>,
        count: usize,
    ) -> usize {
        self.counter += 1;
        let mut inserted = 0;
        for i in range.clone() {
            if columns.insert(i) {
                inserted += 1;
                if inserted == count {
                    break;
                }
            }
        }
        inserted
    }
}

pub fn setup_test_params(
    num_vars: usize,
) -> (
    MultilinearZipParams<N, L, K, M, ZipLinearCode<N, L, K, M>>,
    DenseMultilinearExtension<Int<INT_LIMBS>>,
) {
    let poly_size = 1 << num_vars;
    let num_rows = 1 << num_vars.div_ceil(2);

    let mut transcript = MockTranscript::default();
    let code = ZipLinearCode::<N, L, K, M>::new(&DefaultLinearCodeSpec, poly_size, &mut transcript);
    let pp = MultilinearZipParams::new(num_vars, num_rows, code);

    let evaluations: Vec<_> = (1..=poly_size).map(|v| Int::from(v as i32)).collect();
    let poly = DenseMultilinearExtension::from_evaluations_vec(num_vars, evaluations);

    (pp, poly)
}
