use std::ops::Add;

use ark_ff::BigInt;

pub struct RandomFieldConfig {}

pub struct RandomField<'config> {
    pub config: &'config RandomFieldConfig,
    pub value: BigInt<8>,
}

impl<'a, 'config> Add<&'a RandomField<'config>> for RandomField<'config> {
    type Output = RandomField<'config>;

    fn add(self, rhs: &'a RandomField<'config>) -> RandomField<'config> {
        &self + rhs
    }
}

impl<'a, 'b, 'config> Add<&'a RandomField<'config>> for &'b RandomField<'config> {
    type Output = RandomField<'config>;

    fn add(self, rhs: &'a RandomField<'config>) -> RandomField<'config> {
        todo!()
    }
}
