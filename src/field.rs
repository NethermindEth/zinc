use std::ops::Add;

use ark_ff::BigInt;

pub struct RandomFieldConfig {}

pub struct RandomField<'config> {
    pub config: &'config RandomFieldConfig,
    pub value: BigInt<8>,
}

impl<'config> RandomField<'config> {
    pub fn new(config: &'config RandomFieldConfig, value: BigInt<8>) -> Self {
        RandomField { config, value }
    }
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
        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.

        let config_ptr_lhs: *const RandomFieldConfig = self.config;
        let config_ptr_rhs: *const RandomFieldConfig = rhs.config;

        if config_ptr_lhs != config_ptr_rhs {
            panic!("cannot add field elements of different fields");
        }

        todo!()
    }
}
