use std::ops::Add;

use ark_ff::BigInt;

pub struct FieldConfig {}

pub struct RandomField<'config, const N: usize> {
    pub config: &'config FieldConfig,
    pub value: BigInt<N>,
}

impl<'config, const N: usize> RandomField<'config, N> {
    pub fn new(config: &'config FieldConfig, value: BigInt<N>) -> Self {
        RandomField { config, value }
    }
}

impl<'a, 'config, const N: usize> Add<&'a RandomField<'config, N>> for RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        &self + rhs
    }
}

impl<'a, 'config, const N: usize> Add<&'a RandomField<'config, N>> for &RandomField<'config, N> {
    type Output = RandomField<'config, N>;

    fn add(self, rhs: &'a RandomField<'config, N>) -> RandomField<'config, N> {
        // Here we assume that the elements of a random field are
        // created using the same RandomFieldConfig.

        let config_ptr_lhs: *const FieldConfig = self.config;
        let config_ptr_rhs: *const FieldConfig = rhs.config;

        if config_ptr_lhs != config_ptr_rhs {
            panic!("cannot add field elements of different fields");
        }

        todo!()
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_add() {
        // TODO: fill this in
    }
}
