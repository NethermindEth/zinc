use ark_std::ops::Mul;
use num_traits::{One, Zero};

pub trait Field:
    Zero + One + Clone + Copy + Default + Sync + Send + for<'a> Mul<&'a Self, Output = Self>
{
    type I: Integer;
    type C: Config<Self::I>;
    type Cr: ConfigReference<Self::I, Self::C>;
}

pub trait Integer {}
pub trait Config<I: Integer> {}
pub trait ConfigReference<I: Integer, C: Config<I>> {}
