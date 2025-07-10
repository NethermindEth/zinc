pub trait Field {
    type I: Integer;
    type C: Config<Self::I>;
    type Cr: ConfigReference<Self::I, Self::C>;
}

pub trait Integer {}
pub trait Config<I: Integer> {}
pub trait ConfigReference<I: Integer, C: Config<I>> {}
