use crate::colors::RGB;

#[derive(Copy, Clone)]
pub struct Diffuse {
    pub color: RGB,
}

#[derive(Copy, Clone)]
pub struct Reflective {
    pub color: RGB,
}

#[derive(Copy, Clone)]
pub enum Material {
    Diffuse(Diffuse),
    Reflective(Reflective),
}
