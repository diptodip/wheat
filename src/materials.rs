use crate::colors::RGB;

#[derive(Copy, Clone)]
pub enum Surface {
    Diffuse,
    Reflective,
    FuzzyReflective(f64),
    Refractive(f64),
}

#[derive(Copy, Clone)]
pub struct Material {
    pub color: RGB,
    pub surface: Surface,
}
