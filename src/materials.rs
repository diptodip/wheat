use crate::colors::RGB;

#[derive(Copy, Clone)]
pub enum Surface {
    Diffuse,
    Reflective,
    FuzzyReflective(f32),
    Refractive(f32),
}

#[derive(Copy, Clone)]
pub struct Material {
    pub color: RGB,
    pub surface: Surface,
}
