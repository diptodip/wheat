use crate::linalg::Vec3D;
use crate::rand::PRNG;
// use crate::rand::prelude::*;

#[derive(Copy, Clone)]
pub struct RGB {
    pub r: f32,
    pub g: f32,
    pub b: f32,
}

impl RGB {
    pub fn to_vec3d(self) -> Vec3D {
        Vec3D(self.r, self.g, self.b)
    }

    pub fn random(rng: &mut PRNG) -> RGB {
        let color_vec = Vec3D::random(rng, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        RGB::from_vec(color_vec)
    }

    pub fn from_vec(color_vec: Vec3D) -> RGB {
        let color_vec = color_vec.clamp(0.0, 1.0);
        RGB {
            r: color_vec.0,
            g: color_vec.1,
            b: color_vec.2,
        }
    }
}

pub fn rgb(r: f32, g: f32, b: f32) -> RGB {
    let color_vec = Vec3D(r, g, b).clamp(0.0, 1.0);
    RGB {
        r: color_vec.0,
        g: color_vec.1,
        b: color_vec.2,
    }
}
