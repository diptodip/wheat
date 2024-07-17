use crate::linalg::Vec3D;
use crate::rand::prelude::*;

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

    pub fn random() -> RGB {
        let mut rng = thread_rng();
        let color_vec = Vec3D::random(&mut rng, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        vec_to_rgb(color_vec)
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

pub fn vec_to_rgb(color_vec: Vec3D) -> RGB {
    let color_vec = color_vec.clamp(0.0, 1.0);
    RGB {
        r: color_vec.0,
        g: color_vec.1,
        b: color_vec.2,
    }
}
