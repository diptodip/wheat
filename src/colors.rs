use crate::linalg::Vec3D;

#[derive(Copy, Clone)]
pub struct RGB {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl RGB {
    pub fn to_vec3d(self) -> Vec3D {
        Vec3D(self.r, self.g, self.b)
    }

    pub fn random() -> RGB {
        let color_vec = Vec3D::random(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        vec_to_rgb(color_vec)
    }
}

pub fn rgb(r: f64, g: f64, b: f64) -> RGB {
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
