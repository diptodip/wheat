use crate::linalg::Vec3D;

#[derive(Copy, Clone)]
pub struct RGB {
    r: f64,
    g: f64,
    b: f64,
}

impl RGB {
    pub fn r(self) -> f64 {
        return self.r;
    }

    pub fn g(self) -> f64 {
        return self.g;
    }

    pub fn b(self) -> f64 {
        return self.b;
    }
}

pub fn rgb(r: f64, g: f64, b: f64) -> RGB {
    let colors = Vec3D(r, g, b).max_normalize();
    RGB { r: colors.0, g: colors.1, b: colors.2 }
}
