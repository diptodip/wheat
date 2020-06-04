use crate::linalg::Vec3D;
use crate::linalg::cross;

use crate::geometry::Ray;

#[derive(Copy, Clone)]
pub struct ImagePlane {
    pub origin: Vec3D,
    pub norm: Vec3D,
    pub up: Vec3D,
    pub height: f64,
    pub width: f64,
}

#[derive(Copy, Clone)]
pub struct Camera {
    pub origin: Vec3D,
    pub image_plane: ImagePlane,
}

impl Camera {
    pub fn ray_from_pixel(self, row_frac: f64, col_frac: f64) -> Ray {
        let height = self.image_plane.height as f64;
        let width = self.image_plane.width as f64;
        let basis1 = cross(&self.image_plane.up,
                           &self.image_plane.norm).l2_normalize();
        let basis2 = cross(&self.image_plane.norm, &basis1).l2_normalize();
        let left_top = (self.image_plane.origin
                        - (width * basis1 / 2.0)
                        + (height * basis2 / 2.0));
        return Ray {
            origin: self.origin,
            direction: (left_top
                        + (basis1 * col_frac * width)
                        - (basis2 * row_frac * height)
                        - self.origin),
        };
    }
}

pub fn fov_to_imsize(fov_radians: f64, aspect_ratio: f64) -> (f64, f64) {
    let h = (fov_radians / 2.0).tan();
    let image_plane_height = 2.0 * h;
    let image_plane_width = aspect_ratio * image_plane_height;
    (image_plane_height, image_plane_width)
}

pub fn imsize_to_fov(image_plane_height: f64,
                     image_plane_width: f64) -> (f64, f64) {
    let aspect_ratio = image_plane_width / image_plane_height;
    let h = image_plane_height / 2.0;
    let fov = h.atan() / 2.0;
    (fov, aspect_ratio)
}

pub fn imsize_to_pixels(image_plane_height: f64,
                        image_plane_width: f64,
                        pixel_size: f64) -> (u32, u32, f64) {
    let rows = (image_plane_height / pixel_size) as u32;
    let cols = (image_plane_width / pixel_size) as u32;
    let pixel_size = image_plane_height / (rows as f64);
    (rows, cols, pixel_size)
}
