use crate::linalg::cross;
use crate::linalg::Vec3D;

use crate::ray::Ray;

use crate::colors::RGB;

#[derive(Copy, Clone)]
pub struct Camera {
    pub origin: Vec3D,
    pub norm: Vec3D,
    pub up: Vec3D,
    pub height: f64,
    pub width: f64,
    pub aperture: f64,
    pub focal_distance: f64,
}

impl Camera {
    pub fn prime_ray(self, row_frac: f64, col_frac: f64) -> Ray {
        let height = self.height as f64;
        let width = self.width as f64;
        let focal_distance = self.focal_distance;
        let lens_radius = self.aperture / 2.0;
        let origin = self.origin;
        let basis1 = cross(&self.up, &self.norm).l2_normalize();
        let basis2 = cross(&self.norm, &basis1);
        let left_top = origin - (focal_distance * width * basis1 / 2.0)
            + (focal_distance * height * basis2 / 2.0)
            - focal_distance * self.norm;
        let mut random_lens_offset = lens_radius * Vec3D::random_unit_disk_vector();
        random_lens_offset = basis1 * random_lens_offset.0 + basis2 * random_lens_offset.1;
        Ray {
            origin: origin + random_lens_offset,
            direction: left_top + (focal_distance * basis1 * col_frac * width)
                - (focal_distance * basis2 * row_frac * height)
                - origin
                - random_lens_offset,
        }
    }
}

pub fn fov_to_imsize(fov_radians: f64, aspect_ratio: f64) -> (f64, f64) {
    let h = (fov_radians / 2.0).tan();
    let image_plane_height = 2.0 * h;
    let image_plane_width = aspect_ratio * image_plane_height;
    (image_plane_height, image_plane_width)
}

pub fn imsize_to_fov(image_plane_height: f64, image_plane_width: f64) -> (f64, f64) {
    let aspect_ratio = image_plane_width / image_plane_height;
    let h = image_plane_height / 2.0;
    let fov = h.atan() / 2.0;
    (fov, aspect_ratio)
}

pub fn imsize_to_pixels(
    image_plane_height: f64,
    image_plane_width: f64,
    pixel_size: f64,
) -> (u32, u32, f64) {
    let rows = (image_plane_height / pixel_size).round() as u32;
    let cols = (image_plane_width / pixel_size).round() as u32;
    let pixel_size = image_plane_height / (rows as f64);
    (rows, cols, pixel_size)
}
