use crate::rand::PRNG;

use crate::linalg::cross;
use crate::linalg::Vec3D;

use crate::ray::Ray;

#[derive(Copy, Clone)]
pub struct Camera {
    pub origin: Vec3D,
    pub norm: Vec3D,
    pub up: Vec3D,
    pub height: f32,
    pub width: f32,
    pub aperture: f32,
    pub focal_distance: f32,
}

impl Camera {
    pub fn prime_ray(self, rng: &mut PRNG, row_frac: f32, col_frac: f32) -> Ray {
        let height = self.height as f32;
        let width = self.width as f32;
        let focal_distance = self.focal_distance;
        let lens_radius = self.aperture / 2.0;
        let origin = self.origin;
        let basis1 = cross(&self.up, &self.norm).l2_normalize();
        let basis2 = cross(&self.norm, &basis1);
        let left_top = origin - (focal_distance * width * basis1 / 2.0)
            + (focal_distance * height * basis2 / 2.0)
            - focal_distance * self.norm;
        let mut random_lens_offset = lens_radius * Vec3D::random_unit_disk_vector(rng);
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

pub fn fov_to_imsize(fov_radians: f32, aspect_ratio: f32) -> (f32, f32) {
    let h = (fov_radians / 2.0).tan();
    let image_plane_height = 2.0 * h;
    let image_plane_width = aspect_ratio * image_plane_height;
    (image_plane_height, image_plane_width)
}

pub fn imsize_to_fov(image_plane_height: f32, image_plane_width: f32) -> (f32, f32) {
    let aspect_ratio = image_plane_width / image_plane_height;
    let h = image_plane_height / 2.0;
    let fov = h.atan() / 2.0;
    (fov, aspect_ratio)
}

pub fn imsize_to_pixels(
    image_plane_height: f32,
    image_plane_width: f32,
    pixel_size: f32,
) -> (u32, u32, f32) {
    let rows = (image_plane_height / pixel_size).round() as u32;
    let cols = (image_plane_width / pixel_size).round() as u32;
    let pixel_size = image_plane_height / (rows as f32);
    (rows, cols, pixel_size)
}
