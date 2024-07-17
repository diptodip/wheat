use std::f32::consts::PI;
use std::time::Instant;

use crate::rand::prelude::*;

use crate::linalg::dot;
use crate::linalg::Vec3D;

use crate::colors::rgb;
use crate::colors::vec_to_rgb;
use crate::colors::RGB;

use crate::io::write_ppm;

use crate::geometry::find_intersections;
use crate::geometry::Intersectable;
use crate::geometry::Intersection;
use crate::geometry::Intersects;

use crate::camera::Camera;

use crate::materials::Material;
use crate::materials::Surface;

use rayon::current_num_threads;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Ray {
    pub origin: Vec3D,
    pub direction: Vec3D,
}

impl Ray {
    pub fn at(&self, t: f32) -> Vec3D {
        self.origin + t * self.direction
    }
}

fn diffuse_bounce(
    rng: &mut ThreadRng,
    intersection: &Intersection,
    intersectable: &Intersectable,
) -> Ray {
    let bounce_vector = intersection.local_normal + Vec3D::random_unit_vector(rng);
    Ray {
        origin: intersection.point,
        direction: bounce_vector,
    }
}

fn reflect(intersection: &Intersection, intersectable: &Intersectable, ray: &Ray) -> Ray {
    let direction = ray.direction;
    let normal = intersection.local_normal;
    let reflected = direction - 2.0 * dot(&direction, &normal) * normal;
    Ray {
        origin: intersection.point,
        direction: reflected,
    }
}

fn fuzzy_reflect(
    rng: &mut ThreadRng,
    intersection: &Intersection,
    intersectable: &Intersectable,
    ray: &Ray,
) -> Ray {
    let direction = ray.direction;
    let normal = intersection.local_normal;
    let reflected = direction - 2.0 * dot(&direction, &normal) * normal;
    let mut fuzz: f32 = 0.0;
    let material = intersectable.material();
    let surface = material.surface;
    if let Surface::FuzzyReflective(surface_fuzz) = surface {
        fuzz = surface_fuzz;
    }
    let direction_fuzz = fuzz * Vec3D::random_unit_vector(rng);
    Ray {
        origin: intersection.point,
        direction: reflected + direction_fuzz,
    }
}

#[inline]
fn schlick(cos_theta: f32, index_ratio: f32) -> f32 {
    let mut r0 = (1.0 - index_ratio) / (1.0 + index_ratio);
    r0 = r0 * r0;
    r0 + (1.0 - r0) * (1.0 - cos_theta).powf(5.0)
}

fn refract(
    rng: &mut ThreadRng,
    intersection: &Intersection,
    intersectable: &Intersectable,
    ray: &Ray,
) -> Ray {
    let mut material_index = 1.0;
    let material = intersectable.material();
    let surface = material.surface;
    if let Surface::Refractive(index) = surface {
        material_index = index;
    }
    let index_ratio = if intersection.inside {
        material_index
    } else {
        1.0 / material_index
    };
    let r1 = ray.direction.l2_normalize();
    let cos_theta1 = dot(&(-r1), &intersection.local_normal).min(1.0);
    let sin_theta1 = (1.0 - cos_theta1 * cos_theta1).sqrt();
    if (index_ratio * sin_theta1) > 1.0 {
        return reflect(intersection, intersectable, ray);
    }
    let reflectivity = schlick(cos_theta1, index_ratio);
    // let mut rng = rand::thread_rng();
    let reflect_check: f32 = rng.gen();
    if reflect_check < reflectivity {
        return reflect(intersection, intersectable, ray);
    }
    let r2_per = index_ratio * (r1 + cos_theta1 * intersection.local_normal);
    let r2_par = (-(1.0 - r2_per.length_squared()).sqrt() * intersection.local_normal);
    Ray {
        origin: intersection.point,
        direction: r2_per + r2_par,
    }
}

fn trace(
    rng: &mut ThreadRng,
    ray: &Ray,
    world: &Vec<Intersectable>,
    depth: u64,
    num_casted: &mut u64,
) -> RGB {
    // light enters the void if we hit the depth limit
    if depth <= 0 {
        return rgb(0.0, 0.0, 0.0);
    }
    *num_casted += 1;
    // determine if ray intersects and choose first intersection if so
    // let (intersections, result) = find_intersections(ray, world);
    let result = find_intersections(ray, world);
    match result {
        // calculate color at intersection point
        Some((intersection, intersectable)) => {
            let material = intersectable.material();
            let surface = material.surface;
            let material_color = material.color.to_vec3d();
            let mut traced_color = Vec3D(0.0, 0.0, 0.0);
            if let Surface::Diffuse = surface {
                // light bounces randomly if material is diffuse,
                // so we recurse and trace a randomly bounced ray
                let bounced = &diffuse_bounce(rng, &intersection, intersectable);
                traced_color = trace(rng, bounced, world, depth - 1, num_casted).to_vec3d();
            } else if let Surface::Reflective = surface {
                // light is reflected if material is totally reflective,
                // so we recurse and trace a reflected ray
                let reflected = &reflect(&intersection, intersectable, ray);
                traced_color = trace(rng, reflected, world, depth - 1, num_casted).to_vec3d();
            } else if let Surface::FuzzyReflective(_fuzz) = surface {
                // light is reflected with some random offset
                // if material is fuzzy reflective,
                // so we recurse and trace a fuzzy reflected ray
                // with a check to make sure the reflection is correct
                let reflected = &fuzzy_reflect(rng, &intersection, intersectable, ray);
                if dot(&reflected.direction, &intersection.local_normal) > 0.0 {
                    traced_color = trace(rng, reflected, world, depth - 1, num_casted).to_vec3d();
                }
            } else if let Surface::Refractive(_r) = surface {
                // light is refracted or reflected depending on angle,
                // so we recurse to trace either a refracted/reflected ray
                let refracted = &refract(rng, &intersection, intersectable, ray);
                traced_color = trace(rng, refracted, world, depth - 1, num_casted).to_vec3d();
            }
            let color_vec = material_color * traced_color;
            return vec_to_rgb(color_vec);
        }
        None => {
            let ray_direction = ray.direction.l2_normalize();
            let height = 0.5 * (ray_direction.1 + 1.0);
            return rgb(
                (1.0 - height) + height * 0.5,
                (1.0 - height) + height * 0.7,
                1.0,
            );
        }
    }
}

pub fn render(
    world: &Vec<Intersectable>,
    camera: &Camera,
    rows: usize,
    cols: usize,
    samples_per_pixel: f32,
) -> u64 {
    // construct blank image
    let mut image = vec![vec![0.0; 3]; rows * cols];
    let num_threads = current_num_threads();
    eprintln!(
        "[start] rendering {}px x {}px (width x height) on {} threads",
        cols, rows, num_threads
    );
    eprintln!("[info] {:.2}%", 0.0);
    // need an atomic counter so rust doesn't complain about thread safety
    let completed_rows = AtomicUsize::new(0);
    let num_rays: u64 = 0;
    // iterate through pixels in parallel
    let start = Instant::now();
    image
        .par_chunks_mut(cols)
        .enumerate()
        .for_each(|(row, chunk)| {
            // create RNG
            let mut rng = rand::thread_rng();
            let mut num_traced = 0;
            // loop through columns in row
            for col in 0..cols {
                if (col == (cols - 1)) {
                    // have to use atomic counter updating functions here
                    completed_rows.store(
                        completed_rows.load(Ordering::Relaxed) + 1 as usize,
                        Ordering::Relaxed,
                    );
                    eprintln!(
                        "[info] {:.2}%",
                        completed_rows.load(Ordering::Relaxed) as f32 / rows as f32 * 100.0
                    );
                }
                let mut r: f32 = 0.0;
                let mut g: f32 = 0.0;
                let mut b: f32 = 0.0;
                for _sample in 0..samples_per_pixel as usize {
                    // calculate ray for current pixel
                    // making sure to center ray within pixel
                    // we also perturb the ray direction slightly per sample
                    let row_rand = rng.gen::<f32>();
                    let col_rand = rng.gen::<f32>();
                    let row_frac = (row as f32 + 0.5 + row_rand) / (rows as f32);
                    let col_frac = (col as f32 + 0.5 + col_rand) / (cols as f32);
                    let ray = camera.prime_ray(&mut rng, row_frac, col_frac);
                    // trace ray for current pixel
                    let color = trace(&mut rng, &ray, world, 50, &mut num_traced);
                    r += color.r;
                    g += color.g;
                    b += color.b;
                }
                r /= samples_per_pixel;
                g /= samples_per_pixel;
                b /= samples_per_pixel;
                let pixel = &mut chunk[col];
                pixel[0] = r;
                pixel[1] = g;
                pixel[2] = b;
            }
        });
    let duration = start.elapsed();
    println!("[info] scene rendered in {:.2?}", duration);
    eprintln!("[info] writing image...");
    _ = write_ppm(image, rows, cols);
    eprintln!("[ok] done!");
    num_rays
}
