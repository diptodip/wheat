use std::f64::consts::PI;

use crate::rand::prelude::*;

use crate::linalg::dot;
use crate::linalg::Vec3D;

use crate::colors::rgb;
use crate::colors::vec_to_rgb;
use crate::colors::RGB;

use crate::geometry::first_intersection;
use crate::geometry::Intersectable;
use crate::geometry::Intersection;
use crate::geometry::Intersects;

use crate::materials::Material;
use crate::materials::Surface;

pub struct Ray {
    pub origin: Vec3D,
    pub direction: Vec3D,
}

impl Ray {
    pub fn at(&self, t: f64) -> Vec3D {
        self.origin + t * self.direction
    }
}

pub fn diffuse_bounce(intersection: &Intersection, intersectable: &Intersectable) -> Ray {
    let bounce_vector = intersection.local_normal + Vec3D::random_unit_vector();
    Ray {
        origin: intersection.point,
        direction: bounce_vector,
    }
}

pub fn reflect(intersection: &Intersection, intersectable: &Intersectable, ray: &Ray) -> Ray {
    let direction = ray.direction;
    let normal = intersection.local_normal;
    let reflected = direction - 2.0 * dot(&direction, &normal) * normal;
    Ray {
        origin: intersection.point,
        direction: reflected,
    }
}

pub fn fuzzy_reflect(intersection: &Intersection, intersectable: &Intersectable, ray: &Ray) -> Ray {
    let direction = ray.direction;
    let normal = intersection.local_normal;
    let reflected = direction - 2.0 * dot(&direction, &normal) * normal;
    let mut fuzz: f64 = 0.0;
    let material = intersectable.material();
    let surface = material.surface;
    if let Surface::FuzzyReflective(surface_fuzz) = surface {
        fuzz = surface_fuzz;
    }
    let direction_fuzz = fuzz * Vec3D::random_unit_vector();
    Ray {
        origin: intersection.point,
        direction: reflected + direction_fuzz,
    }
}

#[inline]
fn schlick(cos_theta: f64, index_ratio: f64) -> f64 {
    let mut r0 = (1.0 - index_ratio) / (1.0 + index_ratio);
    r0 = r0 * r0;
    r0 + (1.0 - r0) * (1.0 - cos_theta).powf(5.0)
}

pub fn refract(intersection: &Intersection, intersectable: &Intersectable, ray: &Ray) -> Ray {
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
    let mut rng = rand::thread_rng();
    let reflect_check: f64 = rng.gen();
    if reflect_check < reflectivity {
        return reflect(intersection, intersectable, ray);
    }
    let r2_par = index_ratio * (r1 + cos_theta1 * intersection.local_normal);
    let r2_per = (-(1.0 - r2_par.length_squared()).sqrt() * intersection.local_normal);
    Ray {
        origin: intersection.point,
        direction: r2_par + r2_per,
    }
}

pub fn trace(ray: &Ray, world: &Vec<Intersectable>, depth: u64) -> RGB {
    // light enters the void if we hit the depth limit
    if depth <= 0 {
        return rgb(0.0, 0.0, 0.0);
    }
    // calculate intersection list
    let mut intersections = Vec::new();
    for intersectable in world {
        intersections.push(intersectable.intersects(ray));
    }
    // determine if ray intersects and choose first intersection if so
    let result = first_intersection(intersections, world);
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
                let bounced = &diffuse_bounce(&intersection, intersectable);
                traced_color = trace(bounced, world, depth - 1).to_vec3d();
            } else if let Surface::Reflective = surface {
                // light is reflected if material is totally reflective,
                // so we recurse and trace a reflected ray
                let reflected = &reflect(&intersection, intersectable, ray);
                traced_color = trace(reflected, world, depth - 1).to_vec3d();
            } else if let Surface::FuzzyReflective(fuzz) = surface {
                // light is reflected with some random offset
                // if material is fuzzy reflective,
                // so we recurse and trace a fuzzy reflected ray
                // with a check to make sure the reflection is correct
                let reflected = &fuzzy_reflect(&intersection, intersectable, ray);
                if dot(&reflected.direction, &intersection.local_normal) > 0.0 {
                    traced_color = trace(reflected, world, depth - 1).to_vec3d();
                }
            } else if let Surface::Refractive(r) = surface {
                // light is refracted or reflected depending on angle,
                // so we recurse to trace either a refracted/reflected ray
                let refracted = &refract(&intersection, intersectable, ray);
                traced_color = trace(refracted, world, depth - 1).to_vec3d();
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
                (1.0 - height) + height * 1.0,
            );
        }
    }
}
