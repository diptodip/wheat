use std::option::Option;
use crate::linalg::Vec3D;
use crate::linalg::dot;

pub struct Ray {
    pub origin: Vec3D,
    pub direction: Vec3D,
}

impl Ray {
    pub fn at(&self, t:f64) -> Vec3D {
        t * self.direction
    }
}

pub struct Sphere {
    pub origin: Vec3D,
    pub radius: f64,
}

impl Sphere {
    pub fn intersects(&self, ray: &Ray) -> Option<Vec3D> {
        // using quadratic formula
        let sphere_to_ray = ray.origin - self.origin;
        let a = dot(&ray.direction, &ray.direction);
        let h = dot(&sphere_to_ray, &ray.direction);
        let c = dot(&sphere_to_ray, &sphere_to_ray) - self.radius * self.radius;
        let discriminant = (h * h) - (a * c);
        if discriminant > 0.0 {
            let t = (-h - discriminant.sqrt()) / a;
            if t > 0.0 {
                return Some(ray.at(t));
            }
        }
        return None;
    }
}
