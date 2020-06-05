use std::option::Option;
use std::f64::INFINITY;
use std::f64::consts::PI;

use crate::rand::prelude::*;

use crate::linalg::Vec3D;
use crate::linalg::dot;

use crate::materials::Material;

pub struct Ray {
    pub origin: Vec3D,
    pub direction: Vec3D,
}

impl Ray {
    pub fn at(&self, t:f64) -> Vec3D {
        self.origin + t * self.direction
    }
}

#[derive(Copy, Clone)]
pub struct Intersection {
    pub point: Vec3D,
    pub distance: f64,
    pub local_normal: Vec3D,
    pub inside: bool,
}

pub trait Intersects {
    fn intersects(&self, ray: &Ray) -> Option<Intersection>;
    fn surface_normal(&self, point: Vec3D) -> Vec3D;
}

pub struct Sphere {
    pub origin: Vec3D,
    pub radius: f64,
}

impl Intersects for Sphere {
    fn surface_normal(&self, point: Vec3D) -> Vec3D {
        (point - self.origin).l2_normalize()
    }

    fn intersects(&self, ray: &Ray) -> Option<Intersection> {
        // using quadratic formula
        let sphere_to_ray = ray.origin - self.origin;
        let a = dot(&ray.direction, &ray.direction);
        let h = dot(&sphere_to_ray, &ray.direction);
        let c = dot(&sphere_to_ray, &sphere_to_ray) - self.radius * self.radius;
        let discriminant = (h * h) - (a * c);
        if discriminant < 0.0 {
            return None;
        }
        let discriminant_sqrt = discriminant.sqrt();
        let t0 = (-h - discriminant_sqrt) / a;
        let t1 = (-h + discriminant_sqrt) / a;
        if t0 < 0.0 && t1 < 0.0 {
            return None;
        }
        let point = if t0 < t1 { ray.at(t0) } else { ray.at(t1) };
        let surface_normal = self.surface_normal(point);
        let mut inside = false;
        if dot(&ray.direction, &surface_normal) > 0.0 {
            inside = true;
        }
        let local_normal = if inside { -surface_normal } else { surface_normal };
        let intersection = Intersection {
            point: point,
            local_normal: local_normal,
            inside: inside,
        };
        Some(intersection)
    }
}

pub enum Intersectable {
    Sphere(Sphere),
}

impl Intersects for Intersectable {
    fn surface_normal(&self, point: Vec3D) -> Vec3D {
        match self {
            Intersectable::Sphere(s) => s.surface_normal(point),
        }
    }

    fn intersects(&self, ray: &Ray) -> Option<Intersection> {
        match self {
            Intersectable::Sphere(s) => s.intersects(ray),
        }
    }
}
