use std::f32::INFINITY;
use std::option::Option;

use crate::linalg::dot;
use crate::linalg::Vec3D;

use crate::materials::Material;

use crate::ray::Ray;

#[derive(Copy, Clone)]
pub struct Intersection {
    pub point: Vec3D,
    pub distance: f32,
    pub local_normal: Vec3D,
    pub inside: bool,
}

impl Default for Intersection {
    fn default() -> Intersection {
        Intersection {
            point: Vec3D(0.0, 0.0, 0.0),
            distance: 0.0,
            local_normal: Vec3D(0.0, 0.0, 0.0),
            inside: false,
        }
    }
}

pub trait Intersects {
    fn intersects(&self, ray: &Ray) -> Option<Intersection>;
    fn surface_normal(&self, point: Vec3D) -> Vec3D;
    fn material(&self) -> Material;
}

pub struct Sphere {
    pub origin: Vec3D,
    pub radius: f32,
    pub material: Material,
}

impl Intersects for Sphere {
    fn surface_normal(&self, point: Vec3D) -> Vec3D {
        (point - self.origin) / self.radius
    }

    fn intersects(&self, ray: &Ray) -> Option<Intersection> {
        // using quadratic formula
        let sphere_to_ray = ray.origin - self.origin;
        let a = ray.direction.length_squared();
        let h = dot(&sphere_to_ray, &ray.direction);
        let c = sphere_to_ray.length_squared() - self.radius * self.radius;
        let discriminant = (h * h) - (a * c);
        if discriminant < 0.0 {
            return None;
        }
        let discriminant_sqrt = discriminant.sqrt();
        let t0 = (-h - discriminant_sqrt) / a;
        let t1 = (-h + discriminant_sqrt) / a;
        if t0 < 1e-2 && t1 < 1e-2 {
            return None;
        }
        let t = if t0 >= 1e-4 { t0 } else { t1 };
        let point = ray.at(t);
        let surface_normal = self.surface_normal(point);
        let mut inside = false;
        if dot(&ray.direction, &surface_normal) > 0.0 {
            inside = true;
        }
        let local_normal = if inside {
            -surface_normal
        } else {
            surface_normal
        };
        let distance = (point - ray.origin).length();
        let intersection = Intersection {
            point,
            distance,
            local_normal,
            inside,
        };
        Some(intersection)
    }

    fn material(&self) -> Material {
        self.material
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

    fn material(&self) -> Material {
        match self {
            Intersectable::Sphere(s) => s.material(),
        }
    }
}

pub fn find_intersections<'a, 'b>(
    ray: &'a Ray,
    world: &'b Vec<Intersectable>,
) -> Option<(Intersection, &'b Intersectable)> {
    // calculate intersection list
    let num_objects = world.len() as usize;
    let mut closest_distance = INFINITY;
    let mut closest_intersectable = &world[0];
    let mut closest_intersection = Intersection::default();
    for i in 0..num_objects {
        let result = world[i].intersects(ray);
        match result {
            Some(intersection) => {
                if intersection.distance < closest_distance {
                    closest_distance = intersection.distance;
                    closest_intersectable = &world[i];
                    closest_intersection = intersection;
                }
            }
            None => {}
        }
    }
    if closest_distance == INFINITY {
        None
    } else {
        Some((closest_intersection, closest_intersectable))
    }
}
