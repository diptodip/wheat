extern crate rand;

use rand::prelude::*;

mod linalg;
use linalg::dot;
use linalg::Vec3D;

mod colors;
use colors::rgb;
use colors::vec_to_rgb;
use colors::RGB;

mod io;
use io::write_ppm;

mod geometry;
use geometry::find_intersections;
use geometry::Intersectable;
use geometry::Intersection;
use geometry::Intersects;
use geometry::Sphere;

mod camera;
use camera::fov_to_imsize;
use camera::imsize_to_pixels;
use camera::Camera;

mod materials;
use materials::Material;
use materials::Surface;

mod ray;
use ray::*;

fn test_spheres() {
    // construct camera
    let fov = 20.0_f32.to_radians();
    let aspect_ratio = 16.0 / 9.0;
    let pixel_height = 216.0;
    let (image_plane_height, image_plane_width) = fov_to_imsize(fov, aspect_ratio);
    let (rows, cols, pixel_size) = imsize_to_pixels(
        image_plane_height,
        image_plane_width,
        image_plane_height / pixel_height,
    );
    let origin = Vec3D(-2.0, 2.0, 1.0);
    let target = Vec3D(0.0, 0.0, -1.0);
    let up = Vec3D(0.0, 1.0, 0.0);
    let focal_distance = (origin - target).length();
    let aperture = 0.1;
    let camera = Camera {
        origin: origin,
        norm: (origin - target).l2_normalize(),
        up: up,
        height: image_plane_height,
        width: image_plane_width,
        aperture: aperture,
        focal_distance: focal_distance,
    };
    let rows = rows as usize;
    let cols = cols as usize;
    let background = rgb(0.5, 0.7, 1.0);
    // initialize list of objects, aka the world
    let mut world = Vec::new();
    // construct ground sphere in scene
    let ground = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, -100.5, -1.0),
        radius: 100.0,
        material: Material {
            color: rgb(0.8, 0.8, 0.0),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::Diffuse,
        },
    });
    world.push(ground);
    // construct 3 large spheres in scene center
    let big_glass_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(-1.0, 0.0, -1.0),
        radius: 0.5,
        material: Material {
            color: rgb(1.0, 1.0, 1.0),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::Refractive(1.5),
        },
    });
    let small_glass_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(-1.0, 0.0, -1.0),
        radius: -0.4,
        material: Material {
            color: rgb(1.0, 1.0, 1.0),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::Refractive(1.5),
        },
    });
    let big_diffuse_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, 0.0, -1.0),
        radius: 0.5,
        material: Material {
            color: rgb(0.1, 0.2, 0.5),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::Diffuse,
        },
    });
    let big_reflective_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(1.0, 0.0, -1.0),
        radius: 0.5,
        material: Material {
            color: rgb(0.8, 0.6, 0.2),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::FuzzyReflective(0.3),
        },
    });
    world.push(big_glass_sphere);
    world.push(small_glass_sphere);
    world.push(big_diffuse_sphere);
    world.push(big_reflective_sphere);
    render(&world, &camera, background, rows, cols, 100.0);
}

fn random_spheres() {
    // construct camera
    let fov = 20.0_f32.to_radians();
    let aspect_ratio = 16.0 / 9.0;
    let pixel_height = 360.0;
    let (image_plane_height, image_plane_width) = fov_to_imsize(fov, aspect_ratio);
    let (rows, cols, pixel_size) = imsize_to_pixels(
        image_plane_height,
        image_plane_width,
        image_plane_height / pixel_height,
    );
    let origin = Vec3D(13.0, 2.0, 3.0);
    let target = Vec3D(0.0, 0.0, 0.0);
    let up = Vec3D(0.0, 1.0, 0.0);
    let focal_distance = 10.0;
    let aperture = 0.1;
    let camera = Camera {
        origin: origin,
        norm: (origin - target).l2_normalize(),
        up: up,
        height: image_plane_height,
        width: image_plane_width,
        aperture: aperture,
        focal_distance: focal_distance,
    };
    let rows = rows as usize;
    let cols = cols as usize;
    let background = rgb(0.5, 0.7, 1.0);
    // initialize list of objects, aka the world
    let mut world = Vec::new();
    // construct ground sphere in scene
    let ground = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, -1000.0, 0.0),
        radius: 1000.0,
        material: Material {
            color: rgb(0.5, 0.5, 0.5),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::Diffuse,
        },
    });
    world.push(ground);
    // construct random small spheres in scene
    let mut rng = thread_rng();
    for i in -11..11 {
        for j in -11..11 {
            let material_check = rng.gen::<f32>();
            let center = Vec3D(
                i as f32 + 0.9 * rng.gen::<f32>(),
                0.2,
                j as f32 + 0.9 * rng.gen::<f32>(),
            );
            if (center - Vec3D(4.0, 0.2, 0.0)).length() > 0.9 {
                if material_check < 0.8 {
                    // make diffuse sphere
                    let color = vec_to_rgb(RGB::random().to_vec3d() * RGB::random().to_vec3d());
                    let diffuse_sphere = Intersectable::Sphere(Sphere {
                        origin: center,
                        radius: 0.2,
                        material: Material {
                            color: color,
                            emit: rgb(0.0, 0.0, 0.0),
                            surface: Surface::Diffuse,
                        },
                    });
                    world.push(diffuse_sphere);
                } else if material_check < 0.95 {
                    // make fuzzy reflective sphere
                    let color = vec_to_rgb(Vec3D::random(&mut rng, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0));
                    let fuzz: f32 = rng.gen_range(0.0, 0.5);
                    let reflective_sphere = Intersectable::Sphere(Sphere {
                        origin: center,
                        radius: 0.2,
                        material: Material {
                            color: color,
                            emit: rgb(0.0, 0.0, 0.0),
                            surface: Surface::FuzzyReflective(fuzz),
                        },
                    });
                    world.push(reflective_sphere);
                } else {
                    // make refractive sphere
                    let color = rgb(1.0, 1.0, 1.0);
                    let refractive_sphere = Intersectable::Sphere(Sphere {
                        origin: center,
                        radius: 0.2,
                        material: Material {
                            color: color,
                            emit: rgb(0.0, 0.0, 0.0),
                            surface: Surface::Refractive(1.5),
                        },
                    });
                    world.push(refractive_sphere);
                }
            }
        }
    }
    // construct 3 large spheres in scene center
    let big_glass_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, 1.0, 0.0),
        radius: 1.0,
        material: Material {
            color: rgb(1.0, 1.0, 1.0),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::Refractive(1.5),
        },
    });
    let big_diffuse_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(-4.0, 1.0, 0.0),
        radius: 1.0,
        material: Material {
            color: rgb(0.4, 0.2, 0.1),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::Diffuse,
        },
    });
    let big_reflective_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(4.0, 1.0, 0.0),
        radius: 1.0,
        material: Material {
            color: rgb(0.7, 0.6, 0.5),
            emit: rgb(0.0, 0.0, 0.0),
            surface: Surface::Reflective,
        },
    });
    world.push(big_glass_sphere);
    world.push(big_diffuse_sphere);
    world.push(big_reflective_sphere);
    render(&world, &camera, background, rows, cols, 100.0);
}

fn main() {
    // random_spheres();
    test_spheres();
}
