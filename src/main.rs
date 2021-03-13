extern crate rand;
use rand::prelude::*;

mod linalg;
use linalg::Vec3D;
use linalg::dot;

mod colors;
use colors::RGB;
use colors::rgb;
use colors::vec_to_rgb;

mod io;
use io::output_ppm;

mod geometry;
use geometry::Sphere;
use geometry::Intersects;
use geometry::Intersection;
use geometry::Intersectable;
use geometry::first_intersection;

mod camera;
use camera::Camera;
use camera::fov_to_imsize;
use camera::imsize_to_pixels;

mod materials;
use materials::Material;
use materials::Surface;

mod ray;
use ray::*;

fn test_spheres() {
    // construct camera
    let fov = 20.0_f64.to_radians();
    let aspect_ratio = 16.0 / 9.0;
    let pixel_height = 216.0;
    let (image_plane_height, image_plane_width) = fov_to_imsize(fov,
                                                                aspect_ratio);
    let (rows, cols, pixel_size) = imsize_to_pixels(image_plane_height,
                                                    image_plane_width,
                                                    image_plane_height / pixel_height);
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
    // construct blank image
    let mut image = vec![vec![vec![0.0; 3]; cols]; rows];
    // initialize list of objects, aka the world
    let mut world = Vec::new();
    // construct ground sphere in scene
    let ground = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, -100.5, -1.0),
        radius: 100.0,
        material: Material {
            color: rgb(0.8, 0.8, 0.0),
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
            surface: Surface::Refractive(1.5)
        },
    });
    let small_glass_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(-1.0, 0.0, -1.0),
        radius: -0.45,
        material: Material {
            color: rgb(1.0, 1.0, 1.0),
            surface: Surface::Refractive(1.5),
        },
    });
    let big_diffuse_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, 0.0, -1.0),
        radius: 0.5,
        material: Material {
            color: rgb(0.1, 0.2, 0.5),
            surface: Surface::Diffuse,
        },
    });
    let big_reflective_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(1.0, 0.0, -1.0),
        radius: 0.5,
        material: Material {
            color: rgb(0.8, 0.6, 0.2),
            surface: Surface::FuzzyReflective(0.3),
        },
    });
    world.push(big_glass_sphere);
    world.push(small_glass_sphere);
    world.push(big_diffuse_sphere);
    world.push(big_reflective_sphere);
    // loop over pixels and create rays
    let mut rng = rand::thread_rng();
    eprintln!("[start] processing {}px x {}px (width x height)...", cols, rows);
    eprintln!("[info] remaining scan lines: {}", rows);
    let samples_per_pixel = 100.0;
    for row in 0..rows {
        for col in 0..cols {
            for sample in 0..samples_per_pixel as usize {
                // calculate ray for current pixel
                // making sure to center ray within pixel
                // we also perturb the ray direction slightly per sample
                let row_rand = rng.gen::<f64>();
                let col_rand = rng.gen::<f64>();
                let row_frac = (row as f64 + 0.5 + row_rand) / (rows as f64);
                let col_frac = (col as f64 + 0.5 + col_rand) / (cols as f64);
                let ray = camera.prime_ray(row_frac, col_frac);
                // trace ray for current pixel
                let color = trace(&ray, &world, 50);
                // add observed color from trace to image at current pixel
                image[row][col][0] += (color.r() / samples_per_pixel);
                image[row][col][1] += (color.g() / samples_per_pixel);
                image[row][col][2] += (color.b() / samples_per_pixel);
            }
        }
        eprintln!("[info] remaining scan lines: {}", rows - row - 1);
    }
    eprintln!("[info] saving image...");
    output_ppm(image, rows, cols);
    eprintln!("[ok] done!");
}

fn random_spheres() {
    // construct camera
    let fov = 20.0_f64.to_radians();
    let aspect_ratio = 1.5;
    let pixel_height = 300.0;
    let (image_plane_height, image_plane_width) = fov_to_imsize(fov,
                                                                aspect_ratio);
    let (rows, cols, pixel_size) = imsize_to_pixels(image_plane_height,
                                                    image_plane_width,
                                                    image_plane_height / pixel_height);
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
    // construct blank image
    let mut image = vec![vec![vec![0.0; 3]; cols]; rows];
    // initialize list of objects, aka the world
    let mut world = Vec::new();
    // construct ground sphere in scene
    let ground = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, -1000.0, 0.0),
        radius: 1000.0,
        material: Material {
            color: rgb(0.5, 0.5, 0.5),
            surface: Surface::Diffuse,
        },
    });
    world.push(ground);
    // construct random small spheres in scene
    let mut rng = rand::thread_rng();
    for i in -11..11 {
        for j in -11..11 {
            let material_check = rng.gen::<f64>();
            let center = Vec3D(i as f64 + 0.9 * rng.gen::<f64>(),
                               0.2,
                               j as f64 + 0.9 * rng.gen::<f64>());
            if (center - Vec3D(4.0, 0.2, 0.0)).length() > 0.9 {
                if material_check < 0.8 {
                    // make diffuse sphere
                    let color = vec_to_rgb(RGB::random().to_vec3d()
                                           * RGB::random().to_vec3d());
                    let diffuse_sphere = Intersectable::Sphere(Sphere {
                        origin: center,
                        radius: 0.2,
                        material: Material {
                            color: color,
                            surface: Surface::Diffuse,
                        },
                    });
                    world.push(diffuse_sphere);
                } else if material_check < 0.95 {
                    // make fuzzy reflective sphere
                    let color = vec_to_rgb(Vec3D::random(0.5, 1.0,
                                                         0.5, 1.0,
                                                         0.5, 1.0));
                    let fuzz: f64 = rng.gen_range(0.0, 0.5);
                    let reflective_sphere = Intersectable::Sphere(Sphere {
                        origin: center,
                        radius: 0.2,
                        material: Material {
                            color: color,
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
            surface: Surface::Refractive(1.5),
        },
    });
    let big_diffuse_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(-4.0, 1.0, 0.0),
        radius: 1.0,
        material: Material {
            color: rgb(0.4, 0.2, 0.1),
            surface: Surface::Diffuse,
        },
    });
    let big_reflective_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(4.0, 1.0, 0.0),
        radius: 1.0,
        material: Material {
            color: rgb(0.7, 0.6, 0.5),
            surface: Surface::Reflective,
        },
    });
    world.push(big_glass_sphere);
    world.push(big_diffuse_sphere);
    world.push(big_reflective_sphere);
    // loop over pixels and create rays
    eprintln!("[start] processing {}px x {}px (width x height)...", cols, rows);
    eprintln!("[info] remaining scan lines: {}", rows);
    let samples_per_pixel = 100.0;
    for row in 0..rows {
        for col in 0..cols {
            for sample in 0..samples_per_pixel as usize {
                // calculate ray for current pixel
                // making sure to center ray within pixel
                // we also perturb the ray direction slightly per sample
                let row_rand = rng.gen::<f64>();
                let col_rand = rng.gen::<f64>();
                let row_frac = (row as f64 + 0.5 + row_rand) / (rows as f64);
                let col_frac = (col as f64 + 0.5 + col_rand) / (cols as f64);
                let ray = camera.prime_ray(row_frac, col_frac);
                // trace ray for current pixel
                let color = trace(&ray, &world, 50);
                // add observed color from trace to image at current pixel
                image[row][col][0] += (color.r() / samples_per_pixel);
                image[row][col][1] += (color.g() / samples_per_pixel);
                image[row][col][2] += (color.b() / samples_per_pixel);
            }
        }
        eprintln!("[info] remaining scan lines: {}", rows - row - 1);
    }
    eprintln!("[info] saving image...");
    output_ppm(image, rows, cols);
    eprintln!("[ok] done!");
}

fn main() {
    random_spheres();
}
