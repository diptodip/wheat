use std::f64::consts::PI;

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
use geometry::Ray;
use geometry::Sphere;
use geometry::Intersects;
use geometry::Intersection;
use geometry::Intersectable;
use geometry::first_intersection;

mod camera;
use camera::Camera;
use camera::ImagePlane;
use camera::fov_to_imsize;
use camera::imsize_to_pixels;

mod materials;
use materials::Material;
use materials::Diffuse;

fn diffuse_bounce(intersection: &Intersection,
                  intersectable: &Intersectable) -> Ray {
    let bounce_vector = intersection.local_normal + Sphere::random_unit_vector();
    Ray {
        origin: intersection.point,
        direction: bounce_vector,
    }
}

fn reflected_bounce(intersection: &Intersection, intersectable: &Intersectable, ray: &Ray) -> Ray {
    let normal_perpendicular_component = ray.direction - intersection.local_normal;
    let normal_parallel_component = ray.direction - normal_perpendicular_component;
    Ray {
        origin: intersection.point,
        direction: -normal_perpendicular_component + normal_parallel_component,
    }
}

fn trace(ray: &Ray, world: &Vec<&Intersectable>, depth: u64) -> RGB {
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
        // TODO(dip): calculate color using material color
        Some((intersection, intersectable)) => {
            // uncomment to print intersection locations
            // eprintln!("[info] found intersection at {} {} {}",
            //           intersection.point.0,
            //           intersection.point.1,
            //           intersection.point.2);
            let mut rng = rand::thread_rng();
            if let Material::Diffuse(d) = intersectable.material() {
                // light bounces if material is diffuse,
                // so we recurse and trace a bounced ray
                let bounce_ray = &diffuse_bounce(&intersection, intersectable);
                let color_vec =  d.color.to_vec3d() * trace(bounce_ray, world, depth - 1).to_vec3d();
                return vec_to_rgb(color_vec);
            }
            // TODO(dip): handle other materials
            return rgb(0.0, 0.0, 0.0);
        },
        None => {
            let ray_direction = ray.direction.l2_normalize();
            let height = 0.5 * (ray_direction.1 + 1.0);
            return rgb((1.0 - height) + height * 0.5,
                       (1.0 - height) + height * 0.7,
                       (1.0 - height) + height * 1.0);
        }
    }
}

fn spheres() {
    // construct camera
    let fov = PI / 2.0;
    let aspect_ratio = 16.0 / 9.0;
    let (image_plane_height, image_plane_width) = fov_to_imsize(fov,
                                                                aspect_ratio);
    let (rows, cols, pixel_size) = imsize_to_pixels(image_plane_height,
                                                    image_plane_width,
                                                    2.0 / 216.0);
    let camera = Camera {
        origin: Vec3D(0.0, 0.0, 0.0),
        image_plane: ImagePlane {
            origin: Vec3D(0.0, 0.0, -1.0),
            norm: Vec3D(0.0, 0.0, 1.0),
            up: Vec3D(0.0, 1.0, 0.0),
            height: image_plane_height,
            width: image_plane_width,
        }
    };
    let rows = rows as usize;
    let cols = cols as usize;
    // construct blank image
    let mut image = vec![vec![vec![0.0; 3]; cols]; rows];
    // construct sphere objects in scene
    let big_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, 0.0, -1.0),
        radius: 0.5,
        material: Material::Diffuse(Diffuse { color: rgb(0.7, 0.3, 0.3) }),
    });
    let ground_sphere = Intersectable::Sphere(Sphere {
        origin: Vec3D(0.0, -100.5, -1.0),
        radius: 100.0,
        material: Material::Diffuse(Diffuse { color: rgb(0.8, 0.8, 0.0) }),
    });
    // create pointer to list of objects, aka the world
    let world = &vec![&big_sphere, &ground_sphere];
    // loop over pixels and create rays
    eprintln!("[start] processing {}px x {}px (width x height)...", cols, rows);
    eprintln!("[info] remaining scan lines: {}", rows);
    let samples_per_pixel = 100.0;
    let mut rng = rand::thread_rng();
    for row in 0..rows {
        for col in 0..cols {
            for sample in 0..samples_per_pixel as usize {
                // calculate ray for current pixel
                // making sure to center ray within pixel
                // we also perturb the ray direction slightly per sample
                let row_frac = (row as f64 + 0.5 + rng.gen::<f64>()) / (rows as f64);
                let col_frac = (col as f64 + 0.5 + rng.gen::<f64>()) / (cols as f64);
                let ray = &camera.prime_ray(row_frac, col_frac);
                // trace ray for current pixel
                let color = trace(ray, world, 50);
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
    spheres();
}
