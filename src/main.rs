use std::f64::consts::PI;

mod linalg;
use linalg::Vec3D;

mod colors;
use colors::RGB;
use colors::rgb;

mod io;
use io::output_ppm;

mod geometry;
use geometry::Ray;
use geometry::Sphere;
use geometry::Intersects;
use geometry::Intersectable;

mod camera;
use camera::Camera;
use camera::ImagePlane;
use camera::fov_to_imsize;
use camera::imsize_to_pixels;

fn trace(ray: &Ray, intersectable: &Intersectable) -> RGB {
    // determine if ray intersects
    let result = intersectable.intersects(ray);
    match result {
        // calculate color at intersection point
        // TODO(dip): implement recursion, for now only return color
        Some(intersection) => {
            // uncomment to print intersection locations
            // eprintln!("[info] found intersection at {} {} {}",
            //           intersection.point.0,
            //           intersection.point.1,
            //           intersection.point.2);
            let local_normal = intersection.local_normal;
            let rescaled_normal = 0.5 * (local_normal + 1.0);
            return rgb(rescaled_normal.0, rescaled_normal.1, rescaled_normal.2);
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
                                                    2.0 / 480.0);
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
    let big_sphere = Intersectable::Sphere(Sphere {origin: Vec3D(0.0, 0.0, -1.0), radius: 0.5});
    // loop over pixels and create rays
    eprintln!("[start] processing {}px x {}px (width x height)...", cols, rows);
    eprintln!("[info] remaining scan lines: {}", rows);
    for row in 0..rows {
        for col in 0..cols {
            // calculate ray for current pixel
            // making sure to center ray within pixel
            let ray = camera.prime_ray((row as f64 + 0.5) / (rows as f64),
                                       (col as f64 + 0.5) / (cols as f64));
            // trace ray for current pixel
            let color = trace(&ray, &big_sphere);
            // add observed color from trace to image at current pixel
            image[row][col][0] = color.r();
            image[row][col][1] = color.g();
            image[row][col][2] = color.b();
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
