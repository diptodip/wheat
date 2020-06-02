use std::ops;
use std::f64::consts::PI;
use std::option::Option;

struct RGB {
    r: f64,
    g: f64,
    b: f64,
}

fn rgb(r: f64, g: f64, b: f64) -> RGB {
    RGB{r: r, g: g, b: b}
}

fn save_rgb(image: Vec<Vec<Vec<f64>>>, rows: usize, cols:usize) {
    let rows = rows;
    let cols = cols;
    println!("P3");
    print!("{}", cols);
    print!(" ");
    println!("{}", rows);
    println!("{}", 255);
    for i in 0..rows {
        for j in 0..cols {
            let r = (255.999 * image[i][j][0]) as i32;
            let g = (255.999 * image[i][j][1]) as i32;
            let b = (255.999 * image[i][j][2]) as i32;
            println!("{} {} {} ", r, g, b);
        }
    }
}

#[derive(Copy, Clone)]
struct Vector3D {
    x: f64,
    y: f64,
    z: f64,
}

fn vec3d(x: f64, y: f64, z: f64) -> Vector3D {
    Vector3D{x: x, y: y, z: z}
}

impl Vector3D {
    fn unit_vector(self) -> Vector3D {
        let norm = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        vec3d(self.x / norm, self.y / norm, self.z / norm)
    }

    fn max_norm(self) -> Vector3D {
        let norm = self.max();
        vec3d(self.x / norm, self.y / norm, self.z / norm)
    }

    fn max(self) -> f64 {
        let mut max = self.x;
        if self.y > max {
            max = self.y;
        }
        if self.z > max {
            max = self.z;
        }
        return max;
    }
}

impl ops::Add<Vector3D> for Vector3D {
    type Output = Vector3D;
    fn add(self, _rhs: Vector3D) -> Vector3D {
        return vec3d(
            self.x + _rhs.x,
            self.y + _rhs.y,
            self.z + _rhs.z,
        );
    }
}

impl ops::Add<f64> for Vector3D {
    type Output = Vector3D;
    fn add(self, _rhs: f64) -> Vector3D {
        return vec3d(
            self.x + _rhs,
            self.y + _rhs,
            self.z + _rhs,
        );
    }
}

impl ops::Add<Vector3D> for f64 {
    type Output = Vector3D;
    fn add(self, _rhs: Vector3D) -> Vector3D {
        return vec3d(
            _rhs.x + self,
            _rhs.y + self,
            _rhs.z + self,
        );
    }
}

impl ops::Sub<Vector3D> for Vector3D {
    type Output = Vector3D;
    fn sub(self, _rhs: Vector3D) -> Vector3D {
        return vec3d(
            self.x - _rhs.x,
            self.y - _rhs.y,
            self.z - _rhs.z,
        );
    }
}

impl ops::Sub<f64> for Vector3D {
    type Output = Vector3D;
    fn sub(self, _rhs: f64) -> Vector3D {
        return vec3d(
            self.x - _rhs,
            self.y - _rhs,
            self.z - _rhs,
        );
    }
}

impl ops::Sub<Vector3D> for f64 {
    type Output = Vector3D;
    fn sub(self, _rhs: Vector3D) -> Vector3D {
        return vec3d(
            _rhs.x - self,
            _rhs.y - self,
            _rhs.z - self,
        );
    }
}

impl ops::Mul<Vector3D> for Vector3D {
    type Output = Vector3D;
    fn mul(self, _rhs: Vector3D) -> Vector3D {
        return vec3d(
            self.x * _rhs.x,
            self.y * _rhs.y,
            self.z * _rhs.z,
        );
    }
}

impl ops::Mul<&Vector3D> for &Vector3D {
    type Output = Vector3D;
    fn mul(self, _rhs: &Vector3D) -> Vector3D {
        return vec3d(
            self.x * _rhs.x,
            self.y * _rhs.y,
            self.z * _rhs.z,
        );
    }
}

impl ops::Mul<f64> for Vector3D {
    type Output = Vector3D;
    fn mul(self, _rhs: f64) -> Vector3D {
        return vec3d(
            self.x * _rhs,
            self.y * _rhs,
            self.z * _rhs,
        );
    }
}

impl ops::Mul<Vector3D> for f64 {
    type Output = Vector3D;
    fn mul(self, _rhs: Vector3D) -> Vector3D {
        return vec3d(
            _rhs.x * self,
            _rhs.y * self,
            _rhs.z * self,
        );
    }
}

impl ops::Div<Vector3D> for Vector3D {
    type Output = Vector3D;
    fn div(self, _rhs: Vector3D) -> Vector3D {
        return vec3d(
            self.x / _rhs.x,
            self.y / _rhs.y,
            self.z / _rhs.z,
        );
    }
}

impl ops::Div<f64> for Vector3D {
    type Output = Vector3D;
    fn div(self, _rhs: f64) -> Vector3D {
        return vec3d(
            self.x / _rhs,
            self.y / _rhs,
            self.z / _rhs,
        )
    }
}

impl ops::Div<Vector3D> for f64 {
    type Output = Vector3D;
    fn div(self, _rhs: Vector3D) -> Vector3D {
        return vec3d(
            _rhs.x / self,
            _rhs.y / self,
            _rhs.z / self,
        )
    }
}

fn cross(a: &Vector3D, b: &Vector3D) -> Vector3D {
    vec3d(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn dot(a: &Vector3D, b: &Vector3D) -> f64 {
    let c = a * b;
    c.x + c.y + c.z
}

struct Ray {
    origin: Vector3D,
    direction: Vector3D,
}

impl Ray {
    fn at(&self, t:f64) -> Vector3D {
        t * self.direction
    }
}

#[derive(Copy, Clone)]
struct ImagePlane {
    origin: Vector3D,
    norm: Vector3D,
    up: Vector3D,
    height: f64,
    width: f64,
}

#[derive(Copy, Clone)]
struct Camera {
    origin: Vector3D,
    image_plane: ImagePlane,
}

impl Camera {
    fn ray_from_pixel(self, row_frac: f64, col_frac: f64) -> Ray {
        let height = self.image_plane.height as f64;
        let width = self.image_plane.width as f64;
        let basis1 = cross(&self.image_plane.up, &self.image_plane.norm).unit_vector();
        let basis2 = cross(&self.image_plane.norm, &basis1).unit_vector();
        let left_bottom = self.image_plane.origin - (width * basis1 / 2.0) - (height * basis2 / 2.0);
        return Ray {
            origin: self.origin,
            direction: (left_bottom
                        + (basis1 * col_frac * width)
                        + (basis2 * row_frac * height)
                        - self.origin),
        };
    }
}

struct Sphere {
    origin: Vector3D,
    radius: f64,
}

impl Sphere {
    fn intersects(&self, ray: &Ray) -> Option<Vector3D> {
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

fn fov_to_image_plane_size(fov_radians: f64, aspect_ratio: f64) -> (f64, f64) {
    let h = (fov_radians / 2.0).tan();
    let image_plane_height = 2.0 * h;
    let image_plane_width = aspect_ratio * image_plane_height;
    (image_plane_height, image_plane_width)
}

fn image_plane_size_to_fov(image_plane_height: f64, image_plane_width: f64) -> (f64, f64) {
    let aspect_ratio = image_plane_width / image_plane_height;
    let h = image_plane_height / 2.0;
    let fov = h.atan() / 2.0;
    (fov, aspect_ratio)
}

fn image_plane_size_to_pixel_shape(image_plane_height: f64, image_plane_width: f64, pixel_size: f64) -> (u32, u32, f64) {
    let rows = (image_plane_height / pixel_size) as u32;
    let cols = (image_plane_width / pixel_size) as u32;
    let pixel_size = image_plane_height / (rows as f64);
    (rows, cols, pixel_size)
}

fn trace(ray: &Ray, sphere: &Sphere) -> RGB {
    // determine if ray intersects
    let result = sphere.intersects(ray);
    match result {
        // calculate color at intersection point
        // TODO(dip): implement recursion, for now only return color
        Some(intersection) => {
            //
            // uncomment to print intersection locations
            // eprintln!("[info] found intersection at {} {} {}", intersection.x, intersection.y, intersection.z);
            let surface_normal = (intersection - sphere.origin).unit_vector();
            return rgb(0.5 * (surface_normal.x + 1.0), 0.5 * (surface_normal.y + 1.0), 0.5 * (surface_normal.z + 1.0));
        },
        None => {
            let ray_direction = ray.direction.unit_vector();
            let height = 0.5 * (ray_direction.y + 1.0);
            return rgb((1.0 - height) + height * 0.5, (1.0 - height) + height * 0.7, (1.0 - height) + height * 1.0);
        }
    }
}

fn spheres() {
    // construct camera
    let fov = (PI / 2.0) as f64;
    let aspect_ratio = 16.0 / 9.0;
    let (image_plane_height, image_plane_width) = fov_to_image_plane_size(fov, aspect_ratio);
    let (rows, cols, pixel_size) = image_plane_size_to_pixel_shape(image_plane_height, image_plane_width, 2.0 / 216.0);
    let camera = Camera {
        origin: vec3d(0.0, 0.0, 0.0),
        image_plane: ImagePlane {
            origin: vec3d(0.0, 0.0, -1.0),
            norm: vec3d(0.0, 0.0, 1.0),
            up: vec3d(0.0, 1.0, 0.0),
            height: image_plane_height,
            width: image_plane_width,
        }
    };
    let rows = rows as usize;
    let cols = cols as usize;
    let odd_rows = if rows % 2 == 0 { rows - 1 } else { rows };
    let odd_cols = if cols % 2 == 0 { cols - 1 } else { cols };
    // construct blank image
    let mut image = vec![vec![vec![0.0; 3]; cols]; rows];
    // construct sphere objects in scene
    let big_sphere = Sphere {origin: vec3d(0.0, 0.0, -1.0), radius: 0.5};
    // loop over pixels and create rays
    for row in (0..rows).rev() {
        for col in 0..cols {
            // calculate ray for current pixel
            let ray = camera.ray_from_pixel((row as f64) / (odd_rows as f64), (col as f64) / (odd_cols as f64));
            // trace ray for current pixel
            let color = trace(&ray, &big_sphere);
            // add observed color from trace to image at current pixel
            image[rows - row - 1][col][0] = color.r;
            image[rows - row - 1][col][1] = color.g;
            image[rows - row - 1][col][2] = color.b;
        }
        eprintln!("[info] remaining scan lines: {}", row);
    }
    eprintln!("[info] saving image...");
    save_rgb(image, rows, cols);
    eprintln!("[ok] done!");
}

fn main() {
    spheres();
}
