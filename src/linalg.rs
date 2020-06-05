use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Neg};
use rand::prelude::*;

#[derive(Copy, Clone)]
pub struct Vec3D(pub f64, pub f64, pub f64);

impl Vec3D {
    #[inline]
    pub fn length(self) -> f64 {
        (self.0 * self.0 + self.1 * self.1 + self.2 * self.2).sqrt()
    }
    #[inline]
    fn normalize(self, norm: f64) -> Vec3D {
        Vec3D(self.0 / norm, self.1 / norm, self.2 / norm)
    }

    pub fn l2_normalize(self) -> Vec3D {
        let norm = self.length();
        self.normalize(norm)
    }

    pub fn max_normalize(self) -> Vec3D {
        let mut max = self.0;
        if self.1 > max {
            max = self.1
        }
        if self.2 > max {
            max = self.2
        }
        let norm = max;
        max += 1e-12;
        self.normalize(norm)
    }
}

macro_rules! implement_binary_operation {
    ($trait: ident, $trait_fn: ident, $op: tt) => {
        impl $trait<Vec3D> for Vec3D {
            type Output = Vec3D;
            fn $trait_fn(self, _rhs: Vec3D) -> Vec3D {
                Vec3D(self.0 $op _rhs.0, self.1 $op _rhs.1, self.2 $op _rhs.2)
            }
        }

        impl $trait<&Vec3D> for &Vec3D {
            type Output = Vec3D;
            fn $trait_fn(self, _rhs: &Vec3D) -> Vec3D {
                Vec3D(self.0 $op _rhs.0, self.1 $op _rhs.1, self.2 $op _rhs.2)
            }
        }

        impl $trait<f64> for Vec3D {
            type Output = Vec3D;
            fn $trait_fn(self, _rhs: f64) -> Vec3D {
                Vec3D(self.0 $op _rhs, self.1 $op _rhs, self.2 $op _rhs)
            }
        }

        impl $trait<Vec3D> for f64 {
            type Output = Vec3D;
            fn $trait_fn(self, _rhs: Vec3D) -> Vec3D {
                Vec3D(self $op _rhs.0, self $op _rhs.1, self $op _rhs.2)
            }
        }

        impl $trait<f64> for &Vec3D {
            type Output = Vec3D;
            fn $trait_fn(self, _rhs: f64) -> Vec3D {
                Vec3D(self.0 $op _rhs, self.1 $op _rhs, self.2 $op _rhs)
            }
        }

        impl $trait<&Vec3D> for f64 {
            type Output = Vec3D;
            fn $trait_fn(self, _rhs: &Vec3D) -> Vec3D {
                Vec3D(self $op _rhs.0, self $op _rhs.1, self $op _rhs.2)
            }
        }
    };
}

macro_rules! implement_assign_operation {
    ($trait: ident, $trait_fn: ident, $op: tt) => {
        impl $trait<Vec3D> for Vec3D {
            fn $trait_fn(&mut self, _rhs: Vec3D) {
                self.0 $op _rhs.0;
                self.1 $op _rhs.1;
                self.2 $op _rhs.2;
            }
        }

        impl $trait<f64> for Vec3D {
            fn $trait_fn(&mut self, _rhs: f64) {
                self.0 $op _rhs;
                self.1 $op _rhs;
                self.2 $op _rhs;
            }
        }
    };
}

macro_rules ! implement_unary_operation {
    ($trait: ident, $trait_fn: ident, $op: tt) => {
        impl $trait for Vec3D {
            type Output = Vec3D;
            fn $trait_fn(self) -> Vec3D {
                Vec3D($op self.0, $op self.1, $op self.2)
            }
        }
    };
}

// binary operations: + - * /
implement_binary_operation!(Add, add, +);
implement_binary_operation!(Sub, sub, -);
implement_binary_operation!(Mul, mul, *);
implement_binary_operation!(Div, div, /);

// assign operations: += -= *= /=
implement_assign_operation!(AddAssign, add_assign, +=);
implement_assign_operation!(SubAssign, sub_assign, -=);
implement_assign_operation!(MulAssign, mul_assign, *=);
implement_assign_operation!(DivAssign, div_assign, /=);

// unary operation: -
implement_unary_operation!(Neg, neg, -);

pub fn cross(a: &Vec3D, b: &Vec3D) -> Vec3D {
    Vec3D(
        a.1 * b.2 - a.2 * b.1,
        a.2 * b.0 - a.0 * b.2,
        a.0 * b.1 - a.1 * b.0,
    )
}

pub fn dot(a: &Vec3D, b: &Vec3D) -> f64 {
    let c = a * b;
    c.0 + c.1 + c.2
}
