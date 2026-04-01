//! Planes in 2D and 3D.

use glam::{Vec3A, Vec4};
use thiserror::Error;

/// Error returned when trying to create a plane with a non-unit length normal.
#[derive(Error, Debug, PartialEq)]
#[error("The normal vector must be of unit length.")]
pub struct UnnormalizedNormalError;

/// A plane in 3D defined by a normal and a signed offset from the origin.
///
/// The plane consists of all points `p` such that `dot(normal, point) == offset`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane3d {
    /// The normal vector and a signed offset from the origin.
    normal_d: Vec4,
}

impl Plane3d {
    /// Creates a new [`Plane3d`] from an outward facing normal and a signed offset from the origin.
    ///
    /// The normal vector must be of unit length.
    ///
    /// # Panics
    ///
    /// Panics if the normal vector is not of unit length when debug assertions are enabled.
    #[inline]
    pub const fn new(normal: Vec3A, offset: f32) -> Self {
        let [x, y, z] = normal.to_array();

        debug_assert!(
            (x * x + y * y + z * z - 1.0).abs() < 1e-6,
            "The normal vector must be of unit length."
        );

        Self {
            normal_d: Vec4::new(x, y, z, offset),
        }
    }

    /// Tries to create a new [`Plane3d`] from an outward facing normal and a signed offset from the origin.
    ///
    /// The normal vector must be of unit length.
    ///
    /// # Errors
    ///
    /// Returns an [`UnnormalizedNormalError`] if the normal vector is not of unit length.
    #[inline]
    pub fn try_new(normal: Vec3A, offset: f32) -> Result<Self, UnnormalizedNormalError> {
        if (normal.length_squared() - 1.0).abs() >= 1e-6 {
            return Err(UnnormalizedNormalError);
        }
        Ok(Self {
            normal_d: Vec4::new(normal.x, normal.y, normal.z, offset),
        })
    }

    /// Creates a new [`Plane3d`] from a `Vec4` where the `x`, `y`, and `z` components represent the normal
    /// and the `w` component represents the offset.
    ///
    /// # Panics
    ///
    /// Panics if the normal vector `(x, y, z)` is not of unit length when debug assertions are enabled.
    #[inline]
    pub const fn from_vec4(normal_offset: Vec4) -> Self {
        let [x, y, z, w] = normal_offset.to_array();
        Self::from_coefficients(x, y, z, w)
    }

    /// Tries to create a new [`Plane3d`] from a `Vec4` where the `x`, `y`, and `z` components represent the normal
    /// and the `w` component represents the offset.
    ///
    /// # Errors
    ///
    /// Returns an [`UnnormalizedNormalError`] if the normal vector `(x, y, z)` is not of unit length.
    #[inline]
    pub const fn try_from_vec4(normal_offset: Vec4) -> Result<Self, UnnormalizedNormalError> {
        let [x, y, z, w] = normal_offset.to_array();
        Self::try_from_coefficients(x, y, z, w)
    }

    /// Creates a new [`Plane3d`] from a point on the plane and an outward facing normal.
    #[inline]
    pub fn from_point_and_normal(point: Vec3A, normal: Vec3A) -> Self {
        let offset = -normal.dot(point);
        Self::new(normal, offset)
    }

    /// Creates a new [`Plane3d`] from the coefficients of the plane equation `ax + by + cz + d = 0`.
    ///
    /// The normal vector `(a, b, c)` must be of unit length.
    ///
    /// # Panics
    ///
    /// Panics if the normal vector `(a, b, c)` is not of unit length when debug assertions are enabled.
    #[inline]
    pub const fn from_coefficients(a: f32, b: f32, c: f32, d: f32) -> Self {
        debug_assert!(
            (a * a + b * b + c * c - 1.0).abs() < 1e-6,
            "The normal vector (a, b, c) must be of unit length."
        );

        Self {
            normal_d: Vec4::new(a, b, c, d),
        }
    }

    /// Tries to create a new [`Plane3d`] from the coefficients of the plane equation `ax + by + cz + d = 0`.
    ///
    /// The normal vector `(a, b, c)` must be of unit length.
    ///
    /// # Errors
    ///
    /// Returns an [`UnnormalizedNormalError`] if the normal vector `(a, b, c)` is not of unit length.
    #[inline]
    pub const fn try_from_coefficients(
        a: f32,
        b: f32,
        c: f32,
        d: f32,
    ) -> Result<Self, UnnormalizedNormalError> {
        let len_sq = a * a + b * b + c * c;
        if (len_sq - 1.0).abs() >= 1e-6 {
            return Err(UnnormalizedNormalError);
        }

        Ok(Self {
            normal_d: Vec4::new(a, b, c, d),
        })
    }

    /// Returns the outward facing normal of the plane.
    #[inline]
    pub fn normal(&self) -> Vec3A {
        Vec3A::from_vec4(self.normal_d)
    }

    /// Returns the signed offset from the origin.
    #[inline]
    pub fn offset(&self) -> f32 {
        self.normal_d.w
    }

    /// Returns the normal vector and signed offset as a `Vec4`.
    #[inline]
    pub const fn as_vec4(&self) -> Vec4 {
        self.normal_d
    }

    /// Computes the signed distance from the plane to a point.
    #[inline]
    pub fn signed_distance_to_point(&self, point: Vec3A) -> f32 {
        self.normal().dot(point) + self.offset()
    }
}
