//! # Quickhull
//!
//! A Rust implementation of the Quickhull algorithm for computing [convex hulls] for 2D and 3D point sets.
//!
//! [convex hulls]: https://en.wikipedia.org/wiki/Convex_hull
//!
//! ## Features
//!
//! - [`ConvexHull2d`] for computing 2D convex hulls
//! - [`ConvexTriangleMesh`] for computing 3D convex hulls with triangular faces
//! - [`ConvexHull3d`] for computing 3D convex hulls with polygonal faces
//!
//! ## References
//!
//! - The [`geo`] crate
//! - The [`parry3d`] crate
//! - [Jolt Physics]
//! - C. Bradford Barber et al. 1996. [The Quickhull Algorithm for Convex Hulls](https://www.cise.ufl.edu/~ungor/courses/fall06/papers/QuickHull.pdf) (the original paper)
//! - Dirk Gregorius. GDC 2014. [Physics for Game Programmers: Implementing Quickhull](https://archive.org/details/GDC2014Gregorius)
//!
//! [`geo`]: https://github.com/georust/geo/blob/8940db79aa6aa4ec8820d6328c68a2ae08ac8fdc/geo/src/algorithm/convex_hull/qhull.rs
//! [`parry3d`]: https://github.com/dimforge/parry/blob/9db68641adf69e1f307ac9199d34d82b6d049219/src/transformation/convex_hull3/convex_hull.rs
//! [Jolt Physics]: https://jrouwe.github.io/JoltPhysics/

#![warn(missing_docs)]

mod dim2;
mod dim3;

pub(crate) mod collections;

pub use dim2::ConvexHull2d;
pub use dim3::{
    ConvexHull3d, ConvexHull3dError, ConvexHull3dSettings, ConvexTriangleMesh, HullFace, Plane3d,
};
