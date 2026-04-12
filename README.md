# Quickhull

A Rust implementation of the Quickhull algorithm for computing [convex hulls] for 2D and 3D point sets.

![The Utah teapot and its convex hull](./images/utah_teapot.png)

*The Utah teapot and its convex hull*

[convex hulls]: https://en.wikipedia.org/wiki/Convex_hull

## Features

- `ConvexHull2d` for computing 2D convex hulls
- `ConvexTriangleMesh` for computing 3D convex hulls with triangular faces
- `ConvexHull3d` for computing 3D convex hulls with polygonal faces

## Implementation

`ConvexHull2d` is based on [`geo`], and uses the [`robust`] crate for robust geometric predicates.
It should handle arbitrary point clouds without failure.

`ConvexTriangleMesh` is largely adapted from [`parry3d`], with some modifications. It is designed
to be very fast, while handling degenerate cases such as coplanar inputs reliably.

`ConvexHull3d` builds on top of `ConvexTriangleMesh` by merging coplanar triangular faces
into polygonal faces with an arbitrary number of vertices. It also computes the associated
bounding planes, which can be useful for geometry applications. The hull representation
was inspired by [Jolt Physics].

Only `f32` vector types from [`glam`] are currently supported.

[`geo`]: https://github.com/georust/geo/blob/8940db79aa6aa4ec8820d6328c68a2ae08ac8fdc/geo/src/algorithm/convex_hull/qhull.rs
[`robust`]: https://github.com/georust/robust
[`parry3d`]: https://github.com/dimforge/parry/blob/9db68641adf69e1f307ac9199d34d82b6d049219/src/transformation/convex_hull3/convex_hull.rs
[Jolt Physics]: https://jrouwe.github.io/JoltPhysics/
[`glam`]: https://docs.rs/glam/latest/glam/

## 3D Example

```rust
use glam::Vec3A;
use quickhull::{ConvexHull3d, ConvexHull3dSettings};

// Define points representing the corners of a cube.
let points = vec![
    Vec3A::new( 1.0,  1.0,  1.0),
    Vec3A::new( 1.0,  1.0, -1.0),
    Vec3A::new( 1.0, -1.0,  1.0),
    Vec3A::new( 1.0, -1.0, -1.0),
    Vec3A::new(-1.0,  1.0,  1.0),
    Vec3A::new(-1.0,  1.0, -1.0),
    Vec3A::new(-1.0, -1.0,  1.0),
    Vec3A::new(-1.0, -1.0, -1.0),
    // Additional internal points that should not affect the hull
    Vec3A::new(0.0, 0.0, 0.0),
    Vec3A::new(0.5, 0.5, 0.5),
];

let settings = ConvexHull3dSettings {
    // Merge faces with normals within 3 degrees of each other
    coplanarity_dot_tolerance: 0.9986295, // cos(3 degrees)
    // Allow up to 10,000 iterations for hull construction (should not be hit in practice)
    max_iterations: 10_000,
};

// Compute the convex hull.
let hull = ConvexHull3d::try_from_points(&points, settings).unwrap();

// The cube's hull should have 8 vertices and 6 quad faces.
assert_eq!(hull.vertices().len(), 8);
assert_eq!(hull.faces().len(), 6);
```

## References

- The [`geo`] crate
- The [`parry3d`] crate
- [Jolt Physics]
- C. Bradford Barber et al. 1996. [The Quickhull Algorithm for Convex Hulls](https://www.cise.ufl.edu/~ungor/courses/fall06/papers/QuickHull.pdf) (the original paper)
- Dirk Gregorius. GDC 2014. [Physics for Game Programmers: Implementing Quickhull](https://archive.org/details/GDC2014Gregorius)

## License

This Quickhull crate is free and open source. All code in this repository is dual-licensed under either:

- MIT License ([LICENSE-MIT](/LICENSE-MIT) or <http://opensource.org/licenses/MIT>)
- Apache License, Version 2.0 ([LICENSE-APACHE](/LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)

at your option.
