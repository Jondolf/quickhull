// Adapted primarily from `parry3d`:
// <https://github.com/dimforge/parry/blob/9db68641adf69e1f307ac9199d34d82b6d049219/src/transformation/convex_hull3/convex_hull.rs>

mod error;
mod hull;
mod initial_hull;
mod normalize;
mod plane;
mod triangle_mesh;

pub use error::ConvexHull3dError;
pub use hull::{ConvexHull3d, ConvexHull3dSettings, HullFace};
pub use plane::Plane3d;
pub use triangle_mesh::ConvexTriangleMesh;
