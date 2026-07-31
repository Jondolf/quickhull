use approx::relative_ne;
use glam::{Vec2, Vec3, Vec3A};
use glam_matrix_extras::{SymmetricEigen3, SymmetricMat3};
use thiserror::Error;

#[cfg(debug_assertions)]
use super::triangle_mesh::validation::validate_face_connectivity;
use super::triangle_mesh::{
    triangle_face::{PointId, TriangleFace},
    FaceId,
};
use super::ConvexHull3dError;
use crate::dim2::ConvexHull2d;

/// The initial convex hull structure built from the input points.
pub enum InitialConvexHull3d {
    Point(Vec<Vec3A>, Vec<[u32; 3]>),
    Segment(Vec<Vec3A>, Vec<[u32; 3]>),
    Triangle(Vec<Vec3A>, Vec<[u32; 3]>),
    Polyhedron(Vec<TriangleFace>),
}

/// An error returned during [`InitialConvexHull3d`] construction.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum InitialConvexHull3dError {
    #[error("Could not find initial segment points.")]
    MissingSegmentPoints,
    #[error("Could not find initial triangle point.")]
    MissingTrianglePoint,
}

/// Computes the covariance matrix of the given points.
fn cov(points: &[Vec3A]) -> SymmetricMat3 {
    // Compute the centroid.
    let centroid = points.iter().sum::<Vec3A>() / points.len() as f32;

    // Compute the covariance matrix of the points.
    // TODO: Should we multiply by the reciprocal inside the loop for better numerical stability?
    let mut cov = SymmetricMat3::ZERO;
    for point in points {
        cov += SymmetricMat3::from_outer_product(Vec3::from(*point - centroid));
    }
    cov /= points.len() as f32;
    cov
}

fn degenerate_point_hull(point: Vec3A) -> (Vec<Vec3A>, Vec<[u32; 3]>) {
    (vec![point], vec![[0; 3]; 2])
}

fn degenerate_segment_hull(direction: Vec3A, points: &[Vec3A]) -> (Vec<Vec3A>, Vec<[u32; 3]>) {
    // Find the maximum and minimum projections along the direction.
    let mut min_proj = f32::INFINITY;
    let mut max_proj = f32::NEG_INFINITY;
    let mut min_point = Vec3A::ZERO;
    let mut max_point = Vec3A::ZERO;

    for point in points {
        let proj = direction.dot(*point);
        if proj < min_proj {
            min_proj = proj;
            min_point = *point;
        }
        if proj > max_proj {
            max_proj = proj;
            max_point = *point;
        }
    }

    (vec![min_point, max_point], vec![[0, 1, 0], [1, 0, 0]])
}

fn planar_hull(
    points: &[Vec3A],
    normalized_points: &[Vec3A],
    eig: &SymmetricEigen3,
) -> InitialConvexHull3d {
    // The hull is a triangle.
    // Project the points onto the plane defined by the two largest eigenvectors.
    let u = Vec3A::from(eig.eigenvectors.y_axis.normalize());
    let v = Vec3A::from(eig.eigenvectors.z_axis.normalize());

    let mut subspace_points: Vec<Vec2> = Vec::with_capacity(normalized_points.len());
    for p in normalized_points {
        subspace_points.push(Vec2::new(p.dot(u), p.dot(v)));
    }

    // Compute the 2D convex hull of the projected points.
    let hull_2d_indices = ConvexHull2d::indices_from_points(&subspace_points);

    // Triangulate the 2D hull to form faces.
    let num_points = hull_2d_indices.len();
    let mut faces = Vec::with_capacity(2 * num_points - 4);

    for i in 1..(num_points - 1) {
        faces.push([0, i as u32, (i + 1) as u32]);
    }

    // Note: The bottom face uses a different starting point to avoid bad topology
    //       where an edge is shared by more than two faces.
    let end = num_points - 1;
    for i in 0..(end - 1) {
        faces.push([end as u32, i as u32 + 1, i as u32]);
    }

    InitialConvexHull3d::Triangle(hull_2d_indices.iter().map(|&i| points[i]).collect(), faces)
}

/// Computes the initial convex hull structure from the input points.
///
/// The initial hull is a tetrahedron formed by four points that are not coplanar, if such points exist.
/// If the points are coplanar, collinear, or coincident, the initial hull will be a triangle, line segment,
/// or point, respectively.
pub fn compute_initial_hull(
    points: &[Vec3A],
    normalized_points: &[Vec3A],
) -> Result<InitialConvexHull3d, ConvexHull3dError> {
    // Compute the eigen decomposition to see if the points are on a lower-dimensional subspace.
    let cov = cov(normalized_points);
    let eig = SymmetricEigen3::new(cov);

    // Count the number of non-zero eigenvalues to determine the rank.
    let rank = eig
        .eigenvalues
        .to_array()
        .iter()
        .filter(|v| relative_ne!(**v, 0.0, epsilon = 1e-6))
        .count();

    match rank {
        0 => {
            // The hull is a single point.
            let (vertices, faces) = degenerate_point_hull(points[0]);
            Ok(InitialConvexHull3d::Point(vertices, faces))
        }
        1 => {
            // The hull is a line segment.
            let direction = Vec3A::from(eig.eigenvectors.z_axis.normalize());
            let (vertices, faces) = degenerate_segment_hull(direction, points);
            Ok(InitialConvexHull3d::Segment(vertices, faces))
        }
        2 => Ok(planar_hull(points, normalized_points, &eig)),
        3 => {
            // The hull is a tetrahedron.

            // Find the four points that form the initial tetrahedron.
            let mut point_ids = [PointId::PLACEHOLDER; 4];

            let principal_axis = Vec3A::from(eig.eigenvectors.z_axis);

            // The first two vertices should be the ones farthest apart along the principal axis.
            let mut min_proj = f32::INFINITY;
            let mut max_proj = f32::NEG_INFINITY;
            for (i, point) in points.iter().enumerate() {
                let proj = principal_axis.dot(*point);
                if proj < min_proj {
                    min_proj = proj;
                    point_ids[0] = PointId(i as u32);
                }
                if proj > max_proj {
                    max_proj = proj;
                    point_ids[1] = PointId(i as u32);
                }
            }

            if point_ids[0] == PointId::PLACEHOLDER || point_ids[1] == PointId::PLACEHOLDER {
                return Err(ConvexHull3dError::InternalError(
                    "Could not find initial segment points.",
                ));
            }

            // The third vertex should be the one farthest from the line segment
            // between the first two vertices.
            let unit_01 = (points[point_ids[1].index()] - points[point_ids[0].index()]).normalize();

            let mut max_squared_distance = 0.0;

            for i in 0..points.len() {
                let diff = points[i] - points[point_ids[0].index()];
                let cross = unit_01.cross(diff);
                let distance_squared = cross.length_squared();

                if distance_squared > max_squared_distance
                    && points[i] != points[point_ids[0].index()]
                    && points[i] != points[point_ids[1].index()]
                {
                    max_squared_distance = distance_squared;
                    point_ids[2] = PointId(i as u32);
                }
            }

            if point_ids[2] == PointId::PLACEHOLDER {
                return Err(ConvexHull3dError::InternalError(
                    "Could not find initial triangle point.",
                ));
            }

            // The fourth vertex is the one furthest from the plane of the first three.
            let base = normalized_points[point_ids[0].index()];
            let plane_normal = (normalized_points[point_ids[1].index()] - base)
                .cross(normalized_points[point_ids[2].index()] - base)
                .normalize_or_zero();

            let mut furthest: f32 = 0.0;
            for (i, point) in normalized_points.iter().enumerate() {
                let distance = plane_normal.dot(*point - base);
                if distance.abs() > furthest.abs() {
                    furthest = distance;
                    point_ids[3] = PointId(i as u32);
                }
            }

            if point_ids[3] == PointId::PLACEHOLDER || furthest.abs() < 1e-7 {
                // The points are coplanar, so we return a planar hull instead.
                return Ok(planar_hull(points, normalized_points, &eig));
            }

            if furthest > 0.0 {
                point_ids.swap(0, 1);
            }

            let [a, b, c, d] = point_ids;
            let mut faces = vec![
                TriangleFace::from_triangle(normalized_points, [a, b, c]),
                TriangleFace::from_triangle(normalized_points, [a, d, b]),
                TriangleFace::from_triangle(normalized_points, [b, d, c]),
                TriangleFace::from_triangle(normalized_points, [c, d, a]),
            ];

            // Link the faces as neighbors.
            faces[0].set_neighbors(1, 2, 3, 2, 2, 2);
            faces[1].set_neighbors(3, 2, 0, 1, 0, 0);
            faces[2].set_neighbors(1, 3, 0, 1, 0, 1);
            faces[3].set_neighbors(2, 1, 0, 1, 0, 2);

            // Add outside points to the four faces.
            for (i, &point) in normalized_points.iter().enumerate() {
                if point == normalized_points[point_ids[0].index()]
                    || point == normalized_points[point_ids[1].index()]
                    || point == normalized_points[point_ids[2].index()]
                    || point == normalized_points[point_ids[3].index()]
                {
                    continue;
                }

                let point_id = PointId(i as u32);
                let mut furthest_face = FaceId::PLACEHOLDER;
                let mut max_distance = 0.0;

                for (j, face) in faces.iter().enumerate() {
                    if let Some(distance) =
                        face.distance_to_visible_point(point_id, normalized_points)
                    {
                        if distance > max_distance {
                            furthest_face = FaceId(j as u32);
                            max_distance = distance;
                        }
                    }
                }

                if furthest_face != FaceId::PLACEHOLDER {
                    faces[furthest_face.index()].try_add_outside_point(point_id, normalized_points);
                }
            }

            // TODO: Make this optional behind a feature flag or something.
            #[cfg(debug_assertions)]
            {
                for id in 0..faces.len() as u32 {
                    validate_face_connectivity(FaceId(id), &faces);
                }
            }

            Ok(InitialConvexHull3d::Polyhedron(faces))
        }
        _ => unreachable!("Rank can only be 0, 1, 2, or 3."),
    }
}
