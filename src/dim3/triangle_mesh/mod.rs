pub(super) mod triangle_face;
// TODO: Make validation methods public, and move auto-validation behind a feature flag.
#[cfg(debug_assertions)]
pub(super) mod validation;

use crate::dim3::{
    initial_hull::{compute_initial_hull, InitialConvexHull3d},
    normalize, ConvexHull3dError,
};
pub(super) use triangle_face::{FaceHandle, FaceId, PointId, TriangleFace};

use glam::Vec3A;

/// A 3D [convex hull] with triangular faces, representing the smallest convex set
/// containing all input points in a given point set.
///
/// This can be thought of as a shrink wrapping of a 3D object.
///
/// For a more advanced representation that supports polygonal faces, see [`ConvexHull3d`].
///
/// [convex hull]: https://en.wikipedia.org/wiki/Convex_hull
/// [`ConvexHull3d`]: super::hull::ConvexHull3d
///
/// # Example
///
/// ```
/// use glam::Vec3A;
/// use quickhull::ConvexTriangleMesh;
///
/// // Define a set of 3D points.
/// let points = vec![
///     Vec3A::new(0.0, 0.0, 0.0),
///     Vec3A::new(1.0, 0.0, 0.0),
///     Vec3A::new(0.0, 1.0, 0.0),
///     Vec3A::new(0.0, 0.0, 1.0),
///     // Additional internal points that should not affect the hull
///     Vec3A::new(0.1, 0.1, 0.1),
///     Vec3A::new(0.25, 0.25, 0.25),
/// ];
///
/// // No limit on the number of iterations.
/// let max_iter = None;
///
/// // Compute the convex hull.
/// let hull = ConvexTriangleMesh::try_from_points(&points, max_iter).unwrap();
///
/// // Get the vertices and indices of the convex hull.
/// let (vertices, indices) = hull.into_parts();
///
/// // The hull should be a tetrahedron with 4 vertices and 4 triangular faces.
/// assert_eq!(vertices.len(), 4);
/// assert_eq!(indices.len(), 4);
/// ```
#[derive(Clone, Debug, Default)]
pub struct ConvexTriangleMesh {
    /// The vertices of the convex hull.
    vertices: Vec<Vec3A>,

    /// The triangle faces of the convex hull.
    indices: Vec<[u32; 3]>,
}

impl ConvexTriangleMesh {
    /// Attempts to compute a [`ConvexTriangleMesh`] for the given set of points.
    ///
    /// `max_iter` specifies the maximum number of iterations to perform.
    /// If `None`, the algorithm will run until completion.
    ///
    /// Point sets with fewer than 4 points will produce degenerate hulls
    /// representing a point, line segment, or triangle. If this is not desired,
    /// check for the number of input points before constructing the hull.
    ///
    /// # Errors
    ///
    /// Returns a [`ConvexHull3dError`] if hull construction fails.
    pub fn try_from_points(
        points: &[Vec3A],
        max_iter: Option<usize>,
    ) -> Result<Self, ConvexHull3dError> {
        if points.is_empty() {
            // Empty hull.
            return Ok(ConvexTriangleMesh {
                vertices: Vec::new(),
                indices: Vec::new(),
            });
        }

        // Construct a normalized point cloud for better numerical stability.
        let mut normalized_points = points.to_vec();
        normalize::normalize_point_cloud(&mut normalized_points);

        let mut faces: Vec<TriangleFace>;
        let mut undecided_points: Vec<PointId> = Vec::new();

        // Create the initial simplex.
        match compute_initial_hull(points, &normalized_points)? {
            InitialConvexHull3d::Point(vertices, indices)
            | InitialConvexHull3d::Segment(vertices, indices)
            | InitialConvexHull3d::Triangle(vertices, indices) => {
                return Ok(ConvexTriangleMesh { vertices, indices });
            }
            InitialConvexHull3d::Polyhedron(initial_faces) => {
                faces = initial_faces;
            }
        }

        // Run the main quickhull algorithm.
        Self::update(
            &normalized_points,
            &mut undecided_points,
            &mut faces,
            max_iter,
        )?;

        // Collect valid face indices.
        let mut indices: Vec<[u32; 3]> = Vec::new();
        for face in faces.iter() {
            if face.valid {
                indices.push(face.points.map(|id| id.0));
            }
        }

        // Shrink the hull, removing unused points.
        let mut vertices: Vec<Vec3A> = points.to_vec();
        Self::remove_unused_points(&mut vertices, &mut indices);

        Ok(ConvexTriangleMesh { vertices, indices })
    }

    /// Creates a [`ConvexTriangleMesh`] from the given vertices and indices without performing any checks.
    ///
    /// The caller should ensure that the vertices and indices form a valid convex triangle mesh.
    #[inline]
    pub const fn from_vertices_indices_unchecked(
        points: Vec<Vec3A>,
        indices: Vec<[u32; 3]>,
    ) -> Self {
        ConvexTriangleMesh {
            vertices: points,
            indices,
        }
    }

    /// Returns the vertices of the convex hull.
    #[inline]
    pub fn vertices(&self) -> &[Vec3A] {
        &self.vertices
    }

    /// Returns the indices of the convex hull's faces.
    #[inline]
    pub fn indices(&self) -> &[[u32; 3]] {
        &self.indices
    }

    /// Returns the vertices and indices of the convex hull.
    ///
    /// This consumes the hull and allows taking ownership of the underlying data without cloning.
    #[inline]
    pub fn into_parts(self) -> (Vec<Vec3A>, Vec<[u32; 3]>) {
        (self.vertices, self.indices)
    }

    /// The main quickhull algorithm.
    fn update(
        points: &[Vec3A],
        undecided_points: &mut Vec<PointId>,
        faces: &mut Vec<TriangleFace>,
        max_iter: Option<usize>,
    ) -> Result<(), ConvexHull3dError> {
        let mut horizon: Vec<FaceHandle> = Vec::new();
        let mut horizon_fixing_workspace: Vec<u32> = Vec::with_capacity(points.len());
        let mut removed_faces: Vec<FaceId> = Vec::new();

        let max_iter = max_iter.unwrap_or(usize::MAX);

        // The main algorithm of quickhull.
        //
        // For each face that has outside points:
        //
        // 1. Find the outside point that is farthest from the face, the "eye point".
        // 2. Find the "horizon", the vertices that form the boundary between the visible
        //    and non-visible parts of the current hull from the viewpoint of the eye point.
        // 3. Create new faces connecting the eye point to the horizon.
        // 4. Reassign outside points from removed faces to the new faces.
        //
        // Repeat until no faces with outside points remain or the maximum number of iterations is reached.
        let mut i = 0;
        while i != faces.len() {
            horizon.clear();

            if i >= max_iter {
                break;
            }

            let face = &mut faces[i];

            if !face.valid || face.affinely_dependent {
                i += 1;
                continue;
            }

            // Select the furthest point.
            let Some((furthest_point_index, _)) = face.furthest_outside_point else {
                i += 1;
                continue;
            };

            // Mark the face as removed for now.
            face.valid = false;
            removed_faces.clear();
            removed_faces.push(FaceId(i as u32));

            // Compute the horizon.
            for j in 0..3 {
                let face = &faces[i];
                compute_horizon(
                    face.neighbors[j],
                    furthest_point_index,
                    points,
                    faces,
                    &mut removed_faces,
                    &mut horizon,
                );
            }

            // Due to round-off errors, the horizon may contain self-intersections
            // or multiple disjoint loops. Fix the horizon topology if needed.
            fix_horizon_topology(
                points,
                faces,
                &mut removed_faces,
                &mut horizon,
                &mut horizon_fixing_workspace,
            )?;

            if horizon.is_empty() {
                // Due to round-off errors, the horizon could not be computed.
                let is_any_valid = faces[i + 1..]
                    .iter()
                    .any(|f| f.valid && !f.affinely_dependent);

                if is_any_valid {
                    return Err(ConvexHull3dError::InternalError(
                        "Could not compute horizon.",
                    ));
                }

                // All remaining faces are invalid, so we are done.
                faces[i].valid = true;
                break;
            }

            // Create and attach new faces.
            create_and_attach_faces(
                furthest_point_index,
                points,
                faces,
                &mut removed_faces,
                undecided_points,
                &horizon,
            );

            i += 1;
        }

        Ok(())
    }

    /// Removes unused points from the given vertex buffer and remaps indices accordingly.
    fn remove_unused_points(points: &mut Vec<Vec3A>, faces: &mut [[u32; 3]]) {
        let mut used: Vec<bool> = vec![false; points.len()];
        let mut remap: Vec<usize> = (0..points.len()).collect();

        for i in faces.iter() {
            used[i[0] as usize] = true;
            used[i[1] as usize] = true;
            used[i[2] as usize] = true;
        }

        let mut i = 0;
        while i < points.len() {
            if !used[i] {
                points.swap_remove(i);
                remap[points.len()] = i;
                used[i] = used[points.len()];
            } else {
                i += 1;
            }
        }

        // Remap face indices.
        for i in faces.iter_mut() {
            i[0] = remap[i[0] as usize] as u32;
            i[1] = remap[i[1] as usize] as u32;
            i[2] = remap[i[2] as usize] as u32;
        }
    }

    /// Computes the volume of the convex hull.
    #[inline]
    pub fn volume(&self) -> f32 {
        self.indices
            .iter()
            .map(|triangle| {
                let p0 = self.vertices[triangle[0] as usize];
                let p1 = self.vertices[triangle[1] as usize];
                let p2 = self.vertices[triangle[2] as usize];

                // Volume of the tetrahedron formed by the triangle and the origin.
                (p0.dot(p1.cross(p2))).abs() / 6.0
            })
            .sum()
    }

    /// Computes the point on the convex hull that is furthest in the given direction.
    #[inline]
    pub fn support_point(&self, direction: Vec3A) -> Vec3A {
        let mut max = self.vertices[0].dot(direction);
        let mut index = 0;

        for (i, point) in self.vertices.iter().enumerate().skip(1) {
            let dot_product = point.dot(direction);
            if dot_product > max {
                max = dot_product;
                index = i;
            }
        }

        self.vertices[index]
    }
}

/// Computes the horizon of the convex hull from the viewpoint of the given point.
///
/// The horizon is represented as a list of ridges (edges) and the unvisible face opposite to each ridge.
fn compute_horizon(
    start_face: FaceHandle,
    point: PointId,
    points: &[Vec3A],
    faces: &mut [TriangleFace],
    removed_faces: &mut Vec<FaceId>,
    horizon: &mut Vec<FaceHandle>,
) {
    // Maintain a DFS stack of faces to visit.
    // Note that this could also be implemented using recursion.
    let mut stack: Vec<FaceHandle> = Vec::with_capacity(32);
    stack.push(start_face);

    while let Some(handle) = stack.pop() {
        let face = &mut faces[handle.face.index()];

        // Skip already removed faces.
        if !face.valid {
            continue;
        }

        if !face.can_be_seen_by_point_order_independent(point, points) {
            // Face is not visible, so it's part of the horizon.
            horizon.push(handle);
        } else {
            // Face is visible, so remove it.
            face.valid = false;
            removed_faces.push(handle.face);

            // Push neighbors to stack.
            // Note the order: we want to visit neighbor1 before neighbor2,
            let neighbor2 = face.neighbors[(handle.edge.index() + 2) % 3];
            let neighbor1 = face.neighbors[(handle.edge.index() + 1) % 3];
            stack.push(neighbor2);
            stack.push(neighbor1);
        }
    }
}

/// Fixes the topology of the horizon if it contains self-intersections or multiple loops.
fn fix_horizon_topology(
    points: &[Vec3A],
    faces: &mut [TriangleFace],
    removed_faces: &mut Vec<FaceId>,
    horizon: &mut Vec<FaceHandle>,
    // TODO: We can probably use u16 or even u8 here.
    workspace: &mut Vec<u32>,
) -> Result<(), ConvexHull3dError> {
    // Clear and resize the workspace.
    workspace.clear();
    workspace.resize(points.len(), 0);

    let mut needs_fixing = false;

    for handle in horizon.iter() {
        let point = faces[handle.face.index()].second_point_from_edge(handle.edge);
        let workspace_value = &mut workspace[point.index()];
        *workspace_value += 1;

        // If any edge is used more than once, we need to fix the horizon.
        if *workspace_value > 1 {
            needs_fixing = true;
        }
    }

    if needs_fixing {
        // The horizon has multiple loops or self-intersections.

        // First, find which loop we want to keep.
        let mut loop_start = 0;
        for &FaceHandle { face, edge } in horizon.iter() {
            let point1 = points[faces[face.index()].second_point_from_edge(edge).index()];
            let point2 = points[faces[face.index()].first_point_from_edge(edge).index()];
            let direction = point2 - point1;

            let support_index = support_point_iterator_position(
                direction,
                points,
                horizon
                    .iter()
                    .map(|h| faces[h.face.index()].second_point_from_edge(h.edge).index()),
            )
            .ok_or(ConvexHull3dError::MissingSupportPoint)?;
            let selected = horizon[support_index];

            if workspace[faces[selected.face.index()]
                .second_point_from_edge(selected.edge)
                .index()]
                == 1
            {
                // This edge is only used once, so we can start from here.
                loop_start = support_index;
                break;
            }
        }

        // Now, reconstruct the horizon using only the selected loop.
        let mut removing_index = None;
        let old_horizon = core::mem::take(horizon);

        for i in 0..old_horizon.len() {
            let face_index = (loop_start + i) % old_horizon.len();
            let FaceHandle { face, edge } = old_horizon[face_index];

            match removing_index {
                Some(idx) => {
                    let point = faces[face.index()].second_point_from_edge(edge);
                    if idx == point.index() {
                        removing_index = None;
                    }
                }
                _ => {
                    let point = faces[face.index()].second_point_from_edge(edge);
                    if workspace[point.index()] > 1 {
                        removing_index = Some(point.index());
                    }
                }
            }

            if removing_index.is_some() {
                if faces[face.index()].valid {
                    faces[face.index()].valid = false;
                    removed_faces.push(face);
                }
            } else {
                horizon.push(FaceHandle::new(face, edge));
            }
        }
    }

    Ok(())
}

fn create_and_attach_faces(
    point: PointId,
    points: &[Vec3A],
    faces: &mut Vec<TriangleFace>,
    removed_faces: &mut [FaceId],
    undecided_points: &mut Vec<PointId>,
    horizon: &[FaceHandle],
) {
    // We will add the new faces directly to the faces list, so we need to know the current length.
    let index_offset = faces.len();

    // Reserve space for new faces.
    faces.reserve(horizon.len());

    // Create new faces connecting the horizon to the point.
    for FaceHandle { face, edge } in horizon.iter().copied() {
        faces.push(TriangleFace::from_triangle(
            points,
            [
                point,
                faces[face.index()].second_point_from_edge(edge),
                faces[face.index()].first_point_from_edge(edge),
            ],
        ));
    }

    // Link the new faces to their neighbors.
    for i in 0..horizon.len() {
        let previous_face = if i == 0 {
            index_offset as u32 + (horizon.len() as u32) - 1
        } else {
            index_offset as u32 + (i as u32) - 1
        };

        let middle = horizon[i];
        let next_face = index_offset as u32 + (i as u32 + 1) % (horizon.len() as u32);

        // Link the new face to its neighbors.
        let face = &mut faces[index_offset + i];
        face.set_neighbors(previous_face, middle.face, next_face, 2, middle.edge, 0);
        debug_assert!(
            !faces[faces[middle.face.index()]
                .neighbor(middle.edge)
                .face
                .index()]
            .valid,
            "tried to overwrite a valid face link"
        );

        // Link the middle face to the new face.
        let middle_face_neighbor = faces[middle.face.index()].neighbor_mut(middle.edge);
        middle_face_neighbor.face.0 = (index_offset + i) as u32;
        middle_face_neighbor.edge.0 = 1;
    }

    // Reassign outside points from removed faces to the new faces.
    for face in removed_faces.iter() {
        // TODO: Avoid cloning the outside points.
        let outside_points = faces[face.index()].outside_points.clone();
        for outside_point_index in outside_points {
            if points[outside_point_index.index()] == points[point.index()] {
                continue;
            }

            let mut best_distance = f32::MIN;
            let mut best_face = None;

            for (i, face) in faces.iter().enumerate().skip(index_offset) {
                if face.affinely_dependent {
                    continue;
                }

                let distance = face.distance_to_point(outside_point_index, points);
                if distance > best_distance {
                    best_distance = distance;
                    best_face = Some(FaceId(i as u32));
                }
            }

            if let Some(best_face) = best_face {
                let best_face = &mut faces[best_face.index()];
                best_face.try_add_outside_point(outside_point_index, points);
            }

            // If none of the new faces can see the point, it is implicitly removed,
            // because it won't be assigned to any face.
        }
    }

    // Try to assign collinear points in `undecided_points` to one of the new faces.
    let mut i = 0;
    while i != undecided_points.len() {
        let mut best_distance = f32::MIN;
        let mut best_face_index = None;
        let undecided_point_index = undecided_points[i];

        for (j, face) in faces.iter().enumerate().skip(index_offset) {
            if let Some(distance) = face.distance_to_visible_point(undecided_point_index, points) {
                if distance > best_distance {
                    best_distance = distance;
                    best_face_index = Some(j as u32);
                }
            }
        }

        if let Some(best_face_index) = best_face_index {
            let best_face = &mut faces[best_face_index as usize];
            best_face.add_outside_point(undecided_point_index, points);
            undecided_points.swap_remove(i);
        } else {
            i += 1;
        }
    }
}

/// Finds the support point index in the given iterator of point indices.
fn support_point_iterator_position<I>(direction: Vec3A, points: &[Vec3A], iter: I) -> Option<usize>
where
    I: Iterator<Item = usize>,
{
    let mut max_dot = f32::MIN;
    let mut max_index = None;

    for (k, i) in iter.enumerate() {
        let point = points[i];
        let dot = point.dot(direction);
        if dot > max_dot {
            max_dot = dot;
            max_index = Some(k);
        }
    }

    max_index
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_point_set() {
        let points: Vec<Vec3A> = Vec::new();
        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for empty point set");
        let (vertices, indices) = result.into_parts();
        assert!(vertices.is_empty());
        assert!(indices.is_empty());
    }

    #[test]
    fn single_point() {
        let points = vec![Vec3A::splat(1.0)];
        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for single point");
        let (vertices, indices) = result.into_parts();
        assert_eq!(vertices, vec![Vec3A::splat(1.0)]);
        assert_eq!(indices, vec![[0; 3]; 2]);
    }

    #[test]
    fn two_points() {
        let points = vec![
            Vec3A::splat(1.0),
            Vec3A::splat(2.0),
            Vec3A::splat(1.0),
            Vec3A::splat(2.0),
        ];
        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for two points");
        let (vertices, indices) = result.into_parts();
        assert_eq!(vertices, vec![Vec3A::splat(1.0), Vec3A::splat(2.0)]);
        assert_eq!(indices, vec![[0, 1, 0], [1, 0, 0]]);
    }

    #[test]
    fn three_points() {
        let points = vec![
            Vec3A::new(0.0, 0.0, 0.0),
            Vec3A::new(1.0, 0.0, 0.0),
            Vec3A::new(0.0, 1.0, 0.0),
        ];
        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for three points");
        let (vertices, indices) = result.into_parts();
        assert_eq!(vertices, points);
        assert_eq!(indices, vec![[0, 1, 2], [2, 1, 0]]);
    }

    #[test]
    fn four_points_coincident() {
        let points = (0..4).map(|_| Vec3A::splat(1.0)).collect::<Vec<_>>();

        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for coincident points");
        let (vertices, indices) = result.into_parts();
        assert_eq!(vertices, vec![Vec3A::splat(1.0)]);
        assert_eq!(indices, vec![[0; 3]; 2]);
    }

    #[test]
    fn four_points_collinear() {
        let mut points = (0..4).map(|_| Vec3A::splat(1.0)).collect::<Vec<_>>();
        points[0].x += f32::EPSILON;
        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for collinear points");
        let (vertices, indices) = result.into_parts();
        assert_eq!(vertices, vec![points[1], points[0]]);
        assert_eq!(indices, vec![[0, 1, 0], [1, 0, 0]]);
    }

    #[test]
    fn four_points_coplanar() {
        let mut points = (0..4).map(|_| Vec3A::splat(1.0)).collect::<Vec<_>>();
        points[0].x += f32::EPSILON;
        points[1].y += f32::EPSILON;
        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for coplanar points");
        let (vertices, indices) = result.into_parts();
        assert_eq!(vertices, vec![points[2], points[0], points[1]]);
        assert_eq!(indices, vec![[0, 1, 2], [2, 1, 0]]);
    }

    #[test]
    fn four_points_min_volume() {
        let mut points = (0..4).map(|_| Vec3A::splat(1.0)).collect::<Vec<_>>();
        points[0].x += 3.0 * f32::EPSILON;
        points[1].y += 3.0 * f32::EPSILON;
        points[2].z += 3.0 * f32::EPSILON;
        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for tetrahedron");
        let (vertices, indices) = result.into_parts();

        // Four points that do not lie in a plane have a tetrahedron for a hull.
        assert_eq!(vertices.len(), 4);
        assert_eq!(indices.len(), 4);
    }

    #[test]
    fn volume_should_be_positive() {
        let mut points = (0..4).map(|_| Vec3A::splat(1.0)).collect::<Vec<_>>();
        points[0].x += 1.0 * f32::EPSILON;
        points[1].y += 1.0 * f32::EPSILON;
        points[2].z += 2.0 * f32::EPSILON;
        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for tetrahedron");
        assert!(result.volume() > 0.0);
    }

    #[test]
    fn inner_outer_test() {
        let p1 = Vec3A::new(1.0, 0.0, 0.0);
        let p2 = Vec3A::new(0.0, 1.0, 0.0);
        let p3 = Vec3A::new(0.0, 0.0, 1.0);
        let outer_point = Vec3A::new(0.0, 0.0, 10.0);
        let inner_point = Vec3A::new(0.0, 0.0, 0.0);
        let whithin_point = Vec3A::new(1.0, 0.0, 0.0);
        let points = vec![p1, p2, p3, outer_point, inner_point, whithin_point];
        let face = TriangleFace::from_triangle(&points, [PointId(0), PointId(1), PointId(2)]);
        let can_see_outer = face
            .distance_to_visible_point(PointId(3), &points)
            .is_some();
        assert!(can_see_outer);
        let can_see_inner = face
            .distance_to_visible_point(PointId(4), &points)
            .is_some();
        assert!(!can_see_inner);
    }

    #[test]
    fn octahedron_test() {
        let p1 = Vec3A::new(1.0, 0.0, 0.0);
        let p2 = Vec3A::new(0.0, 1.0, 0.0);
        let p3 = Vec3A::new(0.0, 0.0, 1.0);
        let p4 = Vec3A::new(-1.0, 0.0, 0.0);
        let p5 = Vec3A::new(0.0, -1.0, 0.0);
        let p6 = Vec3A::new(0.0, 0.0, -1.0);
        let (_v, i) = ConvexTriangleMesh::try_from_points(&[p1, p2, p3, p4, p5, p6], None)
            .unwrap()
            .into_parts();
        assert_eq!(i.len(), 8);
    }

    #[test]
    fn octahedron_translation_test() {
        let p1 = Vec3A::new(1.0, 0.0, 0.0);
        let p2 = Vec3A::new(0.0, 1.0, 0.0);
        let p3 = Vec3A::new(0.0, 0.0, 1.0);
        let p4 = Vec3A::new(-1.0, 0.0, 0.0);
        let p5 = Vec3A::new(0.0, -1.0, 0.0);
        let p6 = Vec3A::new(0.0, 0.0, -1.0);
        let points: Vec<_> = [p1, p2, p3, p4, p5, p6]
            .into_iter()
            .map(|p| p + Vec3A::splat(10.0))
            .collect();
        let (_v, i) = ConvexTriangleMesh::try_from_points(&points, None)
            .unwrap()
            .into_parts();
        assert_eq!(i.len(), 8);
    }

    #[test]
    fn cube_test() {
        let p1 = Vec3A::new(1.0, 1.0, 1.0);
        let p2 = Vec3A::new(1.0, 1.0, -1.0);
        let p3 = Vec3A::new(1.0, -1.0, 1.0);
        let p4 = Vec3A::new(1.0, -1.0, -1.0);
        let p5 = Vec3A::new(-1.0, 1.0, 1.0);
        let p6 = Vec3A::new(-1.0, 1.0, -1.0);
        let p7 = Vec3A::new(-1.0, -1.0, 1.0);
        let p8 = Vec3A::new(-1.0, -1.0, -1.0);
        let (_v, i) = ConvexTriangleMesh::try_from_points(&[p1, p2, p3, p4, p5, p6, p7, p8], None)
            .unwrap()
            .into_parts();
        assert_eq!(i.len(), 6 * 2);
    }

    #[test]
    fn cube_volume_test() {
        let p1 = Vec3A::new(2.0, 2.0, 2.0);
        let p2 = Vec3A::new(2.0, 2.0, 0.0);
        let p3 = Vec3A::new(2.0, 0.0, 2.0);
        let p4 = Vec3A::new(2.0, 0.0, 0.0);
        let p5 = Vec3A::new(0.0, 2.0, 2.0);
        let p6 = Vec3A::new(0.0, 2.0, 0.0);
        let p7 = Vec3A::new(0.0, 0.0, 2.0);
        let p8 = Vec3A::new(0.0, 0.0, 0.0);
        let cube =
            ConvexTriangleMesh::try_from_points(&[p1, p2, p3, p4, p5, p6, p7, p8], None).unwrap();
        assert_eq!(cube.volume(), 8.0);
    }

    // Heavy test (~ 0.75s)
    #[test]
    fn sphere_volume_test() {
        let points = sphere_points(50);
        let hull = ConvexTriangleMesh::try_from_points(&points, None).unwrap();
        let volume = hull.volume();
        let expected_volume = 4.0 / 3.0 * std::f32::consts::PI;
        assert!(
            (volume - expected_volume).abs() < 0.1,
            "Expected {expected_volume}, got {volume}"
        );
    }

    #[test]
    fn cube_support_point_test() {
        let p1 = Vec3A::new(1.0, 1.0, 1.0);
        let p2 = Vec3A::new(1.0, 1.0, 0.0);
        let p3 = Vec3A::new(1.0, 0.0, 1.0);
        let p4 = Vec3A::new(1.0, 0.0, 0.0);
        let p5 = Vec3A::new(0.0, 1.0, 1.0);
        let p6 = Vec3A::new(0.0, 1.0, 0.0);
        let p7 = Vec3A::new(0.0, 0.0, 1.0);
        let p8 = Vec3A::new(0.0, 0.0, 0.0);
        let cube =
            ConvexTriangleMesh::try_from_points(&[p1, p2, p3, p4, p5, p6, p7, p8], None).unwrap();
        assert_eq!(cube.support_point(Vec3A::splat(0.5)), p1);
    }

    #[test]
    fn flat_test() {
        let p1 = Vec3A::new(1.0, 1.0, 10.0);
        let p2 = Vec3A::new(1.0, 1.0, 10.0);
        let p3 = Vec3A::new(1.0, -1.0, 10.0);
        let p4 = Vec3A::new(1.0, -1.0, 10.0);
        let p5 = Vec3A::new(-1.0, 1.0, 10.0);
        let p6 = Vec3A::new(-1.0, 1.0, 10.0);
        let p7 = Vec3A::new(-1.0, -1.0, 10.0);
        let p8 = Vec3A::new(-1.0, -1.0, 10.0);

        let result = ConvexTriangleMesh::try_from_points(&[p1, p2, p3, p4, p5, p6, p7, p8], None)
            .expect("could not compute hull for flat points");
        let (vertices, indices) = result.into_parts();

        assert_eq!(vertices, vec![p3, p1, p5, p7]);
        assert_eq!(indices, vec![[0, 1, 2], [0, 2, 3], [3, 1, 0], [3, 2, 1]]);
    }

    #[test]
    fn line_test() {
        let points = (0..10)
            .map(|i| Vec3A::new(i as f32, 1.0, 10.0))
            .collect::<Vec<_>>();

        let result = ConvexTriangleMesh::try_from_points(&points, None)
            .expect("could not compute hull for line points");
        let (vertices, indices) = result.into_parts();

        assert_eq!(vertices, vec![points[0], points[9]]);
        assert_eq!(indices, vec![[0, 1, 0], [1, 0, 0]]);
    }

    #[test]
    fn simplex_may_degenerate_test() {
        let points = vec![
            Vec3A::new(1.0, 0.0, 1.0),
            Vec3A::new(1.0, 1.0, 1.0),
            Vec3A::new(2.0, 1.0, 0.0),
            Vec3A::new(2.0, 1.0, 1.0),
            Vec3A::new(2.0, 0.0, 1.0),
            Vec3A::new(2.0, 0.0, 0.0),
            Vec3A::new(1.0, 1.0, 2.0),
            Vec3A::new(0.0, 1.0, 2.0),
            Vec3A::new(0.0, 0.0, 2.0),
            Vec3A::new(1.0, 0.0, 2.0),
        ];
        let (_v, _i) = ConvexTriangleMesh::try_from_points(&points, None)
            .unwrap()
            .into_parts();
    }

    #[test]
    fn simplex_may_degenerate_test_2() {
        let vertices = vec![
            Vec3A::new(0., 0., 0.),
            Vec3A::new(1., 0., 0.),
            Vec3A::new(1., 0., 1.),
            Vec3A::new(0., 0., 1.),
            Vec3A::new(0., 1., 0.),
            Vec3A::new(1., 1., 0.),
            Vec3A::new(1., 1., 1.),
            Vec3A::new(0., 1., 1.),
            Vec3A::new(2., 1., 0.),
            Vec3A::new(2., 1., 1.),
            Vec3A::new(2., 0., 1.),
            Vec3A::new(2., 0., 0.),
            Vec3A::new(1., 1., 2.),
            Vec3A::new(0., 1., 2.),
            Vec3A::new(0., 0., 2.),
            Vec3A::new(1., 0., 2.),
        ];
        let indices = [4, 5, 1, 11, 1, 5, 1, 11, 10, 10, 2, 1, 5, 8, 11];
        let points = indices.iter().map(|i| vertices[*i]).collect::<Vec<_>>();
        let (_v, _i) = ConvexTriangleMesh::try_from_points(&points, None)
            .unwrap()
            .into_parts();
    }

    #[cfg(test)]
    fn sphere_points(divisions: usize) -> Vec<Vec3A> {
        fn rot_z(point: Vec3A, angle: f32) -> Vec3A {
            let e1 = angle.cos() * point[0] - angle.sin() * point[1];
            let e2 = angle.sin() * point[0] + angle.cos() * point[1];
            let e3 = point[2];
            Vec3A::new(e1, e2, e3)
        }
        fn rot_x(point: Vec3A, angle: f32) -> Vec3A {
            let e1 = point[0];
            let e2 = angle.cos() * point[1] - angle.sin() * point[2];
            let e3 = angle.sin() * point[1] + angle.cos() * point[2];
            Vec3A::new(e1, e2, e3)
        }
        let mut points = Vec::new();
        let unit_y = Vec3A::Y;
        for step_x in 0..divisions {
            let angle_x = 2.0 * std::f32::consts::PI * (step_x as f32 / divisions as f32);
            let p = rot_x(unit_y, angle_x);
            for step_z in 0..divisions {
                let angle_z = 2.0 * std::f32::consts::PI * (step_z as f32 / divisions as f32);
                let p = rot_z(p, angle_z);
                points.push(p);
            }
        }
        points
    }

    #[test]
    fn sphere_test() {
        let points = sphere_points(10);
        let (_v, _i) = ConvexTriangleMesh::try_from_points(&points, None)
            .unwrap()
            .into_parts();
    }

    /// Useful for fuzzing and profiling.
    /// Creates a sea-urchin like point cloud with points distributed arbitrarily within a sphere.
    #[test]
    fn heavy_sea_urchin_test() {
        use rand::prelude::{Distribution, SeedableRng, SliceRandom};

        // increase this to ~1000 to gather more samples for a sampling profiler
        let iterations = 1;

        for s in 0..iterations {
            let mut rng = rand::rngs::StdRng::seed_from_u64(s);
            let dist = rand::distr::StandardUniform;

            fn rot_z(point: Vec3A, angle: f32) -> Vec3A {
                let e1 = angle.cos() * point[0] - angle.sin() * point[1];
                let e2 = angle.sin() * point[0] + angle.cos() * point[1];
                let e3 = point[2];
                Vec3A::new(e1, e2, e3)
            }
            fn rot_x(point: Vec3A, angle: f32) -> Vec3A {
                let e1 = point[0];
                let e2 = angle.cos() * point[1] - angle.sin() * point[2];
                let e3 = angle.sin() * point[1] + angle.cos() * point[2];
                Vec3A::new(e1, e2, e3)
            }
            let mut points = Vec::new();
            let dev = 100;
            let unit_y = Vec3A::Y;
            for step_x in 0..dev {
                let angle_x = 2.0 * std::f32::consts::PI * (step_x as f32 / dev as f32);
                let p = rot_x(unit_y, angle_x);
                for step_z in 0..dev {
                    let angle_z = 2.0 * std::f32::consts::PI * (step_z as f32 / dev as f32);
                    let p = rot_z(p, angle_z);
                    let rand_offset: f32 = dist.sample(&mut rng);
                    points.push(p * rand_offset);
                }
            }

            points.shuffle(&mut rng);
            let (_v, _i) = ConvexTriangleMesh::try_from_points(&points, None)
                .unwrap()
                .into_parts();
        }
    }
}
