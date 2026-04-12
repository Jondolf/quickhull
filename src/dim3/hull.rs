use glam::Vec3A;

use super::plane::Plane3d;
use crate::{collections::HashMap, ConvexHull3dError, ConvexTriangleMesh};

/// A 3D [convex hull] representing the smallest convex set containing
/// all input points in a given point set.
///
/// This can be thought of as a shrink wrapping of a 3D object.
///
/// Unlike [`ConvexTriangleMesh`], which always has triangular faces,
/// [`ConvexHull3d`] supports polygonal faces formed by merging coplanar triangles.
///
/// The maximum number of vertices is 65,535 ([`u16::MAX`]).
///
/// [convex hull]: https://en.wikipedia.org/wiki/Convex_hull
///
/// # Example
///
/// ```
/// use glam::Vec3A;
/// use quickhull::{ConvexHull3d, ConvexHull3dSettings};
///
/// // Define points representing the corners of a cube.
/// let points = vec![
///     Vec3A::new( 1.0,  1.0,  1.0),
///     Vec3A::new( 1.0,  1.0, -1.0),
///     Vec3A::new( 1.0, -1.0,  1.0),
///     Vec3A::new( 1.0, -1.0, -1.0),
///     Vec3A::new(-1.0,  1.0,  1.0),
///     Vec3A::new(-1.0,  1.0, -1.0),
///     Vec3A::new(-1.0, -1.0,  1.0),
///     Vec3A::new(-1.0, -1.0, -1.0),
/// ];
///
/// let settings = ConvexHull3dSettings {
///     // Merge faces with normals within 1e-4 radians of each other
///     coplanarity_tolerance: 1e-4,
///     // Allow up to 10,000 iterations for hull construction (should not be hit in practice)
///     max_iterations: 10_000,
/// };
///
/// // Compute the convex hull.
/// let hull = ConvexHull3d::try_from_points(&points, settings).unwrap();
///
/// // The cube's hull should have 8 vertices and 6 quad faces.
/// assert_eq!(hull.vertices().len(), 8);
/// assert_eq!(hull.faces().len(), 6);
/// ```
#[derive(Clone, Debug, Default)]
pub struct ConvexHull3d {
    /// The vertices of the convex hull.
    vertices: Vec<Vec3A>,

    /// The indices of the convex hull vertices that define the faces of the hull.
    ///
    /// Each face is defined by a contiguous sequence of vertex indices, with the starting index
    /// and number of vertices for each face specified in `faces`.
    indices: Vec<u16>,

    /// The faces of the convex hull, defined by the indices in `indices`.
    faces: Vec<HullFace>,

    /// The bounding planes associated with each face of the convex hull.
    ///
    /// These are one-to-one with `faces`.
    planes: Vec<Plane3d>,
}

/// The index of the first vertex of the face in [`ConvexHull3d::indices`]
/// and the number of vertices in the face. Used to define the faces of the convex hull.
#[derive(Clone, Debug, PartialEq)]
pub struct HullFace {
    /// The index of the first vertex of the face in [`ConvexHull3d::indices`].
    first_vertex: u16,

    /// The number of vertices in the face.
    num_vertices: u16,
}

impl HullFace {
    /// Returns the index of the first vertex of the face in [`ConvexHull3d::indices`].
    #[inline]
    pub fn first_vertex(&self) -> u16 {
        self.first_vertex
    }

    /// Returns the number of vertices in the face.
    #[inline]
    pub fn num_vertices(&self) -> u16 {
        self.num_vertices
    }

    /// Returns the indices of the vertices that make up this face.
    #[inline]
    pub fn vertex_indices<'a>(&self, indices: &'a [u16]) -> &'a [u16] {
        &indices[self.first_vertex as usize..(self.first_vertex + self.num_vertices) as usize]
    }
}

/// Settings for constructing a [`ConvexHull3d`].
#[derive(Clone, Copy, Debug)]
pub struct ConvexHull3dSettings {
    /// The maximum angle (in radians) between the normals of two adjacent triangle faces
    /// for them to be considered coplanar and merged.
    ///
    /// **Default:** `1e-4`
    pub coplanarity_tolerance: f32,

    /// The maximum number of iterations for the Quickhull algorithm.
    ///
    /// This is a safeguard against infinite loops in degenerate cases, but should not be hit in practice.
    ///
    /// **Default:** `10_000`
    pub max_iterations: usize,
}

impl Default for ConvexHull3dSettings {
    fn default() -> Self {
        Self {
            coplanarity_tolerance: 1e-4,
            max_iterations: 10_000,
        }
    }
}

impl ConvexHull3d {
    /// Attempts to compute a [`ConvexHull3d`] for the given set of points.
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
        settings: ConvexHull3dSettings,
    ) -> Result<Self, ConvexHull3dError> {
        let mesh = ConvexTriangleMesh::try_from_points(points, Some(settings.max_iterations))?;
        Ok(Self::from_convex_mesh(
            &mesh,
            settings.coplanarity_tolerance,
        ))
    }

    /// Creates a [`ConvexHull3d`] from a [`ConvexTriangleMesh`] by merging
    /// coplanar triangular faces into polygonal faces.
    ///
    /// `coplanarity_tolerance` is the maximum angle (in radians) between the normals
    /// of two adjacent triangle faces for them to be considered coplanar and merged.
    ///
    /// # Example
    ///
    /// ```
    /// use glam::Vec3A;
    /// use quickhull::{ConvexHull3d, ConvexTriangleMesh};
    ///
    /// // A cube's convex hull has 12 triangles (2 per quad face).
    /// let points = vec![
    ///     Vec3A::new( 1.0,  1.0,  1.0),
    ///     Vec3A::new( 1.0,  1.0, -1.0),
    ///     Vec3A::new( 1.0, -1.0,  1.0),
    ///     Vec3A::new( 1.0, -1.0, -1.0),
    ///     Vec3A::new(-1.0,  1.0,  1.0),
    ///     Vec3A::new(-1.0,  1.0, -1.0),
    ///     Vec3A::new(-1.0, -1.0,  1.0),
    ///     Vec3A::new(-1.0, -1.0, -1.0),
    /// ];
    ///
    /// // Compute the convex hull with triangle faces.
    /// let mesh = ConvexTriangleMesh::try_from_points(&points, None).unwrap();
    /// assert_eq!(mesh.indices().len(), 12);
    ///
    /// // Construct a convex hull with polygonal faces from the triangle mesh,
    /// // merging triangles with normals within 1e-4 radians of each other.
    /// let hull = ConvexHull3d::from_convex_mesh(&mesh, 1e-4);
    ///
    /// // The cube's hull should have 8 vertices and 6 quad faces after merging.
    /// assert_eq!(hull.faces().len(), 6);
    /// ```
    pub fn from_convex_mesh(mesh: &ConvexTriangleMesh, coplanarity_tolerance: f32) -> Self {
        let tri_indices = mesh.indices();
        let points = mesh.vertices();

        if tri_indices.is_empty() {
            return ConvexHull3d::default();
        }

        let cos_tolerance = coplanarity_tolerance
            .clamp(0.0, core::f32::consts::PI)
            .cos();
        let num_tris = tri_indices.len();

        // Compute triangle normals.
        let tri_normals: Vec<Vec3A> = tri_indices
            .iter()
            .map(|tri| {
                // TODO: Should we use the shortest edge? Is the additional precision worth it?
                //       https://box2d.org/posts/2014/01/troublesome-triangle/
                let a = points[tri[0] as usize];
                let b = points[tri[1] as usize];
                let c = points[tri[2] as usize];
                (b - a).cross(c - a).normalize_or_zero()
            })
            .collect();

        struct InternalEdge {
            v0: u32,
            v1: u32,
            tri_a: u32,
            tri_b: u32,
        }

        // Compute edges from triangle adjacency.
        let mut candidates: Vec<InternalEdge> = Vec::new();
        let mut edge_map: HashMap<(u32, u32), u32> = HashMap::new();
        for (face_idx, tri) in tri_indices.iter().enumerate() {
            let face_idx = face_idx as u32;
            for ei in 0..3 {
                let v0 = tri[ei];
                let v1 = tri[(ei + 1) % 3];
                let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };
                if let Some(&other_face_idx) = edge_map.get(&key) {
                    candidates.push(InternalEdge {
                        v0: key.0,
                        v1: key.1,
                        tri_a: other_face_idx,
                        tri_b: face_idx,
                    });
                } else {
                    edge_map.insert(key, face_idx);
                }
            }
        }

        // Initialize each triangle as a separate face.
        let mut face_vertices: Vec<Vec<u32>> = tri_indices
            .iter()
            .map(|tri| vec![tri[0], tri[1], tri[2]])
            .collect();
        let mut face_normals: Vec<Vec3A> = tri_normals;
        let mut face_parent: Vec<u32> = (0..num_tris as u32).collect();

        // Reusable buffers.
        let mut seen = vec![false; points.len()];
        let mut merge_buffer = Vec::new();

        // Greedily merge coplanar faces.
        for edge in &candidates {
            let face_a = Self::find_root(&mut face_parent, edge.tri_a) as usize;
            let face_b = Self::find_root(&mut face_parent, edge.tri_b) as usize;

            if face_a == face_b {
                continue;
            }

            if face_normals[face_a].dot(face_normals[face_b]) < cos_tolerance {
                continue;
            }

            if Self::try_merge_faces(
                &face_vertices[face_a],
                &face_vertices[face_b],
                edge.v0,
                edge.v1,
                points,
                face_normals[face_a],
                &mut seen,
                &mut merge_buffer,
            ) {
                face_normals[face_a] = Self::polygon_normal(&merge_buffer, points);
                core::mem::swap(&mut face_vertices[face_a], &mut merge_buffer);
                face_vertices[face_b].clear();
                face_parent[face_b] = face_a as u32;
            }
        }

        // Collect alive faces into the hull representation.
        let mut indices: Vec<u16> = Vec::new();
        let mut faces: Vec<HullFace> = Vec::new();
        let mut planes: Vec<Plane3d> = Vec::new();

        for face_idx in 0..num_tris {
            if face_parent[face_idx] as usize != face_idx {
                continue;
            }
            let vertices = &face_vertices[face_idx];
            let first = indices.len() as u16;
            for &v in vertices {
                indices.push(v as u16);
            }
            let normal = face_normals[face_idx];
            planes.push(Plane3d::from_point_and_normal(
                points[vertices[0] as usize],
                normal,
            ));
            faces.push(HullFace {
                first_vertex: first,
                num_vertices: vertices.len() as u16,
            });
        }

        ConvexHull3d {
            vertices: points.to_vec(),
            indices,
            faces,
            planes,
        }
    }

    /// Returns the vertices of the convex hull.
    #[inline]
    pub fn vertices(&self) -> &[Vec3A] {
        &self.vertices
    }

    /// Returns the indices of the convex hull vertices that define the faces of the hull.
    ///
    /// Each face is defined by a contiguous sequence of vertex indices, with the starting index
    /// and number of vertices for each face specified in [`faces`](Self::faces).
    #[inline]
    pub fn indices(&self) -> &[u16] {
        &self.indices
    }

    /// Returns the faces of the convex hull, defined by the indices in `indices`.
    #[inline]
    pub fn faces(&self) -> &[HullFace] {
        &self.faces
    }

    /// Returns the bounding planes associated with each face of the convex hull.
    #[inline]
    pub fn planes(&self) -> &[Plane3d] {
        &self.planes
    }

    /// Returns the vertices, indices, faces, and planes of the convex hull.
    ///
    /// This consumes the hull and allows taking ownership of the underlying data without cloning.
    #[inline]
    pub fn into_parts(self) -> (Vec<Vec3A>, Vec<u16>, Vec<HullFace>, Vec<Plane3d>) {
        (self.vertices, self.indices, self.faces, self.planes)
    }

    /// Attempts to merge two convex face loops across a shared edge.
    ///
    /// On success, writes the merged vertex loop into `out` and returns `true`.
    /// On failure (concave or degenerate result), returns `false` without
    /// modifying `out`.
    ///
    /// The `seen` buffer must have a length greater than or equal to
    /// the number of points and be zeroed.
    #[allow(clippy::too_many_arguments)]
    fn try_merge_faces(
        a: &[u32],
        b: &[u32],
        edge_v0: u32,
        edge_v1: u32,
        points: &[Vec3A],
        normal: Vec3A,
        seen: &mut [bool],
        out: &mut Vec<u32>,
    ) -> bool {
        let a_len = a.len();
        let b_len = b.len();

        // Find the shared edge positions in both face loops.
        let (a_idx, b_idx) = match Self::find_shared_edge(a, b, edge_v0, edge_v1) {
            Some(v) => v,
            None => return false,
        };

        // O(1) convexity check: only the 2 seam vertices change angle.
        // - At v0: (a_prev, v0, v1) becomes (a_prev, v0, b_next)
        // - At v1: (v0, v1, a_next) becomes (b_prev, v1, a_next)
        let a_prev = a[(a_idx + a_len - 1) % a_len];
        let a_next = a[(a_idx + 2) % a_len];
        let b_prev = b[(b_idx + b_len - 1) % b_len];
        let b_next = b[(b_idx + 2) % b_len];

        let p_v0 = points[a[a_idx] as usize];
        let p_v1 = points[a[(a_idx + 1) % a_len] as usize];

        // Check that the merged angle at v0 and v1 is convex
        if (p_v0 - points[a_prev as usize])
            .cross(points[b_next as usize] - p_v0)
            .dot(normal)
            < 0.0
        {
            return false;
        }
        if (p_v1 - points[b_prev as usize])
            .cross(points[a_next as usize] - p_v1)
            .dot(normal)
            < 0.0
        {
            return false;
        }

        // Check for shared vertices beyond the shared edge.
        // This avoids creating duplicate vertices in the merged polygon.
        for k in 0..(a_len - 1) {
            seen[a[(a_idx + 2 + k) % a_len] as usize] = true;
        }
        let mut has_duplicates = false;
        for k in 0..(b_len - 1) {
            if seen[b[(b_idx + 2 + k) % b_len] as usize] {
                has_duplicates = true;
                break;
            }
        }
        for k in 0..(a_len - 1) {
            seen[a[(a_idx + 2 + k) % a_len] as usize] = false;
        }
        if has_duplicates {
            return false;
        }

        // Build the merged polygon into the output buffer.
        out.clear();
        out.reserve(a_len + b_len - 2);
        for k in 0..(a_len - 1) {
            out.push(a[(a_idx + 2 + k) % a_len]);
        }
        for k in 0..(b_len - 1) {
            out.push(b[(b_idx + 2 + k) % b_len]);
        }

        true
    }

    /// Finds the positions of a shared edge in two face loops.
    ///
    /// Returns `(a_idx, b_idx)` such that `a[a_idx] -> a[(a_idx+1) % a_len]` and
    /// `b[b_idx] -> b[(b_idx+1) % b_len]` are the same edge in opposite directions.
    fn find_shared_edge(
        a: &[u32],
        b: &[u32],
        edge_v0: u32,
        edge_v1: u32,
    ) -> Option<(usize, usize)> {
        for i in 0..a.len() {
            let a0 = a[i];
            let a1 = a[(i + 1) % a.len()];

            if (a0 == edge_v0 && a1 == edge_v1) || (a0 == edge_v1 && a1 == edge_v0) {
                // Find the matching reverse edge in B.
                for j in 0..b.len() {
                    if b[j] == a1 && b[(j + 1) % b.len()] == a0 {
                        return Some((i, j));
                    }
                }
            }
        }

        None
    }

    /// Finds the root of the union-find structure for face merging.
    fn find_root(parent: &mut [u32], mut i: u32) -> u32 {
        while parent[i as usize] != i {
            parent[i as usize] = parent[parent[i as usize] as usize];
            i = parent[i as usize];
        }
        i
    }

    /// Computes the outward normal of a polygon using [Newell's method].
    ///
    /// [Newell's method]: https://wikis.khronos.org/opengl/Calculating_a_Surface_Normal
    fn polygon_normal(indices: &[u32], vertices: &[Vec3A]) -> Vec3A {
        let n = indices.len();
        if n < 3 {
            return Vec3A::ZERO;
        }

        let v0 = vertices[indices[0] as usize];
        let mut sum = Vec3A::ZERO;

        for i in 1..(n - 1) {
            let edge1 = vertices[indices[i] as usize] - v0;
            let edge2 = vertices[indices[i + 1] as usize] - v0;
            sum += edge1.cross(edge2);
        }

        sum.normalize_or_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ConvexTriangleMesh;

    fn cube_points() -> Vec<Vec3A> {
        vec![
            Vec3A::new(1.0, 1.0, 1.0),
            Vec3A::new(1.0, 1.0, -1.0),
            Vec3A::new(1.0, -1.0, 1.0),
            Vec3A::new(1.0, -1.0, -1.0),
            Vec3A::new(-1.0, 1.0, 1.0),
            Vec3A::new(-1.0, 1.0, -1.0),
            Vec3A::new(-1.0, -1.0, 1.0),
            Vec3A::new(-1.0, -1.0, -1.0),
        ]
    }

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
            let angle_x = 2.0 * core::f32::consts::PI * (step_x as f32 / divisions as f32);
            let p = rot_x(unit_y, angle_x);
            for step_z in 0..divisions {
                let angle_z = 2.0 * core::f32::consts::PI * (step_z as f32 / divisions as f32);
                let p = rot_z(p, angle_z);
                points.push(p);
            }
        }
        points
    }

    #[test]
    fn cube_merge_faces_test() {
        let mesh = ConvexTriangleMesh::try_from_points(&cube_points(), None).unwrap();
        assert_eq!(mesh.indices().len(), 12);

        let hull = ConvexHull3d::from_convex_mesh(&mesh, 1e-3);
        assert_eq!(hull.faces().len(), 6);
        for face in hull.faces() {
            assert_eq!(face.num_vertices(), 4, "each cube face should be a quad");
        }
    }

    #[test]
    fn cube_planes_test() {
        let hull =
            ConvexHull3d::try_from_points(&cube_points(), ConvexHull3dSettings::default()).unwrap();

        // Each plane should have unit normal.
        for plane in hull.planes() {
            let normal_len = plane.normal().length();
            assert!(
                (normal_len - 1.0).abs() < 1e-5,
                "plane normal should be unit length, got {normal_len}"
            );
        }

        // Each vertex of the face should lie on its plane.
        for (face, plane) in hull.faces().iter().zip(hull.planes().iter()) {
            let indices = face.vertex_indices(hull.indices());
            for &idx in indices {
                let point = hull.vertices()[idx as usize];
                let dist = plane.signed_distance_to_point(point).abs();
                assert!(
                    dist < 1e-4,
                    "vertex should lie on face plane, distance = {dist}"
                );
            }
        }
    }

    #[test]
    fn octahedron_merge_faces_test() {
        let points = vec![
            Vec3A::new(1.0, 0.0, 0.0),
            Vec3A::new(0.0, 1.0, 0.0),
            Vec3A::new(0.0, 0.0, 1.0),
            Vec3A::new(-1.0, 0.0, 0.0),
            Vec3A::new(0.0, -1.0, 0.0),
            Vec3A::new(0.0, 0.0, -1.0),
        ];
        let hull = ConvexHull3d::try_from_points(&points, ConvexHull3dSettings::default()).unwrap();
        assert_eq!(hull.faces().len(), 8);
        for face in hull.faces() {
            assert_eq!(
                face.num_vertices(),
                3,
                "each octahedron face should remain a triangle"
            );
        }
    }

    #[test]
    fn sphere_merge_faces_test() {
        let points = sphere_points(10);
        let mesh = ConvexTriangleMesh::try_from_points(&points, None).unwrap();
        let tri_count = mesh.indices().len();
        let hull = ConvexHull3d::from_convex_mesh(&mesh, 1e-3);

        // Merging can only reduce (or maintain) the face count.
        assert!(hull.faces().len() <= tri_count);

        for face in hull.faces() {
            assert!(face.num_vertices() >= 3);
        }
    }

    #[test]
    fn empty_hull_test() {
        let points: Vec<Vec3A> = Vec::new();
        let hull = ConvexHull3d::try_from_points(&points, ConvexHull3dSettings::default()).unwrap();
        assert!(hull.vertices().is_empty());
        assert!(hull.indices().is_empty());
        assert!(hull.faces().is_empty());
        assert!(hull.planes().is_empty());
    }

    // With aggressive merging, a group can surround another group, creating
    // a multi-loop boundary. The algorithm should detect this and fall back
    // to individual triangles rather than producing a disjoint polygon face.
    #[test]
    fn aggressive_merge_no_disjoint_faces() {
        let points = sphere_points(10);
        let mesh = ConvexTriangleMesh::try_from_points(&points, None).unwrap();

        // Test a range of aggressive tolerances.
        for tolerance in [0.5, 1.0, 1.5, core::f32::consts::PI] {
            let hull = ConvexHull3d::from_convex_mesh(&mesh, tolerance);

            for (face_idx, face) in hull.faces().iter().enumerate() {
                let vertex_indices = face.vertex_indices(hull.indices());

                // Every consecutive pair of vertices in the face should be distinct.
                for i in 0..vertex_indices.len() {
                    let v0 = vertex_indices[i];
                    let v1 = vertex_indices[(i + 1) % vertex_indices.len()];
                    assert_ne!(
                        v0, v1,
                        "face {face_idx} has duplicate consecutive vertex {v0}"
                    );
                }

                // No vertex should appear more than once.
                let mut sorted = vertex_indices.to_vec();
                sorted.sort();
                sorted.dedup();
                assert_eq!(
                    sorted.len(),
                    vertex_indices.len(),
                    "face {face_idx} has repeated vertices (tolerance={tolerance})"
                );
            }
        }
    }

    // All merged polygon faces must be convex. With a non-zero tolerance,
    // drift could merge triangles that curve around the hull surface, but the algorithm
    // should detect this and fall back to individual triangles.
    #[test]
    fn merged_faces_are_convex() {
        let points = sphere_points(10);
        let mesh = ConvexTriangleMesh::try_from_points(&points, None).unwrap();

        for tolerance in [1e-4, 0.1, 0.5, 1.0, 1.5, core::f32::consts::PI] {
            let hull = ConvexHull3d::from_convex_mesh(&mesh, tolerance);

            for (face_idx, (face, plane)) in hull.faces().iter().zip(hull.planes()).enumerate() {
                let vertex_indices = face.vertex_indices(hull.indices());
                let n = vertex_indices.len();
                let face_normal = plane.normal();

                for i in 0..n {
                    let a = hull.vertices()[vertex_indices[i] as usize];
                    let b = hull.vertices()[vertex_indices[(i + 1) % n] as usize];
                    let c = hull.vertices()[vertex_indices[(i + 2) % n] as usize];
                    let cross = (b - a).cross(c - b);
                    assert!(
                        cross.dot(face_normal) >= -1e-6,
                        "face {face_idx} is concave at vertex {} (tolerance={tolerance})",
                        (i + 1) % n,
                    );
                }
            }
        }
    }

    // Face count must be monotonically non-increasing as tolerance grows.
    #[test]
    fn face_count_monotonically_decreasing() {
        let points = sphere_points(10);
        let mesh = ConvexTriangleMesh::try_from_points(&points, None).unwrap();

        let tolerances = [0.0, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1];
        let mut prev_count = usize::MAX;

        for tolerance in tolerances {
            let hull = ConvexHull3d::from_convex_mesh(&mesh, tolerance);
            let count = hull.faces().len();
            assert!(
                count <= prev_count,
                "face count increased from {prev_count} to {count} at tolerance={tolerance}"
            );
            prev_count = count;
        }
    }
}
