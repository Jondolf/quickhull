use glam::Vec3A;
use hashbrown::HashMap;

use super::plane::Plane3d;
use crate::{ConvexHull3dError, ConvexTriangleMesh};

/// A 3D [convex hull] representing the smallest convex set containing
/// all input points in a given point set.
///
/// This can be thought of as a shrink wrapping of a 3D object.
///
/// Unlike [`ConvexTriangleMesh`], which always has triangular faces,
/// [`ConvexHull3d`] supports polygonal faces formed by merging coplanar triangles.
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
    /// `angle_tolerance` is the maximum angle (in radians) between the normals
    /// of two adjacent triangle faces for them to be considered coplanar and merged.
    ///
    /// A typical value is `1e-4` for near-exact coplanarity, or `1e-2` for more
    /// aggressive merging.
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
    /// let mesh = ConvexTriangleMesh::try_from_points(&points, None).unwrap();
    /// assert_eq!(mesh.indices().len(), 12);
    ///
    /// // Merge coplanar faces into polygons.
    /// let hull = ConvexHull3d::from_convex_mesh(&mesh, 1e-4);
    /// assert_eq!(hull.faces().len(), 6); // 6 quad faces
    /// ```
    pub fn from_convex_mesh(mesh: &ConvexTriangleMesh, angle_tolerance: f32) -> Self {
        let tri_indices = mesh.indices();
        let points = mesh.points();

        if tri_indices.is_empty() {
            return ConvexHull3d::default();
        }

        // Build triangle adjacency from indices.
        // `tri_adjacency[i][e]` is the index of the neighboring face
        // across edge `e` (0, 1, 2) of face `i`,
        let tri_adjacency = Self::build_triangle_adjacency(tri_indices);

        let cos_tolerance = angle_tolerance.cos();
        let num_tris = tri_indices.len();

        // Compute triangle normals.
        let normals: Vec<Vec3A> = tri_indices
            .iter()
            .map(|tri| {
                let a = points[tri[0] as usize];
                let b = points[tri[1] as usize];
                let c = points[tri[2] as usize];

                // The most accurate normal is calculated by using the two shortest edges
                // from their shared vertex, which avoids precision loss for thin triangles.
                // https://box2d.org/posts/2014/01/troublesome-triangle/
                let ab = b - a;
                let ac = c - a;
                let bc = c - b;
                let ab_sq = ab.length_squared();
                let ac_sq = ac.length_squared();
                let bc_sq = bc.length_squared();

                // Find the longest edge and use the opposite vertex.
                // The cross product order preserves CCW winding.
                let normal = if bc_sq >= ab_sq && bc_sq >= ac_sq {
                    // BC is longest -> use vertex A
                    ab.cross(ac)
                } else if ac_sq >= ab_sq {
                    // AC is longest -> use vertex B
                    bc.cross(-ab)
                } else {
                    // AB is longest -> use vertex C
                    (-ac).cross(-bc)
                };
                normal.normalize_or_zero()
            })
            .collect();

        // Flood-fill coplanar groups using triangle adjacency.
        // Each "group" is a set of coplanar triangles that will be merged into a single face of the hull.
        // `tri_to_group[i]` is the group ID of triangle `i`, or `u32::MAX` if it hasn't been assigned to a group yet.
        let mut tri_to_group = vec![u32::MAX; num_tris];
        let mut num_groups: u32 = 0;
        let mut stack: Vec<usize> = Vec::new();

        // Track one boundary edge per group for boundary loop extraction.
        // `group_starts[g]` is `(tri_index, edge_index)` of a boundary edge for group `g`,
        // or `(usize::MAX, 0)` for degenerate groups.
        let mut group_starts: Vec<(usize, usize)> = Vec::new();

        for start in 0..num_tris {
            if tri_to_group[start] != u32::MAX {
                continue;
            }

            let group_id = num_groups;
            num_groups += 1;
            tri_to_group[start] = group_id;
            stack.push(start);

            let mut found_boundary = false;

            while let Some(tri) = stack.pop() {
                for (edge, &neighbor_tri) in tri_adjacency[tri].iter().enumerate() {
                    let neighbor = neighbor_tri as usize;

                    if tri_to_group[neighbor] == u32::MAX {
                        if normals[tri].dot(normals[neighbor]) >= cos_tolerance {
                            tri_to_group[neighbor] = group_id;
                            stack.push(neighbor);
                        } else if !found_boundary {
                            group_starts.push((tri, edge));
                            found_boundary = true;
                        }
                    } else if tri_to_group[neighbor] != group_id && !found_boundary {
                        group_starts.push((tri, edge));
                        found_boundary = true;
                    }
                }
            }

            if !found_boundary {
                // All faces are in one group (degenerate hull with no boundary edges).
                group_starts.push((usize::MAX, 0));
            }
        }

        // Build the merged hull.
        let mut indices: Vec<u16> = Vec::new();
        let mut faces: Vec<HullFace> = Vec::with_capacity(num_groups as usize);
        let mut planes: Vec<Plane3d> = Vec::with_capacity(num_groups as usize);

        for group in 0..num_groups {
            let (start_tri, start_edge) = group_starts[group as usize];

            if start_tri == usize::MAX {
                // Degenerate: the entire mesh is one group with no boundary edges.
                // Emit each triangle as an individual face.
                for (i, &face_group) in tri_to_group.iter().enumerate() {
                    if face_group == group {
                        let first = indices.len() as u16;
                        let tri = &tri_indices[i];

                        for &vi in tri {
                            indices.push(vi as u16);
                        }

                        let normal = normals[i];
                        planes.push(Plane3d::from_point_and_normal(
                            points[tri[0] as usize],
                            normal,
                        ));
                        faces.push(HullFace {
                            first_vertex: first,
                            num_vertices: 3,
                        });
                    }
                }
                continue;
            }

            let first = indices.len() as u16;

            // Extract boundary loop.
            let loop_vertex_indices = Self::extract_boundary_loop(
                tri_indices,
                &tri_adjacency,
                &tri_to_group,
                group,
                start_tri,
                start_edge,
            );

            for &i in loop_vertex_indices.iter() {
                indices.push(i as u16);
            }

            // Compute plane from the group's normal and a point on the face.
            let normal = normals[start_tri];
            let point_on_face = points[loop_vertex_indices[0] as usize];
            planes.push(Plane3d::from_point_and_normal(point_on_face, normal));

            faces.push(HullFace {
                first_vertex: first,
                num_vertices: loop_vertex_indices.len() as u16,
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

    /// Builds triangle adjacency from indices using an edge map.
    ///
    /// `adjacency[i][e]` is the index of the neighboring face across edge `e` of face `i`,
    /// where edge `e` goes from `indices[i][e]` to `indices[i][(e + 1) % 3]`.
    fn build_triangle_adjacency(indices: &[[u32; 3]]) -> Vec<[u32; 3]> {
        let mut adjacency = vec![[u32::MAX; 3]; indices.len()];
        let mut edge_map: HashMap<(u32, u32), (u32, u32)> = HashMap::new();

        for (face_i, tri) in indices.iter().enumerate() {
            for edge_i in 0..3u32 {
                let v0 = tri[edge_i as usize];
                let v1 = tri[((edge_i + 1) % 3) as usize];
                let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };

                if let Some(&(other_fi, other_ei)) = edge_map.get(&key) {
                    adjacency[face_i][edge_i as usize] = other_fi;
                    adjacency[other_fi as usize][other_ei as usize] = face_i as u32;
                } else {
                    edge_map.insert(key, (face_i as u32, edge_i));
                }
            }
        }

        adjacency
    }

    /// Extracts the ordered boundary loop of a face group by walking around its edges.
    ///
    /// Uses the triangle adjacency to traverse the triangle fan at each vertex
    /// without any additional allocations.
    ///
    /// The time complexity is O(V) where V is the number of vertices in the face,
    /// since each edge is visited at most twice.
    fn extract_boundary_loop(
        indices: &[[u32; 3]],
        adjacency: &[[u32; 3]],
        tri_to_group: &[u32],
        group_id: u32,
        start_tri: usize,
        start_edge: usize,
    ) -> Vec<u32> {
        let mut result = Vec::new();
        let mut tri = start_tri;
        let mut edge = start_edge;

        loop {
            // The current boundary edge goes from indices[tri][edge]
            // to indices[tri][(edge + 1) % 3].
            result.push(indices[tri][edge]);

            // Walk around vertex V = indices[tri][(edge + 1) % 3] through
            // the triangle fan to find the next boundary edge.
            let v = indices[tri][(edge + 1) % 3];
            let mut current_tri = tri;
            let mut current_edge = (edge + 1) % 3;

            loop {
                let neighbor = adjacency[current_tri][current_edge] as usize;
                if tri_to_group[neighbor] != group_id {
                    // This edge crosses a group boundary; it's the next boundary edge.
                    tri = current_tri;
                    edge = current_edge;
                    break;
                }

                // Cross into the neighbor triangle and find the shared edge.
                let w = indices[current_tri][(current_edge + 1) % 3];
                let ntri = indices[neighbor];
                let nedge = if ntri[0] == w && ntri[1] == v {
                    0
                } else if ntri[1] == w && ntri[2] == v {
                    1
                } else {
                    debug_assert!(ntri[2] == w && ntri[0] == v, "broken mesh adjacency");
                    2
                };

                // Continue walking around vertex V from the next edge.
                current_tri = neighbor;
                current_edge = (nedge + 1) % 3;
            }

            if tri == start_tri && edge == start_edge {
                break;
            }
        }

        result
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
    fn cube_merge_faces_test() {
        let mesh = ConvexTriangleMesh::try_from_points(&cube_points(), None).unwrap();
        assert_eq!(mesh.indices().len(), 12);

        let hull = ConvexHull3d::from_convex_mesh(&mesh, 1e-4);
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
        let hull = ConvexHull3d::from_convex_mesh(&mesh, 1e-4);

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
}
