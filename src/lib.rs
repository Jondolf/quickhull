//! # Quickhull
//!
//! A Rust-implementation of the Quickhull algorithm for computing convex hulls for point sets.
//!
//! This is a simplified and cleaned up version of [chull](https://github.com/u65xhd/chull),
//! focusing on making the algorithm robust and efficient for the 2D and 3D cases.
//!
//! ## References
//!
//! - C. Bradford Barber et al. 1996. [The Quickhull Algorithm for Convex Hulls](https://www.cise.ufl.edu/~ungor/courses/fall06/papers/QuickHull.pdf) (the original paper)
//! - Dirk Gregorius. GDC 2014. [Physics for Game Programmers: Implementing Quickhull](https://archive.org/details/GDC2014Gregorius)

#![warn(missing_docs)]

use glam::{DMat4, DVec3, DVec4};

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::error::Error;
use std::fmt;

/// A polygonal face belonging to a [`ConvexHull`].
#[derive(Debug, Clone)]
struct Face {
    /// The indices of the face's points.
    indices: Vec<usize>,
    /// The indices of points in front of the face plane, or the points that can "see" the face,
    /// and the distance to each of those points along the normal.
    outside_points: Vec<(usize, f64)>,
    /// The indices of neighboring faces.
    neighbor_faces: Vec<usize>,
    /// The normal of the face.
    normal: DVec3,
    /// How far away from the origin this face is along its normal.
    distance_from_origin: f64,
}

impl Face {
    fn new(points: &[DVec3], indices: [usize; 3]) -> Self {
        let points_of_face = indices.map(|i| points[i]);
        let normal = face_normal(points_of_face);
        let origin = normal.dot(points_of_face[0]);

        Self {
            indices: indices.to_vec(),
            outside_points: Vec::new(),
            neighbor_faces: Vec::new(),
            normal,
            distance_from_origin: origin,
        }
    }
}

/// The type of error returned during [`ConvexHull`] construction.
#[derive(Debug, Clone, PartialEq)]
pub enum ErrorKind {
    /// The given point set is empty, so no convex hull could be computed.
    Empty,
    /// The convex hull generation algorithm encountered degeneracies.
    Degenerated,
    /// The given point set cannot produce a valid convex hull.
    DegenerateInput(DegenerateInput),
    /// A round-off error.
    RoundOffError(String),
}

/// The type of degeneracy for when attempting to compute a convex hull for a point set.
#[derive(Debug, Clone, PartialEq)]
pub enum DegenerateInput {
    /// The input points are approximately equal.
    Coincident,
    /// The input points are approximately on the same line.
    Collinear,
    /// The input points are approximately on the same plane.
    Coplanar,
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self {
            ErrorKind::Empty => write!(f, "empty"),
            ErrorKind::Degenerated => write!(f, "degenerated"),
            ErrorKind::DegenerateInput(kind) => write!(f, "degenerate input: {:?}", kind),
            ErrorKind::RoundOffError(msg) => {
                write!(f, "erroneous results by roundoff error: {}", msg)
            }
        }
    }
}
impl Error for ErrorKind {}

/// A 3D convex hull representing the smallest convex set containing
/// all input points in a given point set.
///
/// This can be thought of as a shrink wrapping of a 3D object.
#[derive(Clone, Debug)]
pub struct ConvexHull {
    /// The points of the convex hull.
    pub points: Vec<DVec3>,
    /// The faces of the convex hull.
    faces: BTreeMap<usize, Face>,
}

impl ConvexHull {
    /// Attempts to compute a [`ConvexHull`] for the given set of points.
    pub fn try_new(
        points: &[DVec3],
        threshold: impl Into<f64>,
        max_iter: Option<usize>,
    ) -> Result<Self, ErrorKind> {
        let threshold = threshold.into();
        let num_points = points.len();

        if num_points == 0 {
            return Err(ErrorKind::Empty);
        }

        if num_points <= 3 || is_degenerate(points, threshold) {
            return Err(ErrorKind::Degenerated);
        }

        // Create the initial simplex, a tetrahedron in 3D.
        let mut c_hull = Self::create_simplex(points, threshold)?;

        // Run the main quick hull algorithm.
        c_hull.update(threshold, max_iter)?;

        // Shrink the hull, removing unused points.
        c_hull.remove_unused_points();

        if c_hull.points.len() <= 3 {
            return Err(ErrorKind::Degenerated);
        }
        Ok(c_hull)
    }

    fn create_simplex(points: &[DVec3], threshold: f64) -> Result<Self, ErrorKind> {
        let indices_set = Self::select_vertices_for_simplex(points)?;
        let mut faces = BTreeMap::new();

        let mut face_add_count = 0;

        #[allow(clippy::explicit_counter_loop)]
        for i_face in 0..4 {
            let mut face_indices = Vec::new();
            // create face
            for (j, set) in indices_set.iter().enumerate().take(4) {
                if j != i_face {
                    face_indices.push(*set);
                }
            }
            let mut face = Face::new(points, face_indices.try_into().unwrap());

            // Check the order of the face's vertices.
            let rem_point = indices_set[i_face];
            let pos = position_from_face(points, &face, rem_point);
            if pos > threshold {
                face.indices.swap(0, 1);
                face.normal = -face.normal;
                face.distance_from_origin = -face.distance_from_origin;
            }
            if face.indices.len() != 3 {
                return Err(ErrorKind::RoundOffError(
                    "number of face's vertices should be 3".to_string(),
                ));
            }
            faces.insert(face_add_count, face);
            face_add_count += 1;
        }

        // Link neighbors.
        let simplex_face_key: Vec<_> = faces.keys().copied().collect();
        for (key, face) in &mut faces.iter_mut() {
            face.outside_points
                .sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
            for neighbors_key in simplex_face_key
                .iter()
                .filter(|neighbor_key| *neighbor_key != key)
            {
                face.neighbor_faces.push(*neighbors_key);
            }
        }
        let simplex = Self {
            points: points.to_vec(),
            faces,
        };

        Ok(simplex)
    }

    fn compute_min_and_max(points: &[DVec3]) -> ((DVec3, [usize; 3]), (DVec3, [usize; 3])) {
        let mut min = points[0];
        let mut max = points[0];

        let mut min_vertices = [0; 3];
        let mut max_vertices = [0; 3];

        for (i, vtx) in points.iter().enumerate().skip(1) {
            if vtx.x < min.x {
                min.x = vtx.x;
                min_vertices[0] = i;
            } else if vtx.x > max.x {
                max.x = vtx.x;
                max_vertices[0] = i;
            }

            if vtx.y < min.y {
                min.y = vtx.y;
                min_vertices[1] = i;
            } else if vtx.y > max.y {
                max.y = vtx.y;
                max_vertices[1] = i;
            }

            if vtx.z < min.z {
                min.z = vtx.z;
                min_vertices[2] = i;
            } else if vtx.z > max.z {
                max.z = vtx.z;
                max_vertices[2] = i;
            }
        }

        ((min, min_vertices), (max, max_vertices))
    }

    fn tolerance(min: DVec3, max: DVec3) -> f64 {
        3.0 * f64::EPSILON * max.abs().max(min.abs()).element_sum()
    }

    fn update(&mut self, threshold: f64, max_iter: Option<usize>) -> Result<(), ErrorKind> {
        let mut face_add_count = *self.faces.iter().last().map(|(k, _v)| k).unwrap() + 1;
        let mut num_iter = 0;
        let mut assigned_point_indices: HashSet<usize> = HashSet::new();

        for face in self.faces.values() {
            for index in &face.indices {
                assigned_point_indices.insert(*index);
            }
        }

        // Initialize the outside points.
        for (_key, face) in &mut self.faces.iter_mut() {
            for (i, _point) in self.points.iter().enumerate() {
                if assigned_point_indices.contains(&i) {
                    continue;
                }
                let pos = position_from_face(&self.points, face, i);
                if pos > threshold {
                    face.outside_points.push((i, pos));
                }
            }
        }

        let (max_iter, truncate) = if let Some(iter) = max_iter {
            (iter, true)
        } else {
            (0, false)
        };

        // The main algorithm of quick hull.
        while let Some((key, face)) = self
            .faces
            .iter()
            .find(|(_, face)| !face.outside_points.is_empty())
            .map(|(a, b)| (*a, b))
        {
            if truncate && num_iter >= max_iter {
                break;
            }

            num_iter += 1;

            // Select the furthest point.
            let (furthest_point_index, _) = *face.outside_points.last().unwrap();

            // Initialize the visible set.
            let visible_set = initialize_visible_set(
                &self.points,
                furthest_point_index,
                &self.faces,
                key,
                face,
                threshold,
            );

            // Get the horizon.
            let horizon = get_horizon(&visible_set, &self.faces)?;

            // Create new face.
            let mut new_keys = Vec::new();
            for (ridge, unvisible) in horizon {
                let mut new_face = vec![furthest_point_index];

                assigned_point_indices.insert(furthest_point_index);

                for point in ridge {
                    new_face.push(point);
                    assigned_point_indices.insert(point);
                }

                if new_face.len() != 3 {
                    return Err(ErrorKind::RoundOffError(
                        "number of new face's vertices should be 3".to_string(),
                    ));
                }

                let mut new_face = Face::new(&self.points, new_face.try_into().unwrap());
                new_face.neighbor_faces.push(unvisible);

                let new_key = face_add_count;
                face_add_count += 1;

                self.faces.insert(new_key, new_face);
                let unvisible_faset = self.faces.get_mut(&unvisible).unwrap();
                unvisible_faset.neighbor_faces.push(new_key);
                new_keys.push(new_key);
            }

            if new_keys.len() < 3 {
                return Err(ErrorKind::RoundOffError(
                    "number of new faces should be grater than 3".to_string(),
                ));
            }

            // Link the face to its neighbor.
            for (i, key_a) in new_keys.iter().enumerate() {
                let points_of_new_face_a: HashSet<_> = self
                    .faces
                    .get(key_a)
                    .unwrap()
                    .indices
                    .iter()
                    .copied()
                    .collect();

                for key_b in new_keys.iter().skip(i + 1) {
                    let points_of_new_face_b: HashSet<_> = self
                        .faces
                        .get(key_b)
                        .unwrap()
                        .indices
                        .iter()
                        .copied()
                        .collect();

                    let num_intersection_points = points_of_new_face_a
                        .intersection(&points_of_new_face_b)
                        .collect::<Vec<_>>()
                        .len();

                    if num_intersection_points == 2 {
                        let face_a = self.faces.get_mut(key_a).unwrap();
                        face_a.neighbor_faces.push(*key_b);

                        let face_b = self.faces.get_mut(key_b).unwrap();
                        face_b.neighbor_faces.push(*key_a);
                    }
                }

                let face_a = self.faces.get(key_a).unwrap();
                if face_a.neighbor_faces.len() != 3 {
                    return Err(ErrorKind::RoundOffError(
                        "number of neighbors should be 3".to_string(),
                    ));
                }
            }

            // Check the order of the new face's vertices.
            for new_key in &new_keys {
                let new_face = self.faces.get(new_key).unwrap();
                let mut degenerate = true;

                for assigned_point_index in &assigned_point_indices {
                    let position =
                        position_from_face(&self.points, new_face, *assigned_point_index);

                    if position <= threshold && position >= -threshold {
                        continue;
                    } else if position > 0.0 {
                        let new_face = self.faces.get_mut(new_key).unwrap();
                        new_face.indices.swap(0, 1);
                        new_face.normal = -new_face.normal;
                        new_face.distance_from_origin = -new_face.distance_from_origin;
                        degenerate = false;
                        break;
                    }

                    degenerate = false;
                    break;
                }

                if degenerate {
                    return Err(ErrorKind::Degenerated);
                }
            }

            // Set outside points for each new face.
            let mut visible_faces = Vec::new();
            for visible in &visible_set {
                visible_faces.push(self.faces.get(visible).unwrap().clone());
            }

            for new_key in &new_keys {
                let new_face = self.faces.get_mut(new_key).unwrap();
                let mut checked_point_set = HashSet::new();

                for visible_face in &visible_faces {
                    for (outside_point_index, _) in visible_face.outside_points.iter() {
                        if assigned_point_indices.contains(outside_point_index)
                            || checked_point_set.contains(outside_point_index)
                        {
                            continue;
                        }

                        checked_point_set.insert(outside_point_index);

                        let pos = position_from_face(&self.points, new_face, *outside_point_index);
                        if pos > threshold {
                            new_face.outside_points.push((*outside_point_index, pos));
                        }
                    }
                }

                new_face
                    .outside_points
                    .sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
            }

            // Delete the visible faces.
            for visible in visible_set {
                let visible_face = self.faces.get(&visible).unwrap().clone();
                for neighbor_key in visible_face.neighbor_faces {
                    let neighbor = self.faces.get_mut(&neighbor_key).unwrap();
                    let index = neighbor
                        .neighbor_faces
                        .iter()
                        .enumerate()
                        .find(|(_, k)| **k == visible)
                        .map(|(i, _)| i)
                        .unwrap();
                    neighbor.neighbor_faces.swap_remove(index);
                }
                self.faces.remove(&visible);
            }
        }

        if !self.is_convex(threshold) {
            return Err(ErrorKind::RoundOffError("concave".to_string()));
        }

        Ok(())
    }

    /// Adds the given points to the point set, attempting to update the convex hull.
    pub fn add_points(
        &mut self,
        points: &[DVec3],
        threshold: impl Into<f64>,
    ) -> Result<(), ErrorKind> {
        let threshold = threshold.into();

        self.points.append(&mut points.to_vec());
        self.update(threshold, None)?;
        self.remove_unused_points();

        if self.points.len() <= 3 {
            return Err(ErrorKind::Degenerated);
        }

        Ok(())
    }

    /// Returns the vertices and indices of the convex hull.
    pub fn vertices_indices(&self) -> (Vec<DVec3>, Vec<usize>) {
        let mut indices = Vec::new();
        for face in self.faces.values() {
            for i in &face.indices {
                indices.push(*i);
            }
        }
        (self.points.to_vec(), indices)
    }

    pub(crate) fn remove_unused_points(&mut self) {
        let mut indices_list = BTreeSet::new();

        for face in self.faces.values() {
            for i in &face.indices {
                indices_list.insert(*i);
            }
        }

        let indices_list: BTreeMap<usize, usize> = indices_list
            .into_iter()
            .enumerate()
            .map(|(i, index)| (index, i))
            .collect();

        for face in self.faces.values_mut() {
            let mut new_face_indices = Vec::with_capacity(face.indices.len());
            for i in &face.indices {
                new_face_indices.push(*indices_list.get(i).unwrap());
            }
            std::mem::swap(&mut face.indices, &mut new_face_indices);
        }

        let mut vertices = Vec::new();

        for (index, _i) in indices_list.iter() {
            vertices.push(self.points[*index]);
        }

        self.points = vertices;
    }

    /// Computes the volume of the convex hull.
    pub fn volume(&self) -> f64 {
        let (hull_vertices, hull_indices) = self.vertices_indices();
        let reference_point = hull_vertices[hull_indices[0]].extend(1.0);
        let mut volume = 0.0;
        for i in (3..hull_indices.len()).step_by(3) {
            let mut mat = DMat4::ZERO;
            for j in 0..3 {
                let row = hull_vertices[hull_indices[i + j]].extend(1.0);
                *mat.col_mut(j) = row;
            }
            *mat.col_mut(3) = reference_point;
            volume += mat.determinant();
        }
        let factorial = {
            let mut result = 1.0;
            let mut m = 1.0 + 1.0;
            let mut n = 3;
            while n > 1 {
                result *= m;
                n -= 1;
                m += 1.0;
            }
            result
        };
        volume / factorial
    }

    /// Checks if the convex hull is convex with the given threshold.
    fn is_convex(&self, threshold: f64) -> bool {
        for face in self.faces.values() {
            if position_from_face(&self.points, face, 0) > threshold {
                return false;
            }
        }
        true
    }

    /// Computes the point on the convex hull that is furthest in the given direction.
    pub fn support_point(&self, direction: DVec3) -> DVec3 {
        let mut max = self.points[0].dot(direction);
        let mut index = 0;

        for (i, point) in self.points.iter().enumerate().skip(1) {
            let dot_product = point.dot(direction);
            if dot_product > max {
                max = dot_product;
                index = i;
            }
        }

        self.points[index]
    }

    fn select_vertices_for_simplex(points: &[DVec3]) -> Result<[usize; 4], ErrorKind> {
        // Find the minimum and maximum points of the bounding box and the indices
        // of the minimum and maximum points for each individual axis.
        let ((min_point, per_axis_min_indices), (max_point, per_axis_max_indices)) =
            Self::compute_min_and_max(points);
        let mut indices = [0; 4];

        let threshold = Self::tolerance(min_point, max_point);

        let mut max_difference = 0.0;
        let mut max_index = 0;

        for i in 0..3 {
            let difference =
                points[per_axis_max_indices[i]][i] - points[per_axis_min_indices[i]][i];
            if difference > max_difference {
                max_difference = difference;
                max_index = i;
            }
        }

        if max_difference <= threshold {
            return Err(ErrorKind::DegenerateInput(DegenerateInput::Coincident));
        }

        // The first two vertices are the ones farthest apart in one dimension.
        indices[0] = per_axis_max_indices[max_index];
        indices[1] = per_axis_min_indices[max_index];

        // The third vertex should be the one farthest from the line segment
        // between the first two vertices.
        let unit_01 = (points[indices[1]] - points[indices[0]]).normalize();
        let mut normal = DVec3::ZERO;

        let mut max_squared_distance = 0.0;

        for i in 0..points.len() {
            let diff = points[i] - points[indices[0]];
            let cross = unit_01.cross(diff);
            let distance_squared = cross.length_squared();

            if distance_squared > max_squared_distance
                && points[i] != points[indices[0]]
                && points[i] != points[indices[1]]
            {
                max_squared_distance = distance_squared;
                indices[2] = i;
                normal = cross;
            }
        }

        if max_squared_distance <= (100.0 * threshold).powi(2) {
            return Err(ErrorKind::DegenerateInput(DegenerateInput::Collinear));
        }

        normal = normal.normalize();

        // Recompute the normal to make sure it is perpendicular to unit_10.
        normal = (normal - normal.dot(unit_01) * unit_01).normalize();

        let mut max_distance = 0.0;
        let d0 = points[indices[2]].dot(normal);

        for i in 0..points.len() {
            let distance = (points[i].dot(normal) - d0).abs();

            if distance > max_distance
                && points[i] != points[indices[0]]
                && points[i] != points[indices[1]]
                && points[i] != points[indices[2]]
            {
                max_distance = distance;
                indices[3] = i;
            }
        }

        if max_distance.abs() <= 100.0 * threshold {
            return Err(ErrorKind::DegenerateInput(DegenerateInput::Coplanar));
        }

        if indices.len() != 4 {
            return Err(ErrorKind::RoundOffError(
                "number of simplex's vertices should be 4".to_string(),
            ));
        }

        Ok(indices)
    }
}

// get visible face viewed from furthest point
fn initialize_visible_set(
    points: &[DVec3],
    furthest_point_index: usize,
    faces: &BTreeMap<usize, Face>,
    faset_key: usize,
    face: &Face,
    threshold: f64,
) -> HashSet<usize> {
    let mut visible_set = HashSet::new();
    visible_set.insert(faset_key);
    let mut neighbor_stack: Vec<_> = face.neighbor_faces.to_vec();
    let mut visited_neighbor = HashSet::new();
    while let Some(neighbor_key) = neighbor_stack.pop() {
        if visited_neighbor.contains(&neighbor_key) {
            continue;
        }

        visited_neighbor.insert(neighbor_key);

        let neighbor = faces.get(&neighbor_key).unwrap();
        let pos = position_from_face(points, neighbor, furthest_point_index);
        if pos > threshold {
            visible_set.insert(neighbor_key);
            neighbor_stack.append(&mut neighbor.neighbor_faces.to_vec());
        }
    }
    visible_set
}

fn get_horizon(
    visible_set: &HashSet<usize>,
    faces: &BTreeMap<usize, Face>,
) -> Result<Vec<(Vec<usize>, usize)>, ErrorKind> {
    let mut horizon = Vec::new();
    for visible_key in visible_set {
        let visible_face = faces.get(visible_key).unwrap();
        let points_of_visible_face: HashSet<_> = visible_face.indices.iter().copied().collect();
        if points_of_visible_face.len() != 3 {
            return Err(ErrorKind::RoundOffError(
                "number of visible face's vertices should be 3".to_string(),
            ));
        }

        for neighbor_key in &visible_face.neighbor_faces {
            // if neighbor is unvisible
            if !visible_set.contains(neighbor_key) {
                let unvisible_neighbor = faces.get(neighbor_key).unwrap();
                let points_of_unvisible_neighbor: HashSet<_> =
                    unvisible_neighbor.indices.iter().copied().collect();
                if points_of_unvisible_neighbor.len() != 3 {
                    return Err(ErrorKind::RoundOffError(
                        "number of unvisible face's vertices should be 3".to_string(),
                    ));
                }

                let horizon_ridge: Vec<_> = points_of_unvisible_neighbor
                    .intersection(&points_of_visible_face)
                    .copied()
                    .collect();
                if horizon_ridge.len() != 2 {
                    return Err(ErrorKind::RoundOffError(
                        "number of ridge's vertices should be 2".to_string(),
                    ));
                }
                horizon.push((horizon_ridge, *neighbor_key));
            }
        }
    }
    if horizon.len() < 3 {
        return Err(ErrorKind::RoundOffError("horizon len < 3".to_string()));
    }
    Ok(horizon)
}

fn position_from_face(points: &[DVec3], face: &Face, point_index: usize) -> f64 {
    // TODO: Use `robust`
    let origin = face.distance_from_origin;
    let pos = face.normal.dot(points[point_index]);
    pos - origin
}

fn is_degenerate(points: &[DVec3], threshold: f64) -> bool {
    let ex_vec: Vec<DVec4> = points.iter().map(|v| v.extend(1.0)).collect();
    let num = ex_vec.len();
    if num <= 3 {
        return true;
    }
    let mut mat = DMat4::ZERO;
    for i in 0..4 {
        let mut row = DVec4::ZERO;
        for j in 0..4 {
            let mut c = 0.0;
            for vec in ex_vec.iter() {
                c += vec[i] * vec[j];
            }
            row[j] = c;
        }
        mat.x_axis[i] = row.x;
        mat.y_axis[i] = row.y;
        mat.z_axis[i] = row.z;
        mat.w_axis[i] = row.w;
    }
    mat.determinant() <= threshold && mat.determinant() >= -threshold
}

fn face_normal([a, b, c]: [DVec3; 3]) -> DVec3 {
    let ab = b - a;
    let ac = c - a;
    ab.cross(ac)
}

#[test]
fn face_normal_test() {
    let p1 = DVec3::new(-1.0, 0.0, 0.0);
    let p2 = DVec3::new(1.0, 0.0, 0.0);
    let p3 = DVec3::new(0.0, 1.0, 0.0);
    let normal_z = face_normal([p1, p2, p3]);
    assert_eq!(normal_z, DVec3::new(0.0, 0.0, 2.0));

    let p1 = DVec3::new(0.0, -1.0, 0.0);
    let p2 = DVec3::new(0.0, 1.0, 0.0);
    let p3 = DVec3::new(0.0, 0.0, 1.0);
    let normal_x = face_normal([p1, p2, p3]);
    assert_eq!(normal_x, DVec3::new(2.0, 0.0, 0.0));

    let p1 = DVec3::new(0.0, 0.0, -1.0);
    let p2 = DVec3::new(0.0, 0.0, 1.0);
    let p3 = DVec3::new(1.0, 0.0, 0.0);
    let normal_y = face_normal([p1, p2, p3]);
    assert_eq!(normal_y, DVec3::new(0.0, 2.0, 0.0));
}

#[test]
fn inner_outer_test() {
    let p1 = DVec3::new(1.0, 0.0, 0.0);
    let p2 = DVec3::new(0.0, 1.0, 0.0);
    let p3 = DVec3::new(0.0, 0.0, 1.0);
    let outer_point = DVec3::new(0.0, 0.0, 10.0);
    let inner_point = DVec3::new(0.0, 0.0, 0.0);
    let whithin_point = DVec3::new(1.0, 0.0, 0.0);
    let points = vec![p1, p2, p3, outer_point, inner_point, whithin_point];
    let face = Face::new(&points, [0, 1, 2]);
    let outer = position_from_face(&points, &face, 3);
    assert!(outer > 0.0);
    let inner = position_from_face(&points, &face, 4);
    assert!(inner < 0.0);
    let within = position_from_face(&points, &face, 5);
    assert!(within == 0.0);
}

#[test]
fn octahedron_test() {
    let p1 = DVec3::new(1.0, 0.0, 0.0);
    let p2 = DVec3::new(0.0, 1.0, 0.0);
    let p3 = DVec3::new(0.0, 0.0, 1.0);
    let p4 = DVec3::new(-1.0, 0.0, 0.0);
    let p5 = DVec3::new(0.0, -1.0, 0.0);
    let p6 = DVec3::new(0.0, 0.0, -1.0);
    let (_v, i) = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6], 0.001, None)
        .unwrap()
        .vertices_indices();
    assert_eq!(i.len(), 8 * 3);
}

#[test]
fn octahedron_translation_test() {
    let p1 = DVec3::new(1.0, 0.0, 0.0);
    let p2 = DVec3::new(0.0, 1.0, 0.0);
    let p3 = DVec3::new(0.0, 0.0, 1.0);
    let p4 = DVec3::new(-1.0, 0.0, 0.0);
    let p5 = DVec3::new(0.0, -1.0, 0.0);
    let p6 = DVec3::new(0.0, 0.0, -1.0);
    let points: Vec<_> = [p1, p2, p3, p4, p5, p6]
        .into_iter()
        .map(|p| p + DVec3::splat(10.0))
        .collect();
    let (_v, i) = ConvexHull::try_new(&points, 0.001, None)
        .unwrap()
        .vertices_indices();
    assert_eq!(i.len(), 8 * 3);
}

#[test]
fn cube_test() {
    let p1 = DVec3::new(1.0, 1.0, 1.0);
    let p2 = DVec3::new(1.0, 1.0, -1.0);
    let p3 = DVec3::new(1.0, -1.0, 1.0);
    let p4 = DVec3::new(1.0, -1.0, -1.0);
    let p5 = DVec3::new(-1.0, 1.0, 1.0);
    let p6 = DVec3::new(-1.0, 1.0, -1.0);
    let p7 = DVec3::new(-1.0, -1.0, 1.0);
    let p8 = DVec3::new(-1.0, -1.0, -1.0);
    let (_v, i) = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], 0.001, None)
        .unwrap()
        .vertices_indices();
    assert_eq!(i.len(), 6 * 2 * 3);
}

#[test]
fn cube_volume_test() {
    let p1 = DVec3::new(2.0, 2.0, 2.0);
    let p2 = DVec3::new(2.0, 2.0, 0.0);
    let p3 = DVec3::new(2.0, 0.0, 2.0);
    let p4 = DVec3::new(2.0, 0.0, 0.0);
    let p5 = DVec3::new(0.0, 2.0, 2.0);
    let p6 = DVec3::new(0.0, 2.0, 0.0);
    let p7 = DVec3::new(0.0, 0.0, 2.0);
    let p8 = DVec3::new(0.0, 0.0, 0.0);
    let cube = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], 0.001, None).unwrap();
    assert_eq!(cube.volume(), 8.0);
}

#[test]
fn cube_support_point_test() {
    let p1 = DVec3::new(1.0, 1.0, 1.0);
    let p2 = DVec3::new(1.0, 1.0, 0.0);
    let p3 = DVec3::new(1.0, 0.0, 1.0);
    let p4 = DVec3::new(1.0, 0.0, 0.0);
    let p5 = DVec3::new(0.0, 1.0, 1.0);
    let p6 = DVec3::new(0.0, 1.0, 0.0);
    let p7 = DVec3::new(0.0, 0.0, 1.0);
    let p8 = DVec3::new(0.0, 0.0, 0.0);
    let cube = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], 0.001, None).unwrap();
    assert_eq!(cube.support_point(DVec3::splat(0.5)), p1);
}

#[test]
#[should_panic(expected = "Degenerated")]
fn flat_test() {
    let p1 = DVec3::new(1.0, 1.0, 10.0);
    let p2 = DVec3::new(1.0, 1.0, 10.0);
    let p3 = DVec3::new(1.0, -1.0, 10.0);
    let p4 = DVec3::new(1.0, -1.0, 10.0);
    let p5 = DVec3::new(-1.0, 1.0, 10.0);
    let p6 = DVec3::new(-1.0, 1.0, 10.0);
    let p7 = DVec3::new(-1.0, -1.0, 10.0);
    let p8 = DVec3::new(-1.0, -1.0, 10.0);
    let (_v, _i) = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], 0.001, None)
        .unwrap()
        .vertices_indices();
}

#[test]
fn simplex_may_degenerate_test() {
    let points = vec![
        DVec3::new(1.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 1.0),
        DVec3::new(2.0, 1.0, 0.0),
        DVec3::new(2.0, 1.0, 1.0),
        DVec3::new(2.0, 0.0, 1.0),
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(1.0, 1.0, 2.0),
        DVec3::new(0.0, 1.0, 2.0),
        DVec3::new(0.0, 0.0, 2.0),
        DVec3::new(1.0, 0.0, 2.0),
    ];
    let (_v, _i) = ConvexHull::try_new(&points, 0.001, None)
        .unwrap()
        .vertices_indices();
}

#[test]
fn simplex_may_degenerate_test_2() {
    let vertices = vec![
        DVec3::new(0., 0., 0.),
        DVec3::new(1., 0., 0.),
        DVec3::new(1., 0., 1.),
        DVec3::new(0., 0., 1.),
        DVec3::new(0., 1., 0.),
        DVec3::new(1., 1., 0.),
        DVec3::new(1., 1., 1.),
        DVec3::new(0., 1., 1.),
        DVec3::new(2., 1., 0.),
        DVec3::new(2., 1., 1.),
        DVec3::new(2., 0., 1.),
        DVec3::new(2., 0., 0.),
        DVec3::new(1., 1., 2.),
        DVec3::new(0., 1., 2.),
        DVec3::new(0., 0., 2.),
        DVec3::new(1., 0., 2.),
    ];
    let indices = [4, 5, 1, 11, 1, 5, 1, 11, 10, 10, 2, 1, 5, 8, 11];
    let points = indices.iter().map(|i| vertices[*i]).collect::<Vec<_>>();
    let (_v, _i) = ConvexHull::try_new(&points, 0.001, None)
        .unwrap()
        .vertices_indices();
}

#[test]
fn sphere_test() {
    fn rot_z(point: DVec3, angle: f64) -> DVec3 {
        let e1 = angle.cos() * point[0] - angle.sin() * point[1];
        let e2 = angle.sin() * point[0] + angle.cos() * point[1];
        let e3 = point[2];
        DVec3::new(e1, e2, e3)
    }
    fn rot_x(point: DVec3, angle: f64) -> DVec3 {
        let e1 = point[0];
        let e2 = angle.cos() * point[1] - angle.sin() * point[2];
        let e3 = angle.sin() * point[1] + angle.cos() * point[2];
        DVec3::new(e1, e2, e3)
    }
    let mut points = Vec::new();
    let dev = 10;
    let unit_y = DVec3::Y;
    for step_x in 0..dev {
        let angle_x = 2.0 * std::f64::consts::PI * (step_x as f64 / dev as f64);
        let p = rot_x(unit_y, angle_x);
        for step_z in 0..dev {
            let angle_z = 2.0 * std::f64::consts::PI * (step_z as f64 / dev as f64);
            let p = rot_z(p, angle_z);
            points.push(p);
        }
    }
    let (_v, _i) = ConvexHull::try_new(&points, 0.001, None)
        .unwrap()
        .vertices_indices();
}

#[test]
#[ignore]
fn heavy_sphere_test() {
    fn rot_z(point: DVec3, angle: f64) -> DVec3 {
        let e1 = angle.cos() * point[0] - angle.sin() * point[1];
        let e2 = angle.sin() * point[0] + angle.cos() * point[1];
        let e3 = point[2];
        DVec3::new(e1, e2, e3)
    }
    fn rot_x(point: DVec3, angle: f64) -> DVec3 {
        let e1 = point[0];
        let e2 = angle.cos() * point[1] - angle.sin() * point[2];
        let e3 = angle.sin() * point[1] + angle.cos() * point[2];
        DVec3::new(e1, e2, e3)
    }
    let mut points = Vec::new();
    let dev = 100;
    let unit_y = DVec3::Y;
    for step_x in 0..dev {
        let angle_x = 2.0 * std::f64::consts::PI * (step_x as f64 / dev as f64);
        let p = rot_x(unit_y, angle_x);
        for step_z in 0..dev {
            let angle_z = 2.0 * std::f64::consts::PI * (step_z as f64 / dev as f64);
            let p = rot_z(p, angle_z);
            points.push(p);
        }
    }
    let (_v, _i) = ConvexHull::try_new(&points, 0.001, None)
        .unwrap()
        .vertices_indices();
}

#[test]
fn is_degenerate_test() {
    let points = vec![
        DVec3::new(1., 0., 0.),
        DVec3::new(0., 0., 0.),
        DVec3::new(0., 1., 0.),
        DVec3::new(1., 0., 1.),
    ];
    assert!(!is_degenerate(&points, 0.00001));
}