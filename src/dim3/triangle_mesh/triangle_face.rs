use approx::relative_eq;
use glam::Vec3A;

/// The index of a point in the input point set.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PointId(pub u32);

impl PointId {
    /// A placeholder `PointId` that does not correspond to any valid point.
    pub const PLACEHOLDER: PointId = PointId(u32::MAX);

    /// Returns the underlying index of the point as a `usize`.
    #[inline]
    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<u32> for PointId {
    #[inline]
    fn from(value: u32) -> Self {
        PointId(value)
    }
}

/// The index of a [`TriangleFace`] in a [`ConvexHull3d`](crate::dim3::ConvexHull3d).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FaceId(pub u32);

impl FaceId {
    /// A placeholder `FaceId` that does not correspond to any valid face.
    pub const PLACEHOLDER: FaceId = FaceId(u32::MAX);

    /// Returns the underlying index of the face as a `usize`.
    #[inline]
    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<u32> for FaceId {
    #[inline]
    fn from(value: u32) -> Self {
        FaceId(value)
    }
}

/// The index of an edge in a [`TriangleFace`], between 0 and 2.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct EdgeIndex(pub u32);

impl EdgeIndex {
    /// A placeholder `EdgeIndex` that does not correspond to any valid edge.
    pub const PLACEHOLDER: EdgeIndex = EdgeIndex(u32::MAX);

    /// Returns the underlying index of the edge as a `usize`.
    #[inline]
    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<u32> for EdgeIndex {
    #[inline]
    fn from(value: u32) -> Self {
        EdgeIndex(value)
    }
}

/// A handle to a face, consisting of a [`FaceId`] and the [`EdgeIndex`]
/// of an edge that connects the face to its neighbor.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FaceHandle {
    /// The ID of the face.
    pub face: FaceId,
    /// The index of the edge connecting the face to its neighbor.
    pub edge: EdgeIndex,
}

impl FaceHandle {
    /// A placeholder `FaceNeighbor` that does not correspond to any valid neighbor.
    pub const PLACEHOLDER: FaceHandle = FaceHandle {
        face: FaceId::PLACEHOLDER,
        edge: EdgeIndex::PLACEHOLDER,
    };

    /// Creates a new [`FaceNeighbor`] with the given face ID and edge index.
    #[inline]
    pub fn new(face: impl Into<FaceId>, edge: impl Into<EdgeIndex>) -> Self {
        Self {
            face: face.into(),
            edge: edge.into(),
        }
    }
}

/// A triangular face belonging to a [`ConvexHull3d`](crate::dim3::ConvexHull3d).
#[derive(Clone, Debug, PartialEq)]
pub struct TriangleFace {
    /// Whether the face is valid (not removed).
    pub(crate) valid: bool,
    /// The points that make up the face.
    pub(crate) points: [PointId; 3],
    /// The neighboring faces across each edge of the face.
    ///
    /// Each entry contains the [`FaceHandle`] of the neighboring face,
    /// consisting of the face ID and the edge index that connects back to this face.
    pub(crate) neighbors: [FaceHandle; 3],
    /// The points in front of the face plane, or the points that can "see" the face.
    pub(crate) outside_points: Vec<PointId>,
    /// The normal of the face.
    pub(crate) normal: Vec3A,
    /// Whether the points of the face are affinely dependent,
    /// meaning they lie on a single line or point.
    pub(crate) affinely_dependent: bool,
    /// The identifier and distance of the furthest outside point, if any.
    pub(crate) furthest_outside_point: Option<(PointId, f32)>,
}

impl TriangleFace {
    /// Creates a [`Face`] using the `points` with the given `indices`.
    #[inline]
    pub fn from_triangle(points: &[Vec3A], indices: [PointId; 3]) -> Self {
        const EPSILON_SQ: f32 = (f32::EPSILON * 100.0) * (f32::EPSILON * 100.0);

        let [a, b, c] = indices.map(|id| points[id.index()]);
        let ab = b - a;
        let ac = c - a;
        let scaled_normal = ab.cross(ac);

        // Check for affinely dependent points (zero-area triangle).
        let affinely_dependent =
            relative_eq!(scaled_normal.length_squared(), 0.0, epsilon = EPSILON_SQ);

        Self {
            points: indices,
            valid: true,
            affinely_dependent,
            neighbors: [FaceHandle::PLACEHOLDER; 3],
            outside_points: Vec::new(),
            normal: scaled_normal.normalize_or_zero(),
            furthest_outside_point: None,
        }
    }

    /// Returns the neighboring face across the edge at the given index.
    #[inline]
    pub fn neighbor(&self, edge_index: EdgeIndex) -> FaceHandle {
        self.neighbors[edge_index.index()]
    }

    /// Returns a mutable reference to the neighboring face across the edge at the given index.
    #[inline]
    pub fn neighbor_mut(&mut self, edge_index: EdgeIndex) -> &mut FaceHandle {
        &mut self.neighbors[edge_index.index()]
    }

    /// Sets the neighboring faces of the face.
    #[inline]
    pub fn set_neighbors(
        &mut self,
        face1: impl Into<FaceId>,
        face2: impl Into<FaceId>,
        face3: impl Into<FaceId>,
        edge1: impl Into<EdgeIndex>,
        edge2: impl Into<EdgeIndex>,
        edge3: impl Into<EdgeIndex>,
    ) {
        self.neighbors[0] = FaceHandle::new(face1, edge1);
        self.neighbors[1] = FaceHandle::new(face2, edge2);
        self.neighbors[2] = FaceHandle::new(face3, edge3);
    }

    /// Returns the first point index of the edge at the given index.
    #[inline]
    pub fn first_point_from_edge(&self, edge_index: EdgeIndex) -> PointId {
        self.points[edge_index.index()]
    }

    /// Returns the second point index of the edge at the given index.
    #[inline]
    pub fn second_point_from_edge(&self, edge_index: EdgeIndex) -> PointId {
        self.points[(edge_index.index() + 1) % 3]
    }

    /// Determines whether the face can be "seen" by the given point, in an order-independent manner,
    /// meaning that the result does not depend on the order of the face's points.
    #[inline]
    pub fn can_be_seen_by_point_order_independent(
        &self,
        point_id: PointId,
        points: &[Vec3A],
    ) -> bool {
        // If the face is degenerate, it can be seen by any point.
        if self.affinely_dependent {
            return true;
        }

        for i in 0..3 {
            let p0 = points[self.points[i].index()];
            let point = points[point_id.index()];
            let distance = (point - p0).dot(self.normal);
            if distance >= 0.0 {
                return true;
            }
        }

        false
    }

    /// Returns the distance to the given point.
    #[inline]
    pub fn distance_to_point(&self, point_id: PointId, points: &[Vec3A]) -> f32 {
        let p0 = points[self.points[0].index()];
        let point = points[point_id.index()];
        (point - p0).dot(self.normal)
    }

    /// Returns the distance to the given point if the face can be "seen" from it.
    #[inline]
    pub fn distance_to_visible_point(&self, point_id: PointId, points: &[Vec3A]) -> Option<f32> {
        // If the face is degenerate, it cannot see any points.
        if self.affinely_dependent {
            return None;
        }

        let distance = self.distance_to_point(point_id, points);
        (distance >= f32::EPSILON * 100.0).then_some(distance)
    }

    /// Adds the given point to the list of outside points for the face.
    #[inline]
    pub fn add_outside_point(&mut self, point_id: PointId, points: &[Vec3A]) {
        let distance = self.distance_to_point(point_id, points);
        debug_assert!(distance >= f32::EPSILON);

        if self
            .furthest_outside_point
            .is_none_or(|(_, d)| distance > d)
        {
            self.furthest_outside_point = Some((point_id, distance));
        }

        self.outside_points.push(point_id);
    }

    /// Tries to add the given point to the list of outside points for the face.
    #[inline]
    pub fn try_add_outside_point(&mut self, point_id: PointId, points: &[Vec3A]) {
        if let Some(distance) = self.distance_to_visible_point(point_id, points) {
            if self
                .furthest_outside_point
                .is_none_or(|(_, d)| distance > d)
            {
                self.furthest_outside_point = Some((point_id, distance));
            }

            self.outside_points.push(point_id);
        }
    }
}
