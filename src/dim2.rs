// Adapted from geo:
// https://github.com/georust/geo/blob/8940db79aa6aa4ec8820d6328c68a2ae08ac8fdc/geo/src/algorithm/convex_hull/qhull.rs

use glam::Vec2;

/// A 2D [convex hull] representing the smallest convex set containing
/// all input points in a given point set.
///
/// This can be thought of as a shrink wrapping of a 2D object.
///
/// [convex hull]: https://en.wikipedia.org/wiki/Convex_hull
///
/// # Example
///
/// ```
/// use glam::Vec2;
/// use quickhull::ConvexHull2d;
///
/// let points = vec![
///     Vec2::new(0.0, 0.0),
///     Vec2::new(1.0, 0.0),
///     Vec2::new(0.0, 1.0),
///     Vec2::new(1.0, 1.0),
///     Vec2::new(0.5, 0.5),
/// ];
///
/// // Create the convex hull.
/// let hull = ConvexHull2d::from_points(&points);
///
/// // Get the points of the convex hull in counterclockwise order.
/// let points = hull.into_vertices();
///
/// assert_eq!(
///     points,
///     &[
///         // The points may not be in the same order as the input.
///         Vec2::new(1.0, 0.0),
///         Vec2::new(1.0, 1.0),
///         Vec2::new(0.0, 1.0),
///         Vec2::new(0.0, 0.0),
///     ],
/// );
/// ```
pub struct ConvexHull2d {
    /// The vertices of the convex hull in counterclockwise order.
    vertices: Vec<Vec2>,
}

impl ConvexHull2d {
    /// Computes a [`ConvexHull2d`] for the given set of 2D points.
    ///
    /// This allocates a new vector for the points. To avoid this allocation,
    /// consider using [`from_mut_points`](Self::from_mut_points).
    #[inline]
    pub fn from_points(points: &[Vec2]) -> Self {
        Self::from_mut_points(&mut points.to_vec())
    }

    /// Computes a [`ConvexHull2d`] for the given mutable set of 2D points.
    ///
    /// The input slice may be reordered and truncated during hull construction.
    /// This is more efficient than [`from_points`](Self::from_points), as it avoids an allocation.
    #[inline]
    pub fn from_mut_points(mut points: &mut [Vec2]) -> Self {
        if points.len() <= 3 {
            return Self::trivial_hull(points.to_vec());
        }

        let mut hull = Vec::new();

        // Find the points with minimum and maximum `x` coordinates.
        let (min, max) = {
            let (min_index, mut max_index) =
                points
                    .iter()
                    .enumerate()
                    .fold((0, 0), |(min_index, max_index), (i, point)| {
                        let min_ordering = lexicographic_cmp(point, &points[min_index]);
                        let min_index = if min_ordering == std::cmp::Ordering::Less {
                            i
                        } else {
                            min_index
                        };

                        let max_ordering = lexicographic_cmp(point, &points[max_index]);
                        let max_index = if max_ordering == std::cmp::Ordering::Greater {
                            i
                        } else {
                            max_index
                        };

                        (min_index, max_index)
                    });

            // Move the min and max points from `points` to the hull.
            let min = *swap_with_first_and_remove(&mut points, min_index);

            // If the max point was at index 0, it was just swapped to `min_index`.
            if max_index == 0 {
                max_index = min_index;
            }

            // If `min_index == max_index`, then an point could be chosen as the maximum.
            // But based on the check above, it could now be 0, in which case we should not decrement it.
            max_index = max_index.saturating_sub(1);

            // Move the max point from `points` to the hull.
            let max = *swap_with_first_and_remove(&mut points, max_index);

            (min, max)
        };

        // Recursively find hull points on either side of the line segment `min, max`.
        {
            let (points, _) = partition_slice(points, |point| is_ccw(max, min, *point));
            Self::hull_set(max, min, points, &mut hull);
        }
        hull.push(max);
        let (points, _) = partition_slice(points, |point| is_ccw(min, max, *point));
        Self::hull_set(min, max, points, &mut hull);
        hull.push(min);

        Self { vertices: hull }
    }

    /// Computes the indices of the convex hull vertices within the original point set.
    ///
    /// The returned indices are in counterclockwise order.
    #[inline]
    pub fn indices_from_points(points: &[Vec2]) -> Vec<usize> {
        // TODO: Optimize this.
        // TODO: Handle duplicate points.
        let hull = Self::from_points(points);
        hull.vertices()
            .iter()
            .map(|hull_point| {
                points
                    .iter()
                    .position(|point| *point == *hull_point)
                    .unwrap()
            })
            .collect()
    }

    /// Returns the vertices of the convex hull in counterclockwise order.
    #[inline]
    pub fn vertices(&self) -> &[Vec2] {
        &self.vertices
    }

    /// Returns the vertices of the convex hull in counterclockwise order.
    ///
    /// This consumes the hull and allows taking ownership of the underlying data without cloning.
    #[inline]
    pub fn into_vertices(self) -> Vec<Vec2> {
        self.vertices
    }

    /// Constructs the convex hull for a point set of size 3 or less.
    ///
    /// The points are sorted lexicographically and oriented counterclockwise.
    /// If two points are collinear, the middle point is removed.
    ///
    /// # Panics
    ///
    /// Panics with `debug_assertions` enabled if `points.len() > 3`.
    fn trivial_hull(mut points: Vec<Vec2>) -> Self {
        debug_assert!(points.len() <= 3);

        points.sort_unstable_by(lexicographic_cmp);

        let orientation = orient2d(points[0], points[1], points[2]);

        // Remove the middle point if all points are collinear.
        if points.len() == 3 && orientation == 0.0 {
            points.remove(1);
        }

        // Ensure CCW orientation.
        if orientation < 0.0 {
            points.swap(1, 2);
        }

        Self { vertices: points }
    }

    // Recursively computes the convex hull of a subset of points.
    fn hull_set(a: Vec2, b: Vec2, mut points: &mut [Vec2], hull: &mut Vec<Vec2>) {
        if points.is_empty() {
            return;
        }

        if points.len() == 1 {
            hull.push(points[0]);
            return;
        }

        let ab = b - a;
        let orthogonal = ab.perp();

        // Find the point furthest from the line segment `ab`.
        let furthest_index = points
            .iter()
            .map(|point| orthogonal.dot(*point - a))
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;

        // Move the furthest point from `points` to the hull.
        let furthest_point = *swap_with_first_and_remove(&mut points, furthest_index);

        // Recursively find hull points on either side of the triangle `a, b, furthest_point`.
        {
            let (points, _) = partition_slice(points, |point| is_ccw(furthest_point, b, *point));
            Self::hull_set(furthest_point, b, points, hull);
        }
        hull.push(furthest_point);
        let (points, _) = partition_slice(points, |point| is_ccw(a, furthest_point, *point));
        Self::hull_set(a, furthest_point, points, hull);
    }
}

/// Gives the orientation of the triangle formed by `a`, `b`, `c`.
///
/// - `orientation > 0`: counterclockwise
/// - `orientation < 0`: clockwise
/// - `orientation == 0`: collinear
#[inline]
fn orient2d(a: Vec2, b: Vec2, c: Vec2) -> f64 {
    use robust::Coord;
    robust::orient2d(
        Coord { x: a.x, y: a.y },
        Coord { x: b.x, y: b.y },
        Coord { x: c.x, y: c.y },
    )
}

/// Returns `true` if the points `a`, `b`, `c` are oriented counterclockwise.
#[inline]
fn is_ccw(a: Vec2, b: Vec2, c: Vec2) -> bool {
    orient2d(a, b, c) > 0.0
}

/// Compares two 2D points first by `x`, then by `y`.
#[inline]
fn lexicographic_cmp(a: &Vec2, b: &Vec2) -> std::cmp::Ordering {
    a.x.partial_cmp(&b.x)
        .unwrap()
        .then(a.y.partial_cmp(&b.y).unwrap())
}

/// Partitions a mutable slice in-place so that it contains all elements for
/// which `predicate(e)` is `true`, followed by all elements for which
/// `predicate(e)` is `false`. Returns sub-slices to all predicated and
/// non-predicated elements, respectively.
///
/// https://github.com/llogiq/partition/blob/master/src/lib.rs
fn partition_slice<T, P>(data: &mut [T], predicate: P) -> (&mut [T], &mut [T])
where
    P: Fn(&T) -> bool,
{
    let len = data.len();

    if len == 0 {
        return (&mut [], &mut []);
    }

    let (mut left, mut right) = (0, len - 1);

    loop {
        while left < len && predicate(&data[left]) {
            left += 1;
        }

        while right > 0 && !predicate(&data[right]) {
            right -= 1;
        }

        if left >= right {
            return data.split_at_mut(left);
        }

        data.swap(left, right);
    }
}

/// Swaps the element at `index` with the first element of `slice`, removes it from the slice,
/// and returns a mutable reference to it.
#[inline]
fn swap_with_first_and_remove<'a, T>(slice: &mut &'a mut [T], index: usize) -> &'a mut T {
    // Temporarily replace `slice` with an empty value.
    let tmp = std::mem::take(slice);
    tmp.swap(0, index);
    let (h, t) = tmp.split_first_mut().unwrap();
    *slice = t;
    h
}

#[cfg(test)]
mod test {
    use glam::vec2;

    use super::*;

    #[test]
    fn hull_correct() {
        let points = vec![
            vec2(0.0, 10.0),
            vec2(1.0, 1.0),
            vec2(10.0, 0.0),
            vec2(1.0, -1.0),
            vec2(0.0, -10.0),
            vec2(-1.0, -1.0),
            vec2(-10.0, 0.0),
            vec2(-1.0, 1.0),
            vec2(0.0, 10.0),
        ];
        let expected = vec![
            vec2(0.0, -10.0),
            vec2(10.0, 0.0),
            vec2(0.0, 10.0),
            vec2(-10.0, 0.0),
        ];
        let result = ConvexHull2d::from_points(&points).into_vertices();
        assert_eq!(result, expected);
    }

    #[test]
    fn ccw() {
        let points = vec![
            vec2(1.0, 0.0),
            vec2(2.0, 1.0),
            vec2(1.75, 1.1),
            vec2(1.0, 2.0),
            vec2(0.0, 1.0),
            vec2(1.0, 0.0),
        ];
        let expected = [
            vec2(1.0, 0.0),
            vec2(2.0, 1.0),
            vec2(1.0, 2.0),
            vec2(0.0, 1.0),
        ];
        let result = ConvexHull2d::from_points(&points).into_vertices();
        assert_eq!(result, expected);
    }

    #[test]
    fn trivial_collinear() {
        let points = vec![vec2(0.0, 0.0), vec2(1.0, 1.0), vec2(2.0, 2.0)];
        let expected = vec![vec2(0.0, 0.0), vec2(2.0, 2.0)];
        let result = ConvexHull2d::from_points(&points).into_vertices();
        assert_eq!(result, expected);
    }

    #[test]
    fn non_trivial_collinear() {
        let points = vec![
            vec2(0.0, 0.0),
            vec2(1.0, 1.0),
            vec2(2.0, 2.0),
            vec2(0.0, 2.0),
            vec2(2.0, 0.0),
        ];
        let expected = vec![
            vec2(2.0, 0.0),
            vec2(2.0, 2.0),
            vec2(0.0, 2.0),
            vec2(0.0, 0.0),
        ];
        let result = ConvexHull2d::from_points(&points).into_vertices();
        assert_eq!(result, expected);
    }
}
