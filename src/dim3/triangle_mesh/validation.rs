use super::{
    triangle_face::{EdgeIndex, TriangleFace},
    FaceId,
};

pub fn validate_face_connectivity(face_id: FaceId, faces: &[TriangleFace]) {
    let face = &faces[face_id.index()];

    for i in 0..3 {
        let edge = EdgeIndex(i as u32);
        let neighbor = face.neighbor(edge);
        let neighbor_face = &faces[neighbor.face.index()];

        assert!(neighbor_face.valid, "Neighbor face is not valid.");

        // Check that the neighbor face points back to this face.
        assert_eq!(
            neighbor_face.neighbor(face.neighbor(edge).edge).face,
            face_id,
        );
        assert_eq!(
            neighbor_face.neighbor(face.neighbor(edge).edge).edge,
            EdgeIndex(i as u32),
        );

        // Check that the points across the shared edge match.
        assert_eq!(
            neighbor_face.first_point_from_edge(face.neighbor(edge).edge),
            face.second_point_from_edge(EdgeIndex(i as u32)),
        );
        assert_eq!(
            neighbor_face.second_point_from_edge(face.neighbor(edge).edge),
            face.first_point_from_edge(EdgeIndex(i as u32)),
        );
    }
}

// TODO: Validate that the hull is properly formed and convex.
