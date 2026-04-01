use crate::dim3::initial_hull::InitialConvexHull3dError;
use thiserror::Error;

/// An error returned during 3D convex hull construction.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum ConvexHull3dError {
    /// Could not find a support point in the given direction.
    #[error("Input points are either invalid (NaN/Inf) or nearly coplanar.")]
    MissingSupportPoint,

    /// An error during initial hull construction.
    #[error("Initial hull construction failed.")]
    InitialHullError(#[from] InitialConvexHull3dError),

    /// An error in the algorithm itself. Please report is as a bug
    /// with a minimal reproducible example.
    #[error("Internal error: {0}")]
    InternalError(&'static str),
}
