//! Collections optimized for performance and determinism.

mod fixed_hasher;
mod hash_map;

pub use fixed_hasher::FixedHasher;
pub use hash_map::HashMap;
