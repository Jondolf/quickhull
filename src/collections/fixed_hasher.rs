use core::hash::BuildHasher;
use foldhash::fast::{FixedState, FoldHasher};

/// For when you want a deterministic hasher.
///
/// Seed was randomly generated with a fair dice roll. Guaranteed to be random:
/// <https://github.com/bevyengine/bevy/pull/1268/files#r560918426>
const FIXED_HASHER: FixedState =
    FixedState::with_seed(0b1001010111101110000001001100010000000011001001101011001001111000);

/// Deterministic hasher based upon a random but fixed state.
#[derive(Copy, Clone, Default, Debug)]
pub struct FixedHasher;
impl BuildHasher for FixedHasher {
    type Hasher = FoldHasher<'static>;

    #[inline]
    fn build_hasher(&self) -> Self::Hasher {
        FIXED_HASHER.build_hasher()
    }
}
