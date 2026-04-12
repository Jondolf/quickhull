use crate::collections::FixedHasher;
use core::ops::{Deref, DerefMut};

/// A new-type for [`HashMap`](hashbrown::HashMap) with [`FixedHasher`] as the hashing provider.
/// Can be trivially converted to and from a [hashbrown] [`HashMap`](hashbrown::HashMap) using [`From`].
///
/// This hasher is not DoS-resistant and is not suitable for cryptographic purposes,
/// but it has high performance and determinism.
///
/// A new-type is used instead of a type alias due to critical methods like [`new`](hashbrown::HashMap::new)
/// being incompatible with a non-default hasher.
#[derive(Clone, Debug, Default)]
pub struct HashMap<K, V>(hashbrown::HashMap<K, V, FixedHasher>);

impl<K, V> Deref for HashMap<K, V> {
    type Target = hashbrown::HashMap<K, V, FixedHasher>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<K, V> DerefMut for HashMap<K, V> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<K, V, T> FromIterator<T> for HashMap<K, V>
where
    hashbrown::HashMap<K, V, FixedHasher>: FromIterator<T>,
{
    #[inline]
    fn from_iter<U: IntoIterator<Item = T>>(iter: U) -> Self {
        Self(FromIterator::from_iter(iter))
    }
}

impl<K, V> HashMap<K, V> {
    /// Creates an empty [`HashMap`].
    #[inline]
    pub const fn new() -> Self {
        Self(hashbrown::HashMap::with_hasher(FixedHasher))
    }
}
