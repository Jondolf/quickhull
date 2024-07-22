# Quickhull

A Rust-implementation of the Quickhull algorithm for computing convex hulls for point sets.

This is a simplified and cleaned up version of [chull](https://github.com/u65xhd/chull),
focusing on making the algorithm robust and efficient for the 2D (note: not implemented yet!) and 3D cases.

## Warning ⚠️

This is a work-in-progress, and currently has some robustness issues. A crate will be released if/when
both 2D and 3D are properly supported and robustness is at an acceptable level.

## References

- C. Bradford Barber et al. 1996. [The Quickhull Algorithm for Convex Hulls](https://www.cise.ufl.edu/~ungor/courses/fall06/papers/QuickHull.pdf) (the original paper)
- Dirk Gregorius. GDC 2014. [Physics for Game Programmers: Implementing Quickhull](https://archive.org/details/GDC2014Gregorius)

## License

This Quickhull crate is free and open source. All code in this repository is dual-licensed under either:

- MIT License ([LICENSE-MIT](/LICENSE-MIT) or <http://opensource.org/licenses/MIT>)
- Apache License, Version 2.0 ([LICENSE-APACHE](/LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)

at your option.
