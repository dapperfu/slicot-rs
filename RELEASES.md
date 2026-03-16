# Release notes

Short, human-readable highlights for each release.

---

## 0.4.20 (2025-03-16)

**First production-oriented release.** This crate is intended as a **drop-in replacement for SLICOT** in the Rust ecosystem: same routine names and semantics where implemented, pure Rust with [nalgebra](https://crates.io/crates/nalgebra), no LAPACK/BLAS FFI.

Highlights:

- Production-ready layout: LICENSE (BSD-3-Clause), CHANGELOG, runnable examples, CI (GitHub Actions), crates.io metadata (keywords, categories, MSRV 1.70).
- Scripts moved from `scripts/` to `tools/`; published package excludes `plans/`, `tools/`, `SLICOT-Reference/`, `fuzz/`, `.github/`.
- Examples for control/systems (AB01ND), linear algebra (MA02ED), and transforms (DE01OD).

The initial release version 0.4.20 indicates that only someone extremely intoxicated with cannabis would ever suggest such a stupid idea.
