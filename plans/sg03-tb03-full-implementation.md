---
name: SG03 TB03 full implementation
overview: "Fully implement 15 SLICOT routines: SG03 (13) and TB03 (2) in pure Rust with no stubs, following existing module patterns and bottom-up dependency order for SG03. Each routine gets a dedicated `.rs` file and tests (unit tests plus example-based validation where SLICOT .dat/.res exist)."
todos: []
isProject: false
---

# Full implementation: SG03 (13) and TB03 (2) with tests

## Scope

- **SG03** (13 routines): Generalized Lyapunov equations and helpers. Sources: [SLICOT-Reference/src/SG03*.f](SLICOT-Reference/src/). Doc links in [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md) (lines 639–651).
- **TB03** (2 routines): Left/right polynomial matrix representation from state-space. Sources: [TB03AD.f](SLICOT-Reference/src/TB03AD.f), [TB03AY.f](SLICOT-Reference/src/TB03AY.f).
- **No stubs**: Every routine is a full implementation (pure Rust, nalgebra; no FFI).
- **Tests**: Every routine has at least one unit test; SG03AD, SG03BD, TB03AD also have example-based tests using [SLICOT-Reference/examples/data/](SLICOT-Reference/examples/data/) and [results/](SLICOT-Reference/examples/results/).

## SG03 dependency order (bottom-up)

Implement in this order so callers exist before routines that use them:

| Phase | Routines | Notes |
|-------|----------|--------|
| 1 | SG03BR, SG03BY | Leaf: complex Givens/rotation in real arithmetic (~222 + ~89 lines Fortran) |
| 2 | SG03BX | 2×2 generalized Lyapunov (uses SG03BR); ~861 lines |
| 3 | SG03BW | Generalized Sylvester (quasi-triangular A); ~424 lines |
| 4 | SG03BS, SG03BT | Cholesky-factor solvers (discrete / continuous); ~648, ~577 lines |
| 5 | SG03BU, SG03BV | Block solvers (use SG03BX, SG03BW); ~692, ~632 lines |
| 6 | SG03AX, SG03AY | Discrete / continuous Lyapunov solvers; ~672, ~671 lines |
| 7 | SG03AD, SG03BD, SG03BZ | Driver + Cholesky drivers; ~663, ~1012, ~943 lines |

SG03 Fortran uses LAPACK (e.g. DGGES for generalized Schur, DLAMCH, DLABAD). Replace with nalgebra: generalized Schur and machine constants via f64::MIN_POSITIVE, f64::EPSILON, etc.

## TB03 dependency order

- **TB03AY** first (~169 lines): builds polynomial matrix V(s) block-by-block.
- **TB03AD** second (~731 lines): calls AB07MD, TB01ID, TB01UD, TC01OD, TB03AY, TB01YD, MA02GD, etc.

## File layout

- **New modules**: `src/sg03/` and `src/tb03/`.
- **sg03**: mod.rs plus sg03ad.rs, sg03ax.rs, sg03ay.rs, sg03bd.rs, sg03br.rs, sg03bs.rs, sg03bt.rs, sg03bu.rs, sg03bv.rs, sg03bw.rs, sg03bx.rs, sg03by.rs, sg03bz.rs.
- **tb03**: mod.rs, tb03ad.rs, tb03ay.rs.
- **lib.rs**: add `pub mod sg03;` and `pub mod tb03;`.

## Tests

- Per-routine: Each file contains `#[cfg(test)] mod tests` with at least one meaningful test.
- Example-based: SG03AD, SG03BD, TB03AD use SLICOT example data/results where available.

## Implementation summary

| Step | Action |
|------|--------|
| 1 | Save this plan under `plans/` (first TODO). |
| 2 | Add src/sg03/mod.rs and src/tb03/mod.rs; register in lib.rs. |
| 3 | Implement SG03BR, SG03BY (with tests). |
| 4 | Implement SG03BX, SG03BW (with tests). |
| 5 | Implement SG03BS, SG03BT (with tests). |
| 6 | Implement SG03BU, SG03BV (with tests). |
| 7 | Implement SG03AX, SG03AY (with tests). |
| 8 | Implement SG03AD, SG03BD, SG03BZ (with example-based tests for AD, BD). |
| 9 | Implement TB03AY then TB03AD (with example-based test for TB03AD). |
| 10 | Update SLICOT_MAPPING.md and FEATURES.md; run cargo test. |

All new code must be pure Rust (no extern/FFI). Follow existing style and workspace git/commit rules.
