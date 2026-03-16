---
name: Implement Unimplemented SLICOT Functions
overview: Implement the 9 remaining SLICOT routines (7 SG03 + 2 TB03) with full pure-Rust implementations, following the existing sg03/tb03 dependency order and atomic single-file commits (new .rs + mod.rs update per routine).
todos: []
isProject: false
---

# Implement All Unimplemented SLICOT Functions

## Unimplemented SLICOT functions (authoritative list)

Source: [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md) — all rows with status **not started**. The mapping lists **15** routines; **6** already have full implementations in the repo (sg03br, sg03bu, sg03bv, sg03bw, sg03bx, sg03by). The **9** that still need implementation:


| SLICOT | Rust module | Rust function | Fortran source |
| ------ | ----------- | ------------- | -------------- |
| SG03AD | sg03        | sg03ad        | SG03AD.f       |
| SG03AX | sg03        | sg03ax        | SG03AX.f       |
| SG03AY | sg03        | sg03ay        | SG03AY.f       |
| SG03BD | sg03        | sg03bd        | SG03BD.f       |
| SG03BS | sg03        | sg03bs        | SG03BS.f       |
| SG03BT | sg03        | sg03bt        | SG03BT.f       |
| SG03BZ | sg03        | sg03bz        | SG03BZ.f       |
| TB03AD | tb03        | tb03ad        | TB03AD.f       |
| TB03AY | tb03        | tb03ay        | TB03AY.f       |


Doc links and paths are in the [Unimplemented list table](docs/SLICOT_MAPPING.md) (lines 635–654).

## Current state

- **sg03**: [src/sg03/mod.rs](src/sg03/mod.rs) declares 6 modules (br, bw, bx, by, bu, bv); all six are full implementations. Missing: ad, ax, ay, bd, bs, bt, bz.
- **tb03**: [src/tb03/mod.rs](src/tb03/mod.rs) exists but has no submodules; both tb03ad and tb03ay need to be added.

## Implementation order (dependencies)

**SG03 (bottom-up):**

1. **SG03BS, SG03BT** — Cholesky-factor solvers (discrete / continuous); use existing SG03BX, SG03BW.
2. **SG03AX, SG03AY** — Discrete / continuous Lyapunov solvers.
3. **SG03AD, SG03BD, SG03BZ** — Driver and Cholesky drivers (top-level entry points).

**TB03:**

1. **TB03AY** first — builds polynomial matrix V(s) block-by-block (~169 lines Fortran).
2. **TB03AD** second — calls AB07MD, TB01ID, TB01UD, TC01OD, TB03AY, TB01YD, MA02GD, etc. (~731 lines Fortran).

## File and commit strategy

- **New files**: One `.rs` file per routine under `src/sg03/` or `src/tb03/`.
- **mod.rs**: Each new routine requires one new `pub mod <name>;` in the corresponding mod.rs.
- **Atomic commits** (Exception 1b): For each new routine, commit **together** the new `.rs` file and the `mod.rs` update.
- **After all implementations**: Update SLICOT_MAPPING.md (set the 9 routines to `done`); optionally FEATURES.md.

## Implementation requirements (no stubs)

- **Pure Rust**: No `extern`, FFI, or Fortran linkage; use nalgebra for linear algebra.
- **1:1 mapping**: Same semantics and INFO codes as the SLICOT Fortran routine.
- **Tests**: Each new file must contain `#[cfg(test)] mod tests` with at least one meaningful test.
