# Implement MC01, MC03, MD03, NF01

## Scope

- **MC01** (15 routines): MC01MD, MC01ND, MC01OD, MC01PD, MC01PY, MC01QD, MC01RD, MC01SD, MC01SW, MC01SX, MC01SY, MC01TD, MC01VD, MC01WD, MC01XD — polynomial operations (real/complex).
- **MC03** (4 routines): MC03MD, MC03ND, MC03NX, MC03NY — matrix pencil operations.
- **MD03** (7 routines): MD03AD, MD03BA, MD03BB, MD03BD, MD03BF, MD03BX, MD03BY — data fitting / nonlinear least squares.
- **NF01** (16 routines): NF01AD, NF01AY, NF01BA–NF01BY (BA, BB, BD, BE, BF, BP, BQ, BR, BS, BU, BV, BW, BX, BY) — nonlinear optimization (NF01).

Reference: [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md) (lines 372–413), [plans/slicot-100-percent-parity-roadmap.md](plans/slicot-100-percent-parity-roadmap.md) Phase 4 and 5.

## Conventions (from existing code)

- One file per routine: `src/<module>/<routine_lower>.rs` (e.g. `src/mc01/mc01md.rs`).
- Module crate: `src/<module>/mod.rs` with `pub mod <routine>;` for each routine.
- Register in [src/lib.rs](src/lib.rs): `pub mod mc01;` (and mc03, md03, nf01).
- Pure Rust only: no FFI; use [nalgebra](https://crates.io/crates/nalgebra) where needed (see e.g. [src/mb01/mb01ld.rs](src/mb01/mb01ld.rs)).
- Each routine file: doc comment with SLICOT name, public function matching SLICOT interface (arguments/return), and `#[cfg(test)] mod tests` with at least one `#[test]`.
- After implementation: set status to `done` in [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md); add/check off in [plans/remaining-slicot-functions-todo.md](plans/remaining-slicot-functions-todo.md); run `./scripts/validate_slicot_done.sh` and `./scripts/gen_features_table.sh`.

## Dependency and order

- **MC01** and **MC03** can be implemented in any order (no cross-dependency between them).
- **MD03** may call or share patterns with MB01/MB02; no hard dependency on MC01/MC03.
- **NF01** is standalone (nonlinear optimization); can be done in parallel or after others.

Recommended execution order: MC01 → MC03 → MD03 → NF01 (or MC01/MC03 in parallel, then MD03, then NF01).

## Section 1: MC01 (15 routines)

- Create `src/mc01/mod.rs` and 15 files `src/mc01/mc01md.rs` … `src/mc01/mc01xd.rs`.
- Implement each routine from SLICOT semantics (polynomial evaluation, composition, derivatives, etc.). If official Fortran or .dat references are available in the repo or from SLICOT release, use them; otherwise use SLICOT documentation or existing literature.
- Add `pub mod mc01;` to `src/lib.rs`.
- Update [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md): set all 15 MC01* rows to `done`.
- Add MC01 section to [plans/remaining-slicot-functions-todo.md](plans/remaining-slicot-functions-todo.md) with 15 checkboxes and check them when done.
- Run validation and feature-table scripts.

**Section 1 is complete only when all 15 MC01 routines are implemented, tested, documented in the mapping, and validated.**

## Section 2: MC03 (4 routines)

- Create `src/mc03/mod.rs` and 4 files: `mc03md.rs`, `mc03nd.rs`, `mc03nx.rs`, `mc03ny.rs`.
- Implement matrix pencil routines (e.g. equivalence, reduction) in pure Rust with nalgebra.
- Add `pub mod mc03;` to `src/lib.rs`.
- Update [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md): set MC03MD, MC03ND, MC03NX, MC03NY to `done`.
- Add MC03 section to [plans/remaining-slicot-functions-todo.md](plans/remaining-slicot-functions-todo.md) with 4 checkboxes.
- Run validation and feature-table scripts.

**Section 2 is complete only when all 4 MC03 routines are implemented, tested, documented, and validated.**

## Section 3: MD03 (7 routines)

- Create `src/md03/mod.rs` and 7 files: `md03ad.rs`, `md03ba.rs`, `md03bb.rs`, `md03bd.rs`, `md03bf.rs`, `md03bx.rs`, `md03by.rs`.
- Implement data fitting / nonlinear least-squares routines (driver and building blocks) in pure Rust.
- Add `pub mod md03;` to `src/lib.rs`.
- Update [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md): set all 7 MD03* rows to `done`.
- Add MD03 section to [plans/remaining-slicot-functions-todo.md](plans/remaining-slicot-functions-todo.md) with 7 checkboxes.
- Run validation and feature-table scripts.

**Section 3 is complete only when all 7 MD03 routines are implemented, tested, documented, and validated.**

## Section 4: NF01 (16 routines)

- Create `src/nf01/mod.rs` and 16 files: `nf01ad.rs`, `nf01ay.rs`, `nf01ba.rs` … `nf01by.rs` (BA, BB, BD, BE, BF, BP, BQ, BR, BS, BU, BV, BW, BX, BY).
- Implement nonlinear optimization routines (residual/jacobian interfaces and solvers) in pure Rust.
- Add `pub mod nf01;` to `src/lib.rs`.
- Update [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md): set all 16 NF01* rows to `done`.
- Add NF01 section to [plans/remaining-slicot-functions-todo.md](plans/remaining-slicot-functions-todo.md) with 16 checkboxes.
- Run validation and feature-table scripts.

**Section 4 is complete only when all 16 NF01 routines are implemented, tested, documented, and validated.**

## Plan completion

The **overall plan is complete** only when **all four sections** (MC01, MC03, MD03, NF01) are complete. Each section has its own TODO; do not mark the plan done until every section TODO is checked off and validation passes for all new routines.

## Notes

- No Fortran sources for MC01/MC03/MD03/NF01 were found under the repo; implementations will rely on SLICOT documentation or external references (e.g. SLICOT release docs or NAG/SLICOT user guides).
- [docs/FEATURES.md](docs/FEATURES.md) is regenerated via `./scripts/gen_features_table.sh`; run after mapping updates.
