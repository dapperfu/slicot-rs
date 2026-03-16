# Full implementation requirement

## Requirement

**No minimal or stub implementations.** Every SLICOT routine in the mapping SHALL be fully implemented according to SLICOT specification and behavior. The task is not complete until all 625 routines are fully implemented.

## Current status

- **Fully implemented (spec-compliant):** MA01CD, MB01YD, MB01ZD, MA02MD (and the 150 routines that were done before the Phase 1/2 stubs).
- **Stubbed (must be replaced):** The 55 routines added in Phase 1 and Phase 2 (AB13AD/AX/BD/CD/DX/ID, MA02HD/HZ/ID/IZ/JD/JZ/MD/MZ/OD/OZ/SD, MB02CD–MB02YD) are currently minimal stubs. Each must be replaced with a full implementation per SLICOT documentation.
- **Not started (must be fully implemented):** The remaining 420 routines (Phases 3–10). No stubs; each new routine must be implemented fully from the start.

## Approach

1. **Source of truth:** SLICOT routine documentation at https://www.slicot.org/objects/software/shared/doc/ and, where available, Fortran source in a SLICOT-Reference tree.
2. **Per routine:** (1) Read SLICOT doc (and Fortran if present); (2) Implement algorithm in pure Rust (nalgebra, no FFI); (3) Add tests that verify numerical behavior, not just return codes; (4) Set status to `done` in SLICOT_MAPPING.md only when the implementation is complete and correct.
3. **Dependencies:** Some routines (e.g. AB13AD) depend on others (Schur form, reordering, Gramians). Implement in dependency order or implement required building blocks first.

## Completion criteria

100% feature parity is achieved only when:

- All 625 routines have status `done` in docs/SLICOT_MAPPING.md.
- Every routine has a full, spec-compliant implementation (no stubs).
- Every routine has at least one test that checks meaningful behavior.
- ./tools/validate_slicot_done.sh passes.

Until then, the task is not complete.
