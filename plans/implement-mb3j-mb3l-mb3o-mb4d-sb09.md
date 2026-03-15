# Implement mb3j, mb3l, mb3o, mb3p, mb4d, sb09 SLICOT routines

## Scope

| Module | Count | Routines |
|--------|-------|----------|
| mb3j | 1 | MB3JZP |
| mb3l | 1 | MB3LZP |
| mb3o | 1 | MB3OYZ |
| mb3p | 1 | MB3PYZ |
| mb4d | 3 | MB4DBZ, MB4DLZ, MB4DPZ |
| sb09 | 1 | SB09MD |

**Total:** 8 routines in 6 new modules.

## Reference behavior

- **SB09MD** — Full specification from SLICOT: compares two multivariable sequences M1(k), M2(k), outputs SS, SE, PRE. Real (double) only.
- **MB3JZP, MB3LZP, MB3OYZ, MB3PYZ** — Implemented as real dense fallbacks (solve A*X = B), mirroring MB03JZ/MB03LZ/MB03OY/MB03PY.
- **MB4DBZ, MB4DLZ, MB4DPZ** — Real dense fallbacks mirroring MB04DB/MB04DL/MB04DP.

## Implementation pattern

- New module layout: one directory per module under `src/`, with `mod.rs` and one `.rs` file per routine.
- MB3*/MB4D* reuse `crate::mb02::common::solve_ax_b`.
- SB09MD: full implementation from spec (SS, SE, PRE, TOL handling, INFO).
