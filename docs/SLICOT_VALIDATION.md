# SLICOT validation report: Pure Rust and tested

This document reports the validation that every routine marked **done** in [SLICOT_MAPPING.md](SLICOT_MAPPING.md) is implemented in pure Rust (no FFI) and has at least one automated test.

## Validation date

Validation performed as part of the plan "Validate All SLICOT Functions: Pure Rust and Tested".

## Scope

- **Source of truth**: [docs/SLICOT_MAPPING.md](SLICOT_MAPPING.md)
- **Routines validated**: All rows with status `done` (150 routines)
- **Criteria**: (1) No FFI/extern/Fortran linkage in implementation; (2) at least one `#[cfg(test)]` module with at least one `#[test]` in the routine's `.rs` file (or equivalent coverage via integration/fuzz).

## Results summary

| Check            | Result |
|------------------|--------|
| Pure Rust (no FFI) | **PASS** — No `extern`, `ffi`, `libslicot`, or Fortran linkage found in `src/**/*.rs`. |
| All done have tests | **PASS** — After fixes: every done routine's implementation file contains at least one test. |
| Total done routines | 150 |
| Gaps found         | 1 (DG01NY lacked tests; fixed). |

## Fixes applied

1. **DG01NY** (`src/dg01/dg01ny.rs`): Had no `#[cfg(test)]` or `#[test]`. Added a `#[cfg(test)] mod tests` with three tests: `test_dg01ny_direct_n2`, `test_dg01ny_inverse_n2`, `test_dg01ny_n_too_small`.

## Per-routine summary

All 150 done routines were checked:

- **Pure Rust**: Each implementation file under `src/<module>/<function>.rs` was searched for `extern`, `ffi`, `libslicot`, `.so`, and similar; none were found.
- **Has test**: Each implementation file was checked for the presence of `#[cfg(test)]`; one file (dg01ny.rs) initially lacked it and was updated as above.

Routine list (SLICOT name | module | rust function): see [SLICOT_MAPPING.md](SLICOT_MAPPING.md) — all rows with status `done`.

## Automation

The script [../scripts/validate_slicot_done.sh](../scripts/validate_slicot_done.sh) can be run to re-validate that every done routine remains pure Rust and has tests. Run from project root: `./scripts/validate_slicot_done.sh`.

## Completion

Validation is complete when: (1) every done routine passes both checks, (2) all identified gaps are fixed, and (3) this report exists. All conditions are satisfied.
