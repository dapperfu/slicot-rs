# Fortran 1:1 validation

This folder contains validation results comparing each Rust SLICOT routine against the FORTRAN reference. Each routine that has SLICOT example data (`.dat`/`.res`) is run with the same input on both Fortran and Rust; outputs are compared with relative tolerance 1e-10.

## How to run

- **Canonical full run** (build FORTRAN, build Rust, run 1:1 validation; fails on any mismatch):  
  `./run.sh`  
  From project root. This builds Fortran (unless `--no-fortran`), builds Rust tests, runs the validation test, and exits non-zero if any validated routine does not match FORTRAN. See [../README.md](../README.md#one-command-to-build-and-validate).

- **Validation only** (assumes Fortran and Rust already built):  
  `./tools/validation/run_all.sh`  
  Builds Fortran if needed, then runs the validation runner and overwrites `validation/*.md`.

- **Subset (cargo test)** (skips if Fortran missing):  
  `SLICOT_EXAMPLES_DIR=/path/to/SLICOT-Reference/examples cargo test --test fortran_validation`  
  If `SLICOT_EXAMPLES_DIR` is unset or the driver is missing, the test skips and exits 0.

## Summary

| Module | Pass | Fail | No ref | Doc |
|--------|------|------|--------|-----|
| ab01 | 0 | 1 | 0 | [ab01](ab01.md) |
| ab04 | 0 | 0 | 0 | [ab04](ab04.md) |
| ab05 | 0 | 0 | 0 | [ab05](ab05.md) |
| ab07 | 0 | 0 | 0 | [ab07](ab07.md) |
| ab08 | 0 | 0 | 0 | [ab08](ab08.md) |
| ab09 | 0 | 0 | 0 | [ab09](ab09.md) |
| ab13 | 0 | 0 | 0 | [ab13](ab13.md) |
| ag08 | 0 | 0 | 0 | [ag08](ag08.md) |
| bb01 | 0 | 0 | 0 | [bb01](bb01.md) |
| bb02 | 0 | 0 | 0 | [bb02](bb02.md) |
| bb03 | 0 | 0 | 0 | [bb03](bb03.md) |
| bb04 | 0 | 0 | 0 | [bb04](bb04.md) |
| bd01 | 0 | 0 | 0 | [bd01](bd01.md) |
| bd02 | 0 | 0 | 0 | [bd02](bd02.md) |
| de01 | 0 | 0 | 0 | [de01](de01.md) |
| df01 | 0 | 0 | 0 | [df01](df01.md) |
| dg01 | 0 | 0 | 0 | [dg01](dg01.md) |
| dk01 | 0 | 0 | 0 | [dk01](dk01.md) |
| fb01 | 0 | 0 | 0 | [fb01](fb01.md) |
| fd01 | 0 | 0 | 0 | [fd01](fd01.md) |
| ib01 | 0 | 0 | 0 | [ib01](ib01.md) |
| ib03 | 0 | 0 | 0 | [ib03](ib03.md) |
| mb01 | 0 | 0 | 0 | [mb01](mb01.md) |
| mb02 | 0 | 0 | 0 | [mb02](mb02.md) |
| mb03 | 0 | 0 | 0 | [mb03](mb03.md) |
| mb04 | 0 | 0 | 0 | [mb04](mb04.md) |
| mb05 | 0 | 0 | 0 | [mb05](mb05.md) |
| mb4d | 0 | 0 | 0 | [mb4d](mb4d.md) |
| mc01 | 0 | 0 | 0 | [mc01](mc01.md) |
| mc03 | 0 | 0 | 0 | [mc03](mc03.md) |
| md03 | 0 | 0 | 0 | [md03](md03.md) |
| sb01 | 0 | 0 | 0 | [sb01](sb01.md) |
| sb02 | 0 | 0 | 0 | [sb02](sb02.md) |
| sb03 | 0 | 0 | 0 | [sb03](sb03.md) |
| sb04 | 0 | 0 | 0 | [sb04](sb04.md) |
| sb06 | 0 | 0 | 0 | [sb06](sb06.md) |
| sb08 | 0 | 0 | 0 | [sb08](sb08.md) |
| sb09 | 0 | 0 | 0 | [sb09](sb09.md) |
| sb10 | 0 | 0 | 0 | [sb10](sb10.md) |
| sb16 | 0 | 0 | 0 | [sb16](sb16.md) |
| sg02 | 0 | 0 | 0 | [sg02](sg02.md) |
| sg03 | 0 | 0 | 0 | [sg03](sg03.md) |
| tb01 | 0 | 0 | 0 | [tb01](tb01.md) |
| tb03 | 0 | 0 | 0 | [tb03](tb03.md) |
| tb04 | 0 | 0 | 0 | [tb04](tb04.md) |
| tb05 | 0 | 0 | 0 | [tb05](tb05.md) |
| tc01 | 0 | 0 | 0 | [tc01](tc01.md) |
| tc04 | 0 | 0 | 0 | [tc04](tc04.md) |
| tc05 | 0 | 0 | 0 | [tc05](tc05.md) |
| td03 | 0 | 0 | 0 | [td03](td03.md) |
| td04 | 0 | 0 | 0 | [td04](td04.md) |
| td05 | 0 | 0 | 0 | [td05](td05.md) |
| tf01 | 0 | 0 | 0 | [tf01](tf01.md) |
| tg01 | 0 | 0 | 0 | [tg01](tg01.md) |
| ud01 | 0 | 0 | 0 | [ud01](ud01.md) |

## Module docs

One markdown file per Rust module that has at least one routine with reference data:

- [ab01](ab01.md)
- *(links added by run_all.sh)*

## Tolerance

Float comparison uses a single relative tolerance: `max(|a|, |b|, 1.0) * 1e-10`.

## See also

- [FORTRAN_BUILD.md](../docs/FORTRAN_BUILD.md) — build Fortran SLICOT and examples
- [SLICOT_MAPPING.md](../docs/SLICOT_MAPPING.md) — SLICOT → Rust mapping
