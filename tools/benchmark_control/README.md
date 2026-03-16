# Control routines benchmark (with plotting)

This folder contains a **benchmark pipeline** for common SLICOT control routines: run **Rust** timings (30s cap per size) and the **Fortran** SLICOT benchmark, then plot **Rust vs Fortran** with seaborn.

## What it does

- **Rust** (`cargo run --release --bin bench_csv`): Runs MA02ED, MA02ES, TB01MD, DLACPY_SLC, MB01MD, DE01OD at sizes 8..512 (up to 30s per (routine, size)), outputs CSV `routine,n,time_us`.
- **Fortran** (`tools/slicot-fortran/run_fortran_benchmarks.sh`): Builds and runs the Fortran SLICOT benchmark for MA02ED, MA02ES, DLACPY_SLC at n=32,64,128,256,512,1024, outputting the same CSV format.
- **Python script** (`run_benchmarks.py`): Runs both, merges data, and produces a seaborn line plot (n vs time per call) with **Rust vs Fortran** per routine, saved as `benchmark_control.png`. If the Fortran build or run fails, only Rust data is plotted.

## Run from repo root

```bash
python3 tools/benchmark_control/run_benchmarks.py
```

**Requirements:** `pandas`, `seaborn`, `matplotlib`. For Fortran comparison: Fortran SLICOT must build (see [FORTRAN_BUILD.md](../../docs/FORTRAN_BUILD.md)); the script uses `tools/slicot-fortran/run_fortran_benchmarks.sh`.

**Note:** The full run can take several minutes (Rust: 30s per size; Fortran: build once then run).
