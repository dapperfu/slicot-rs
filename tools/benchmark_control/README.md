# Control routines benchmark (with plotting)

This folder contains a **benchmark pipeline** for common SLICOT control routines: run timings with a **30-second cap per (routine, size)**, then plot results with **seaborn**.

## What it does

- **Rust binary** (`cargo run --release --bin bench_csv`): Runs MA02ED, MA02ES, TB01MD, DLACPY_SLC, MB01MD, and DE01OD at sizes 8, 16, 32, 64, 128, 256, 512 (power-of-two only for DE01OD). For each (routine, size) it iterates for up to **30 seconds** and prints one CSV row: `routine,n,time_us`.
- **Python script** (`run_benchmarks.py`): Invokes the binary, parses the CSV, and produces a seaborn line plot (n vs time per call, one line per routine) saved as `benchmark_control.png`.

## Run from repo root

```bash
python3 tools/benchmark_control/run_benchmarks.py
```

**Requirements:** `pandas`, `seaborn`, `matplotlib` (e.g. `pip install pandas seaborn matplotlib`). The script will print an error if they are missing.

**Note:** The full run can take several minutes because each (routine, size) is allowed up to 30s of compute.
