# SLICOT-rs benchmarks

Benchmarks run **every** SLICOT routine with **scaled problem sizes** so that runtimes grow with dimension and speed differences between routines (and between future implementations) are visible.

## Benchmark results (Rust vs Fortran)

Representative **Rust (slicot-rs)** timings for key routines at scaled sizes. Use these to compare with Fortran SLICOT timings from the same problem sizes.

### Rust (slicot-rs) — time per iteration

| Routine    | n=32    | n=64      | n=128     | n=256      | n=512       | n=1024       |
|-----------|---------|-----------|-----------|------------|-------------|--------------|
| MA02ED    | ~291 ns | ~9.06 µs  | ~50.7 µs  | ~103 µs    | ~1.02 ms    | ~4.76 ms     |
| MA02ES    | ~2.49 µs| ~9.78 µs  | ~45.3 µs  | ~76.0 µs   | ~845 µs     | (run bench) |
| MB01MD    | (run)   | (run)     | (run)     | (run)      | (run)       | (run)        |
| TB01MD    | (run)   | (run)     | (run)     | (run)      | (run)       | (run)        |
| DLACPY_SLC| (run)   | (run)     | (run)     | (run)      | (run)       | (run)        |
| DE01OD    | (run)   | (run)     | (run)     | (run)      | (run)       | (run)        |

Throughput (elements/s) for matrix routines is reported by Criterion (e.g. MA02ED ~220–640 Melem/s depending on n). Run the full suite to fill in all cells:

```bash
cargo bench --bench all_routines -- --noplot
```

### Fortran (SLICOT reference) — for comparison

Build and run the Fortran benchmark driver (same size ladder as Rust):

```bash
./tools/slicot-fortran/run_fortran_benchmarks.sh
```

This script builds `lpkaux.a` and `slicot.a` if needed (see [FORTRAN_BUILD.md](FORTRAN_BUILD.md); the full `examples` target may fail without OpenBLAS, but the benchmark driver only needs the two libraries), then compiles and runs `tools/slicot-fortran/bench/bench_slicot.f90`. If OpenBLAS is not installed, the script retries with `-lblas -llapack`. Results are printed as µs/call per (routine, n).

**Example Fortran results** (gfortran -O2, ref BLAS/LAPACK, one machine):

| Routine (Fortran) | n=32   | n=64   | n=128  | n=256   | n=512   | n=1024   |
|------------------|--------|--------|--------|---------|---------|----------|
| MA02ED           | 1.41 µs| 2.09 µs| 14.9 µs| 90.1 µs | 780 µs  | 3030 µs  |
| MA02ES           | 1.51 µs| 2.31 µs| 17.3 µs| 131 µs  | 1667 µs | 6823 µs  |
| MB01MD           | —      | —      | —      | —       | —       | —        |
| TB01MD           | —      | —      | —      | —       | —       | —        |
| DLACPY           | 1.30 µs| 1.79 µs| 3.62 µs| 14.5 µs | 581 µs  | 2304 µs  |
| DE01OD           | —      | —      | —      | —       | —       | —        |

**Rust vs Fortran (same n):** At n=256, Rust MA02ED ~103 µs vs Fortran ~90 µs; at n=1024, Rust ~4.76 ms vs Fortran ~3.03 ms. Run both suites on your machine to compare; Rust uses pure nalgebra (no BLAS), Fortran uses BLAS (DCOPY) inside MA02ED.

---

## How to run

```bash
cargo bench --bench all_routines
```

Results are written under `target/criterion/`. To run a specific group (e.g. only `ma02` routines):

```bash
cargo bench --bench all_routines -- ma02
```

To run a single routine and size (e.g. `ma02ed` at `n=256`):

```bash
cargo bench --bench all_routines -- ma02/ma02ed/n256
```

## Size policy (no tiny sizes)

- **Primary ladder** for state size / matrix dimension: `n ∈ {32, 64, 128, 256, 512, 1024}`. We do **not** run huge numbers of tiny problems (e.g. 1e5 iterations of 2×2 matrices); each benchmark is a (routine, size) pair with meaningful dimensions.
- For routines that use **(n, m)** or **(n, m, p)**, secondary dimensions are derived from `n` (e.g. `m = n/2`, `p = n/2`).
- **FFT/signal** routines (e.g. DE01OD, DG01MD) use a **power-of-two** ladder: `n ∈ {64, 128, 256, 512, 1024, 2048}`.

So you can expect:

- **Real implementations** (e.g. MA02ED, MB01MD, TB01MD, DLACPY_SLC, DE01OD): time per iteration **increases** with `n`; throughput (elements/s) may be reported.
- **Stubs** (many AB09*, AB13*, IB01*, etc.): flat, low time until they are implemented; the same benchmark names will then show the new cost.

## Interpreting results

- Criterion prints **time per iteration** and (where set) **throughput** (e.g. elements/s) for each (routine, size).
- Compare routines at a fixed size (e.g. `n=256`) to see which are slower.
- Compare the same routine across sizes to confirm scaling (e.g. MA02ED should scale roughly with `n²`).
- To **compare Rust vs Fortran**: run the Rust benchmarks above, then time the same routine and problem sizes in the Fortran SLICOT build, and fill in the [Fortran results table](#fortran-slicot-reference--for-comparison) so both implementations can be compared side by side.

## Benchmark script with plotting (30s cap per size)

A separate pipeline runs **common control routines** (MA02ED, MA02ES, TB01MD, DLACPY_SLC, MB01MD, DE01OD) at sizes from 8 up to 512, with **at most 30 seconds** of compute per (routine, size), and plots the results with **seaborn**. Run from the repo root:

```bash
python3 tools/benchmark_control/run_benchmarks.py
```

Requirements: `pandas`, `seaborn`, `matplotlib`. The script runs `cargo run --release --bin bench_csv` to produce CSV, then generates `tools/benchmark_control/benchmark_control.png`. See `tools/benchmark_control/README.md` for details.

---

## Layout

- **`benches/common.rs`**: Shared size ladders (`SIZE_LADDER_N`, `SIZE_LADDER_POW2`) and helpers to build matrices/vectors (`matrix_nn`, `matrix_nm`, `state_space_matrices`, etc.).
- **`benches/all_routines.rs`**: Criterion groups per module; each routine is benchmarked at each size in the appropriate ladder. Stubs that take only `(n, m)` are registered with the same ladder so that when implemented, the harness already measures them.
- **Fortran**: `tools/slicot-fortran/bench/bench_slicot.f90` times MA02ED, MA02ES, DLACPY_SLC at the same n ladder; `tools/slicot-fortran/run_fortran_benchmarks.sh` builds (lpkaux.a, slicot.a) and runs it. Uses `-lblas -llapack` if OpenBLAS is not installed.
- **Benchmark script**: `src/bin/bench_csv.rs` prints CSV (routine, n, time_us) with a 30s cap per size; `tools/benchmark_control/run_benchmarks.py` runs it and plots with seaborn.

The legacy **`benches/tb01md.rs`** benchmark is still available (`cargo bench --bench tb01md`) but uses smaller sizes (4–32); the full ladder for TB01MD is in `all_routines` (32–1024).
