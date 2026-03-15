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

To compare with Fortran implementations:

1. **Build SLICOT Fortran** as in [FORTRAN_BUILD.md](FORTRAN_BUILD.md) (gfortran, OpenBLAS, `./scripts/slicot-fortran/build_fortran.sh`).
2. **Time a routine** by either:
   - Running the corresponding example driver (e.g. `TMA02ED`) in a loop with the same problem size and averaging, or
   - Adding a small Fortran program that calls the SLICOT routine in a loop with matrices of size n×n (or n×m) and reports elapsed time.
3. **Record results** in the table below (same sizes as above) so Rust vs Fortran can be compared directly.

| Routine (Fortran) | n=32   | n=64   | n=128  | n=256  | n=512  | n=1024 |
|------------------|--------|--------|--------|--------|--------|--------|
| MA02ED           | —      | —      | —      | —      | —      | —      |
| MA02ES           | —      | —      | —      | —      | —      | —      |
| MB01MD           | —      | —      | —      | —      | —      | —      |
| TB01MD           | —      | —      | —      | —      | —      | —      |
| DLACPY           | —      | —      | —      | —      | —      | —      |
| DE01OD           | —      | —      | —      | —      | —      | —      |

*Fill in with Fortran timings (e.g. µs or ms per call) from your build to compare implementations.*

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

## Layout

- **`benches/common.rs`**: Shared size ladders (`SIZE_LADDER_N`, `SIZE_LADDER_POW2`) and helpers to build matrices/vectors (`matrix_nn`, `matrix_nm`, `state_space_matrices`, etc.).
- **`benches/all_routines.rs`**: Criterion groups per module; each routine is benchmarked at each size in the appropriate ladder. Stubs that take only `(n, m)` are registered with the same ladder so that when implemented, the harness already measures them.

The legacy **`benches/tb01md.rs`** benchmark is still available (`cargo bench --bench tb01md`) but uses smaller sizes (4–32); the full ladder for TB01MD is in `all_routines` (32–1024).
