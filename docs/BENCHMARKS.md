# SLICOT-rs benchmarks

Benchmarks run **every** SLICOT routine with **scaled problem sizes** so that runtimes grow with dimension and speed differences between routines (and between future implementations) are visible.

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

## Layout

- **`benches/common.rs`**: Shared size ladders (`SIZE_LADDER_N`, `SIZE_LADDER_POW2`) and helpers to build matrices/vectors (`matrix_nn`, `matrix_nm`, `state_space_matrices`, etc.).
- **`benches/all_routines.rs`**: Criterion groups per module; each routine is benchmarked at each size in the appropriate ladder. Stubs that take only `(n, m)` are registered with the same ladder so that when implemented, the harness already measures them.

The legacy **`benches/tb01md.rs`** benchmark is still available (`cargo bench --bench tb01md`) but uses smaller sizes (4–32); the full ladder for TB01MD is in `all_routines` (32–1024).
