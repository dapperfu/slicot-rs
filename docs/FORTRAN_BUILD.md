# Building SLICOT Fortran reference

The Fortran SLICOT library and example drivers are used to validate Rust implementations (see [plans/fortran-fuzz-validation.md](../plans/fortran-fuzz-validation.md)).

## Prerequisites

- **gfortran** (GNU Fortran)
- **OpenBLAS** (system library; link with `-lopenblas`). No automatic install is provided; install via your system package manager (e.g. `libopenblas-dev` on Debian/Ubuntu).

## Build order

1. **lpkaux.a** — built from `SLICOT-Reference/src/lapack_aux/*.f` into `SLICOT-Reference/lpkaux.a`
2. **slicot.a** — built from `SLICOT-Reference/src/*.f` into `SLICOT-Reference/slicot.a` (depends on lpkaux.a)
3. **Example drivers (T*)** — built in `SLICOT-Reference/examples/` (e.g. `TAB01ND`, `TMB01TD`), linked with slicot.a, lpkaux.a, and OpenBLAS

## How to build

From the project root (with a local `SLICOT-Reference/` tree present):

```bash
./scripts/slicot-fortran/build_fortran.sh
```

Or run make explicitly (the directory contains a Windows `makefile`, so use `-f` to select the Unix Makefile):

```bash
make -C SLICOT-Reference -f scripts/slicot-fortran/Makefile
```

This will:

- Create `SLICOT-Reference/lpkaux.a`
- Create `SLICOT-Reference/slicot.a`
- Build all example executables in `SLICOT-Reference/examples/` and run their default data (producing `.exa` result files)

To clean:

```bash
./scripts/slicot-fortran/build_fortran.sh clean
```
or
```bash
make -C SLICOT-Reference -f scripts/slicot-fortran/Makefile clean
```

## File-based I/O for fuzzer

For the fuzz pipeline, example drivers must accept optional input/output file paths. A copy of `TAB01ND` with this behaviour is in the repo:

- **Convention:** `T* [input.dat [output.res]]` — no args = stdin/stdout; one arg = input file; two args = input and output files.
- **Copy** the driver into your tree and rebuild:  
  `cp scripts/slicot-fortran/drivers/TAB01ND.f SLICOT-Reference/examples/TAB01ND.f`  
  then run the Fortran build again. See `scripts/slicot-fortran/drivers/README.md` for .dat/.res layout.

## Benchmark driver (no OpenBLAS required for timing)

To time Fortran SLICOT routines (MA02ED, MA02ES, DLACPY_SLC) at the same sizes as the Rust benchmarks, you only need `lpkaux.a` and `slicot.a`. The example drivers (step 3) require OpenBLAS; the benchmark driver can link with reference BLAS/LAPACK instead.

From the project root:

```bash
./scripts/slicot-fortran/run_fortran_benchmarks.sh
```

This builds `lpkaux.a` and `slicot.a` if missing (run `./scripts/slicot-fortran/build_fortran.sh lpkaux.a slicot` if the full build fails at the examples step), then compiles and runs `scripts/slicot-fortran/bench/bench_slicot.f90`. If OpenBLAS is not installed, the script uses `-lblas -llapack`. See [BENCHMARKS.md](BENCHMARKS.md) for the results layout and Rust vs Fortran comparison.

## Configuration

Build options (compiler, flags, BLAS) are set in `SLICOT-Reference/make_Unix.inc`. By default: `FORTRAN = gfortran`, `OPTS = -O2 -fPIC -g`, `BLASLIB = -lopenblas`, `LAPACKLIB = -lopenblas`. Adjust if your OpenBLAS or LAPACK is installed elsewhere.
