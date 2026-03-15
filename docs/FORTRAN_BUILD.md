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

## Configuration

Build options (compiler, flags, BLAS) are set in `SLICOT-Reference/make_Unix.inc`. By default: `FORTRAN = gfortran`, `OPTS = -O2 -fPIC -g`, `BLASLIB = -lopenblas`, `LAPACKLIB = -lopenblas`. Adjust if your OpenBLAS or LAPACK is installed elsewhere.
