# Fuzzing Rust vs Fortran (SLICOT)

This directory contains cargo-fuzz targets that compare Rust SLICOT implementations against the Fortran reference using SLICOT-style `.dat` input and `.res` output.

## Prerequisites

- **Rust** with LLVM (for libFuzzer): `rustup default nightly` or install the appropriate toolchain.
- **cargo-fuzz**: `cargo install cargo-fuzz`
- **Fortran executables** (optional): Build SLICOT and the file-I/O drivers (see [../docs/FORTRAN_BUILD.md](../docs/FORTRAN_BUILD.md) and [../tools/slicot-fortran/drivers/README.md](../tools/slicot-fortran/drivers/README.md)). Copy the built `TAB01ND` (and others) into a directory and set `SLICOT_EXAMPLES_DIR` to that directory. If unset, the fuzzer only runs the Rust implementation (no comparison).

## Targets

- **ab01nd_compare**: Parses fuzz input as AB01ND `.dat`; runs Rust `ab01nd`; if `SLICOT_EXAMPLES_DIR` is set, runs Fortran `TAB01ND` and compares outputs with relative tolerance. When Rust returns `INFO != 0` (e.g. stub), comparison is skipped (expect-fail).

## Corpus

Initial corpus is in `corpus/ab01nd_compare/` (e.g. `seed.dat`). You can add more `.dat`-style files. The fuzzer mutates these to generate new inputs.

## Run

From the project root:

```bash
cargo fuzz run ab01nd_compare
```

To use a corpus directory and run for a limited time:

```bash
cargo fuzz run ab01nd_compare -- -max_total_time=60
```

## Interpretation

- **No crash**: Either input was invalid (parse failed or dimensions out of range), or Rust ran and (if Fortran was run) outputs matched within relative tolerance, or Rust returned non-zero INFO (stub) and comparison was skipped.
- **Panic**: Rust and Fortran outputs differed (when both ran successfully) or another bug. Fix the implementation or adjust tolerance.
