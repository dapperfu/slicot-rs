# Fuzzing Rust vs Fortran

Rust SLICOT implementations are validated against the Fortran reference by fuzzing with SLICOT-style `.dat` inputs and comparing outputs (see [plans/fortran-fuzz-validation.md](../plans/fortran-fuzz-validation.md)).

## Quick start

1. **Build Fortran** (optional, for comparison):  
   [FORTRAN_BUILD.md](FORTRAN_BUILD.md) and copy the file-I/O driver per [scripts/slicot-fortran/drivers/README.md](../scripts/slicot-fortran/drivers/README.md). Set `SLICOT_EXAMPLES_DIR` to the directory containing `TAB01ND`, etc.

2. **Install cargo-fuzz** (if needed):  
   `cargo install cargo-fuzz`

3. **Run the fuzzer** (from project root):  
   `cargo fuzz run ab01nd_compare`

See [fuzz/README.md](../fuzz/README.md) for target list, corpus location, and interpretation.

## Expect-fail for stubs

When a Rust routine is still a stub (returns `INFO != 0`), the fuzzer does not compare outputs; the run is treated as expect-fail. Once the implementation is complete, comparison runs and any mismatch causes a panic.
