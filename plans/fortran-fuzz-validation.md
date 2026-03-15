# Fortran compile and fuzz validation plan

## Goal

- Compile the Fortran SLICOT reference (SLICOT-Reference/) and use **plain text files** (SLICOT-style `.dat` input / `.res` output) as the only interface between Fortran and the fuzzer.
- Use a **fuzzer** (LLVM-based cargo-fuzz preferred) with a **corpus** of existing `.dat` files plus mutations to validate that **Rust implementations match Fortran output 1:1** (Fortran = reference). Comparison uses **relative tolerance** for floats. **Stubs**: include all implemented routines; when Rust returns `INFO != 0` (e.g. stub), treat as expect-fail (no output comparison; optional logging).

No CI in scope; assume **system OpenBLAS** for Fortran build.

---

## 1. Save plan and document decisions

- Save this plan under plans/ with a descriptive name (e.g. `fortran-fuzz-validation.md`).
- Document: reference = Fortran; format = SLICOT .dat/.res; routines = all implemented (stubs expect-fail); fuzzer = corpus + mutations; tolerance = relative; BLAS = system OpenBLAS; no CI.

---

## 2. Compile Fortran SLICOT and dependencies

**Build order** (from SLICOT-Reference/make_Unix.inc and makefiles):

1. **lpkaux**  
   - Build `lpkaux.a` from SLICOT-Reference/src/lapack_aux/ (e.g. `dlacpy_slc.f`, `dlatzm.f`, `dgegs.f`, `dgegv.f`).  
   - May require a small Makefile or script in SLICOT-Reference/ that compiles these and produces `lpkaux.a` in a known location (e.g. `../` relative to `src/` as in current `LPKAUXLIB`).

2. **SLICOT library**  
   - Use SLICOT-Reference/src/makefile_Unix to build `slicot.a` (include path to `lpkaux.a` and `-lopenblas`).  
   - Document in README or docs/: need gfortran, OpenBLAS; no automatic install of OpenBLAS.

3. **Example drivers (T*)**  
   - Keep building drivers via SLICOT-Reference/examples/makefile_Unix so each T*.f is linked with slicot.a, lpkaux, LAPACK/BLAS.  
   - Executables can stay in SLICOT-Reference/examples/ (e.g. TAB01ND, TMB01TD, …).

**Deliverable**: Script or Makefile at repo root or under SLICOT-Reference/ that: builds lpkaux → slicot → example drivers; documented prerequisite: system OpenBLAS.

---

## 3. Plain text I/O: SLICOT-style .dat / .res

- **Input**: Reuse existing SLICOT example format (e.g. SLICOT-Reference/examples/data/AB01ND.dat): optional title line, then numeric lines (N, M, TOL, JOBZ, then matrix rows). No new format; same layout as current `.dat` files.
- **Output**: Reuse the same style as current `.res` (e.g. SLICOT-Reference/examples/results/AB01ND.res): text lines with scalars and matrices so that a parser can extract INFO, NCONT, INDCON, and numeric blocks (A, B, Z, etc.) in a deterministic order.
- **Adapt existing drivers** to **file-based I/O**:
  - Instead of fixed NIN=5 / NOUT=6, take two arguments: **input file path**, **output file path** (or read from stdin and write to stdout when no args, to keep backward compatibility with make).
  - Example: `TAB01ND [input.dat [output.res]]` — if no args, use stdin/stdout; if one or two args, use files. This allows the fuzzer to pass corpus/ab01nd/input.dat and fortran_out.res.
- **Rust side**: Add a small I/O layer (in a crate or tests/ helper) that:
  - **Writes** a `.dat`-style file from Rust (dimensions, scalars, matrices) for a given routine and input.
  - **Reads** a `.res`-style file and parses INFO, scalars, and matrices into a canonical structure (e.g. struct per routine family).
  - **Writes** Rust output in the same `.res`-style so one can diff or compare field-by-field with the Fortran-generated `.res` using relative tolerance.

**Deliverable**: (1) Driver convention and one example driver (e.g. TAB01ND) adapted to optional file args. (2) Spec or comment describing the exact .dat/.res layout for that routine. (3) Rust helper to read/write .dat and .res for at least one routine (pilot).

---

## 4. Fuzzer: cargo-fuzz + corpus + comparison

- **Tool**: Use **cargo-fuzz** (libFuzzer) for in-process, fast fuzzing and good Rust integration. Add `[package.metadata.cargo-fuzz]` and a fuzz target.
- **Corpus**:  
  - **Initial corpus**: Copy or symlink existing SLICOT-Reference/examples/data/*.dat into `fuzz/corpus/<routine_name>/` (e.g. `fuzz/corpus/ma01ad/`, `fuzz/corpus/ab01nd/`).  
  - **Mutation**: Fuzzer receives bytes; treat them as **content of a .dat file** (or a length-prefixed blob). If parsing fails or dimensions are out of range (e.g. N > 20), return early. Otherwise run the pipeline. This gives structure-aware mutations from a .dat-style corpus.
- **Pipeline per run** (for a single routine, e.g. MA01AD or AB01ND):
  1. **Parse** fuzz input as .dat (or skip if invalid).
  2. **Run Fortran**: write parsed input to a temp .dat, invoke T* with input path and output path, read Fortran .res → reference output.
  3. **Run Rust**: call the corresponding Rust routine with the same inputs; build Rust output in the same .res-like structure.
  4. **Compare**:  
     - If Rust returns INFO != 0 (e.g. stub): **expect-fail** — do not assert equality; optionally log "Rust stub, Fortran ref = ...".  
     - If Rust returns INFO == 0: compare all relevant outputs (INFO, scalars, matrices) with **relative tolerance** (e.g. 1e-10 * max(|a|, |b|, 1)).
- **Routines**: Include **all** routines that have both a Fortran driver and a Rust implementation (including stubs). One fuzz target per routine (or a single target that dispatches by routine id from a byte in the input) to keep corpus and coverage per routine.
- **Environment**: Fuzzer must know the path to the Fortran executables (e.g. SLICOT_EXAMPLES_DIR or CARGO_MANIFEST_DIR/../SLICOT-Reference/examples). If the binary is missing, skip the Fortran run and only run Rust (or skip the comparison).

**Deliverable**: fuzz/ directory with at least one fuzz target (e.g. fuzz_targets/ma01ad.rs or fuzz_targets/slicot_compare.rs), corpus dirs seeded with existing .dat files, and a README in fuzz/ explaining how to build Fortran, run the fuzzer, and interpret results.

---

## 5. Relative tolerance and output comparison

- Implement a small **comparison** module (e.g. in tests/ or fuzz/):  
  - For each scalar: `|a - b| <= rel_tol * max(|a|, |b|, 1.0)` (or use existing crate like approx with relative tolerance).  
  - For matrices: same per element; optionally report first index where tolerance is exceeded.
- Use this in the fuzz target when Rust INFO == 0 and Fortran ran successfully.

---

## 6. Documentation and scripts

- **README or docs**:  
  - How to build Fortran (lpkaux → slicot → examples); OpenBLAS requirement.  
  - How to run the fuzzer: `cargo fuzz run <target>`, corpus location, and that Fortran binaries must be present for comparison.  
  - Meaning of expect-fail for stubs and how to add new routines/corpus.
- **Optional**: A small script that builds Fortran, then runs a short fuzz run for one routine to sanity-check the pipeline.

---

## Architecture (high level)

- Corpus: .dat files → Parse .dat → Run Fortran T* and Run Rust impl → Compare with relative tolerance.

---

## Files to add or touch (summary)

| Area | Files / locations |
|------|-------------------|
| Plan | plans/fortran-fuzz-validation.md |
| Fortran build | SLICOT-Reference/ — ensure lpkaux build, then make in src and examples; optional top-level script |
| Driver I/O | SLICOT-Reference/examples/TAB01ND.f (and later others) — optional file args |
| Rust I/O | New helper in tests/ or fuzz/ to read/write .dat and .res |
| Fuzzer | fuzz/Cargo.toml, fuzz/fuzz_targets/*.rs, fuzz/corpus/<routine>/ |
| Docs | README section or docs/ for Fortran build and fuzz usage |

---

## Open points (to resolve during implementation)

- Exact .res parsing rules for every routine family (pilot with one routine, then generalize).
- Whether to use one fuzz target per routine or a single dispatcher target (trade-off: many targets vs. one target with routine selector).
- Upper bounds for dimensions (N, M) when parsing fuzz bytes to avoid OOM (e.g. cap at 20 to match SLICOT NMAX/MMAX).

---

## Decisions (recorded)

| Decision | Choice |
|----------|--------|
| Reference implementation | Fortran is reference; validate Rust matches Fortran. |
| Text format | SLICOT .dat / .res style (free-form numbers, one line per row). |
| Routines in scope | All implemented (including stubs). |
| Fortran build | Adapt existing T* example drivers to accept input/output file paths. |
| Fuzzer | cargo-fuzz (libFuzzer); any backend that fits. |
| Float comparison | Relative tolerance (e.g. 1e-10 * max(\|a\|, \|b\|, 1)). |
| CI | No CI in scope; local/optional only. |
| BLAS/LAPACK | System OpenBLAS; document only, no automatic install. |
| Stub handling | Expect-fail: when Rust returns INFO != 0, no output comparison; optional logging. |
| Input generation | Shared corpus of .dat files + fuzz mutations. |
