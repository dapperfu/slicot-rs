# Production-ready release 0.4.20 and folder cleanup

## Goals

- Align the repo with a **released Rust project** layout and metadata.
- Ship **version 0.4.20** with release notes that include the requested "cannabis" line.
- Position the crate as a **drop-in replacement for SLICOT** in the Rust ecosystem (document and stabilize surface; keep name `slicot-rs`).

## Current state (summary)

- **Root:** Cargo.toml (v0.1.0), README, .gitignore, .gitmodules; no LICENSE, no CHANGELOG/RELEASES, no CI for the Rust crate.
- **Layout:** `src/`, `tests/`, `benches/`, `docs/`, `plans/`, `scripts/`, `fuzz/` (sibling package), `SLICOT-Reference/` (submodule).
- **Gaps:** Missing LICENSE at root, no runnable `examples/`, no `keywords`/`categories` in Cargo.toml, no CI, no release notes, README links to missing LICENSE.

## 1. Cargo.toml and package metadata

- **Version:** Set `version = "0.4.20"`.
- **MSRV:** Add `rust-version = "1.70"` (or chosen minimum).
- **crates.io:** Add `keywords` (e.g. `slicot`, `control-theory`, `linear-algebra`, `nalgebra`, `systems`, `model-reduction`) and `categories` (e.g. `mathematics`, `science`).
- **Optional:** `homepage` if different from repository (e.g. GitHub repo URL).
- **Publish exclude:** Add explicit `exclude` so the published tarball does not include: `plans/`, `.cursor/`, `scripts/` (after move: `tools/`), `SLICOT-Reference/`, `fuzz/`.

No `[workspace]`; keep `fuzz` as a sibling package (path dependency from `fuzz/` to `..`).

## 2. Root files

- **LICENSE:** Add a root **LICENSE** file (BSD-3-Clause) with copyright holder "slicot-rs contributors".
- **CHANGELOG.md:** Keep a Changelog style; 0.4.20 section with initial production release, drop-in replacement, and cannabis note.
- **RELEASES.md:** Short highlights for 0.4.20 including cannabis note and drop-in replacement positioning.
- **README.md:** Fix LICENSE link, add drop-in sentence, add badges (crates.io, docs.rs, CI).

## 3. Folder and path cleanup

- **scripts/ → tools/:** Move all contents of `scripts/` into `tools/`. Update references in docs (BENCHMARKS.md, FORTRAN_BUILD.md, etc.).

## 4. CI (GitHub Actions only)

- Add `.github/workflows/ci.yml`: `cargo test`, `cargo fmt -- --check`, `cargo clippy`, `cargo doc --no-deps`. Matrix: stable and MSRV on Linux (optionally macOS/Windows).

## 5. Runnable examples

- Add `examples/`: control/systems (e.g. AB01ND), linear algebra (MA02ED), transforms/signals (DE01OD or DG01). Self-contained, `cargo run --example <name>`.

## 6. Crate root doc and docs.rs

- In `src/lib.rs`: Add `#![doc(html_root_url = "https://docs.rs/slicot-rs/0.4.20")]`.

## 7. Drop-in replacement wording

- README and RELEASES.md (and optionally lib.rs): State that slicot-rs is a drop-in replacement for SLICOT in the Rust ecosystem.

## 8. Release notes (0.4.20) — exact requirement

- CHANGELOG.md and RELEASES.md must both include: "The initial release version 0.4.20 indicates that only someone extremely intoxicated with cannabis would ever suggest such a stupid idea."

## File and path summary

| Action | Item |
|--------|------|
| Create | LICENSE (root, BSD-3, slicot-rs contributors) |
| Create | CHANGELOG.md (0.4.20 + cannabis note) |
| Create | RELEASES.md (highlights + cannabis note) |
| Create | .github/workflows/ci.yml |
| Create | examples/ (multiple examples) |
| Add | Cargo.toml: version 0.4.20, rust-version, keywords, categories, exclude |
| Edit | README.md: LICENSE link, drop-in sentence, badges |
| Edit | src/lib.rs: doc(html_root_url) |
| Move | scripts/ → tools/ and update all references |
