#!/usr/bin/env bash
# Unified build and validation runner: build FORTRAN, build Rust (tests), run 1:1
# Rust-vs-FORTRAN validation. Fails if any validated routine does not match FORTRAN.
# Run from project root. Prerequisites: gfortran, OpenBLAS (for Fortran). See docs/FORTRAN_BUILD.md.
#
# Usage:
#   ./run.sh                    # full: Fortran + Rust + validation
#   ./run.sh --no-fortran       # skip Fortran build (e.g. already built)
#   ./run.sh --done-check       # also run validate_slicot_done.sh (pure Rust + has tests)
#   ./run.sh --all-targets      # also run cargo build --release --all-targets (lib, bins, benches)

set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

SLICOT_REF="${SLICOT_REF:-SLICOT-Reference}"
export SLICOT_EXAMPLES_DIR="${SLICOT_EXAMPLES_DIR:-$ROOT/$SLICOT_REF/examples}"

SKIP_FORTRAN=false
DONE_CHECK=false
ALL_TARGETS=false
for arg in "$@"; do
    case "$arg" in
        --no-fortran)  SKIP_FORTRAN=true ;;
        --done-check)  DONE_CHECK=true ;;
        --all-targets) ALL_TARGETS=true ;;
        *) echo "Unknown option: $arg" >&2; echo "Usage: $0 [--no-fortran] [--done-check] [--all-targets]" >&2; exit 1 ;;
    esac
done

echo "=== 1. Build FORTRAN (SLICOT reference and example drivers) ==="
if [ "$SKIP_FORTRAN" = true ]; then
    echo "Skipping Fortran build (--no-fortran)."
else
    if [ ! -d "$SLICOT_REF" ]; then
        echo "Error: $SLICOT_REF not found. Set SLICOT_REF or clone SLICOT-Reference." >&2
        echo "See docs/FORTRAN_BUILD.md." >&2
        exit 1
    fi
    if ! ./tools/slicot-fortran/build_fortran.sh; then
        echo "Error: Fortran build failed. Need gfortran and OpenBLAS. See docs/FORTRAN_BUILD.md." >&2
        exit 1
    fi
fi

echo "=== 2. Build Rust (library and tests for validation) ==="
cargo build --tests
if [ "$ALL_TARGETS" = true ]; then
    echo "Building release and all targets (bins, benches)..."
    cargo build --release --all-targets
fi

echo "=== 3. Run 1:1 validation (Rust vs FORTRAN) ==="
cargo test --test fortran_validation -- --nocapture
VALIDATION_EXIT=$?
if [ $VALIDATION_EXIT -ne 0 ]; then
    echo "Validation failed. See validation/FAILURES.md for failed routines." >&2
    exit $VALIDATION_EXIT
fi

if [ "$DONE_CHECK" = true ]; then
    echo "=== 4. Validate done routines (pure Rust + has tests) ==="
    ./tools/validate_slicot_done.sh
fi

echo "=== Done: Fortran built, Rust built, 1:1 validation passed. ==="
