#!/bin/sh
# Run full Fortran 1:1 validation: build Fortran if needed, run Rust validation runner,
# write validation/*.md (module docs, FAILURES.md, README summary).
# Run from project root. Prerequisites: gfortran, OpenBLAS (for full Fortran build).
# See validation/README.md and docs/FORTRAN_BUILD.md.

set -e
cd "$(dirname "$0")/../.."
SLICOT_REF="${SLICOT_REF:-SLICOT-Reference}"
EXAMPLES_DIR="${EXAMPLES_DIR:-$SLICOT_REF/examples}"

if [ ! -d "$SLICOT_REF" ]; then
	echo "Error: $SLICOT_REF not found. Set SLICOT_REF or clone SLICOT-Reference." >&2
	exit 1
fi

# Build Fortran (lpkaux, slicot, example drivers) if not already built
if [ ! -f "$EXAMPLES_DIR/TAB01ND" ] && [ ! -f "$EXAMPLES_DIR/tab01nd" ]; then
	echo "Building Fortran SLICOT and examples..."
	./tools/slicot-fortran/build_fortran.sh || true
fi

export SLICOT_EXAMPLES_DIR="$EXAMPLES_DIR"
cargo test --test fortran_validation -- --nocapture 2>&1 | tee validation_run.log
echo "Validation run complete. See validation/*.md and validation_run.log"
