#!/bin/sh
# Optional sanity check: build Fortran (if SLICOT-Reference present) and run a short fuzz.
# From project root: ./scripts/sanity_check_fortran_fuzz.sh
# Requires: cargo-fuzz (cargo install cargo-fuzz), and optionally gfortran + OpenBLAS for Fortran.

set -e
cd "$(dirname "$0")/.."

echo "=== Sanity check: Fortran build + fuzz ==="
if [ -d "SLICOT-Reference" ]; then
  echo "Building Fortran..."
  ./scripts/slicot-fortran/build_fortran.sh lpkaux.a || true
  echo "Copy file-I/O driver: cp scripts/slicot-fortran/drivers/TAB01ND.f SLICOT-Reference/examples/"
  echo "Then: ./scripts/slicot-fortran/build_fortran.sh"
else
  echo "SLICOT-Reference not found; skip Fortran build."
fi

echo "Running fuzz target ab01nd_compare (10 runs)..."
cargo fuzz run ab01nd_compare -- -runs=10 2>&1 || true

echo "Done. For full fuzz: cargo fuzz run ab01nd_compare"
