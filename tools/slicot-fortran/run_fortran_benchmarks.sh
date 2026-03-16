#!/usr/bin/env sh
# Build Fortran SLICOT (if needed), build and run the Fortran benchmark driver.
# Run from project root. Uses -lblas -llapack if OpenBLAS is not available.
# Output: time per call (us) for MA02ED, MA02ES, DLACPY at n=32,64,...,1024.

set -e
cd "$(dirname "$0")/../.."
SLICOT_REF="${SLICOT_REF:-SLICOT-Reference}"
BENCH_DIR="scripts/slicot-fortran/bench"

# Build slicot.a and lpkaux.a if missing
if [ ! -f "$SLICOT_REF/slicot.a" ] || [ ! -f "$SLICOT_REF/lpkaux.a" ]; then
  echo "Building Fortran SLICOT (lpkaux, slicot)..."
  ./scripts/slicot-fortran/build_fortran.sh lpkaux.a slicot
  # examples target may fail without OpenBLAS; we only need slicot.a
fi

# Build benchmark driver (prefer OpenBLAS; fallback to ref BLAS/LAPACK)
echo "Building Fortran benchmark driver..."
if ! make -C "$BENCH_DIR" SLICOT_REF="$(pwd)/$SLICOT_REF" 2>/dev/null; then
  echo "Retrying with -lblas -llapack..."
  make -C "$BENCH_DIR" SLICOT_REF="$(pwd)/$SLICOT_REF" \
    BLASLIB="-lblas" LAPACKLIB="-llapack"
fi

echo "Running Fortran benchmarks..."
"$BENCH_DIR/bench_slicot"
