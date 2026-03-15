#!/bin/sh
# Build SLICOT Fortran: lpkaux.a -> slicot.a -> example drivers.
# Run from project root. Prerequisites: gfortran, system OpenBLAS.
# See docs/FORTRAN_BUILD.md. Usage: ./scripts/slicot-fortran/build_fortran.sh [targets...]

set -e
cd "$(dirname "$0")/../.."
SLICOT_REF="${SLICOT_REF:-SLICOT-Reference}"
if [ ! -d "$SLICOT_REF" ]; then
	echo "Error: $SLICOT_REF not found. Set SLICOT_REF or clone SLICOT-Reference." >&2
	exit 1
fi
exec make -C "$SLICOT_REF" -f "$(pwd)/scripts/slicot-fortran/Makefile" "$@"
