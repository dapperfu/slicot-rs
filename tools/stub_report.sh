#!/usr/bin/env bash
# Produce a report of SLICOT routines that are stubbed or partially implemented.
# Scans docs/SLICOT_MAPPING.md (Implementation column) and src/**/*.rs for
# doc/source keywords: Stub, not implemented, placeholder.
# Run from project root: ./tools/stub_report.sh

set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MAPPING="${ROOT}/docs/SLICOT_MAPPING.md"
SRC="${ROOT}/src"

echo "=== SLICOT stub / partial implementation report ==="
echo ""

# 1. From SLICOT_MAPPING.md: list stub and partial rows
echo "--- From docs/SLICOT_MAPPING.md (Implementation column) ---"
awk -F'|' '
/^\| [A-Z0-9_]+ \| [a-z0-9]+ \|/ {
  gsub(/^ +| +$/, "", $2)   # SLICOT
  gsub(/^ +| +$/, "", $4)   # Rust function
  gsub(/^ +| +$/, "", $6)   # Implementation (after adding column)
  if ($2 == "" || $4 == "") next
  impl = $6
  if (impl == "stub" || impl == "partial") print impl, $2, $4
}
' "$MAPPING" | sort -k1,1 -k2,2 | while read impl slicot rustfn; do
  echo "  $impl: $slicot ($rustfn)"
done
echo ""

# 2. Grep src for Stub / not implemented / placeholder in docs and comments
echo "--- Source files with 'Stub' / 'not implemented' / 'placeholder' (doc/comments) ---"
grep -r -l -i -E 'Stub|not implemented|placeholder' "$SRC" --include='*.rs' 2>/dev/null | sort -u | while read f; do
  echo "  $f"
  grep -n -i -E 'Stub|not implemented|placeholder' "$f" 2>/dev/null | sed 's/^/    /'
done
echo ""

# 3. Count by Implementation status from mapping
echo "--- Counts (from SLICOT_MAPPING.md) ---"
awk -F'|' '
/^\| [A-Z0-9_]+ \| [a-z0-9]+ \|/ {
  gsub(/^ +| +$/, "", $6)
  if ($6 == "stub") stub++
  else if ($6 == "partial") partial++
  else if ($6 == "full") full++
}
END {
  print "  full:", full + 0
  print "  partial:", partial + 0
  print "  stub:", stub + 0
}
' "$MAPPING"
