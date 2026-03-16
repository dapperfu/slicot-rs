#!/usr/bin/env bash
# Validate that every SLICOT routine marked "done" in docs/SLICOT_MAPPING.md
# is (1) pure Rust (no FFI) and (2) has at least one #[cfg(test)] with #[test].
# Run from project root: ./scripts/validate_slicot_done.sh
# Exit 0 if all pass; non-zero otherwise.

set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
MAPPING="${ROOT}/docs/SLICOT_MAPPING.md"
FAIL=0

echo "=== SLICOT done routines: pure Rust and tested ==="

# 1. Global FFI check in src/
if grep -rE 'extern\s|\bffi\b|libslicot|\.so\b|cdylib|staticlib' --include='*.rs' src/ 2>/dev/null; then
    echo "FAIL: FFI or foreign linkage found in src/"
    FAIL=1
fi

# 2. Per-routine: file exists, no FFI in file, has #[cfg(test)] and #[test]
while read -r mod fn; do
    f="src/${mod}/${fn}.rs"
    if [[ ! -f "$f" ]]; then
        echo "FAIL: missing $f"
        FAIL=1
    elif grep -qE 'extern\s|\bffi\b|libslicot|\.so\b' "$f" 2>/dev/null; then
        echo "FAIL: $f uses FFI"
        FAIL=1
    elif ! grep -q '#\[cfg(test)\]' "$f" 2>/dev/null || ! grep -q '#\[test\]' "$f" 2>/dev/null; then
        echo "FAIL: $f lacks #[cfg(test)] or #[test]"
        FAIL=1
    fi
done < <(awk -F'|' '/\| done \|/ { gsub(/^ +| +$/,"",$3); gsub(/^ +| +$/,"",$4); if ($3!="" && $4!="") print $3, $4 }' "$MAPPING")

if [[ $FAIL -eq 0 ]]; then
    echo "OK: All done routines are pure Rust and have tests."
else
    echo "Validation failed."
fi
exit $FAIL
