#!/usr/bin/env bash
# Generate docs/FEATURES.md from docs/SLICOT_MAPPING.md.
# Summary table by Rust module and full table with Implemented/Unimplemented status.
# Run from project root: ./scripts/gen_features_table.sh

set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MAPPING="${ROOT}/docs/SLICOT_MAPPING.md"
OUT="${ROOT}/docs/FEATURES.md"

# Parse table: skip header/separator, output "module status" per data row
awk -F'|' '
/^\| [A-Z0-9_]+ \| [a-z0-9]+ \|/ {
  gsub(/^ +| +$/, "", $2);  # SLICOT
  gsub(/^ +| +$/, "", $3);  # Rust module
  gsub(/^ +| +$/, "", $4);  # Rust function
  gsub(/^ +| +$/, "", $5); # Status
  if ($2 == "" || $3 == "") next
  mod = $3
  status = ($5 == "done" ? "done" : "unimpl")
  slicot[NR] = $2
  rust_mod[NR] = mod
  rust_fn[NR] = $4
  status_line[NR] = status
  done_count[mod] += (status == "done" ? 1 : 0)
  unimpl_count[mod] += (status == "unimpl" ? 1 : 0)
  modules[mod] = 1
  order[NR] = NR
}
END {
  for (m in modules) all_mods[++k] = m
  n = asort(all_mods)
  for (i = 1; i <= n; i++) {
    m = all_mods[i]
    d = done_count[m] + 0
    u = unimpl_count[m] + 0
    t = d + u
    pct = (t > 0) ? sprintf("%.0f", 100 * d / t) : "0"
    print m, d, u, t, pct
  }
}
' "$MAPPING" > "${ROOT}/.features_summary.txt"

# Build summary table (header + rows sorted by module name)
{
  echo "# SLICOT-rs feature status"
  echo ""
  echo "Implemented vs unimplemented SLICOT routines. The authoritative per-routine list is [docs/SLICOT_MAPPING.md](SLICOT_MAPPING.md)."
  echo ""
  echo "## Summary by module"
  echo ""
  echo "| Module | Implemented | Unimplemented | Total | % done |"
  echo "|--------|-------------|----------------|-------|--------|"
  sort -t' ' -k1,1 .features_summary.txt | while read -r mod impl unimpl total pct; do
    printf "| %s | %s | %s | %s | %s |\n" "$mod" "$impl" "$unimpl" "$total" "$pct"
  done
  echo ""
  impl_total=$(awk '{s+=$2} END {print s+0}' .features_summary.txt)
  unimpl_total=$(awk '{s+=$3} END {print s+0}' .features_summary.txt)
  total_total=$((impl_total + unimpl_total))
  pct_total=$(( total_total > 0 ? impl_total * 100 / total_total : 0 ))
  echo "| **Total** | **$impl_total** | **$unimpl_total** | **$total_total** | **$pct_total** |"
  echo ""
  echo "## Full table (all routines)"
  echo ""
  echo "| SLICOT | Rust module | Rust function | Status |"
  echo "|--------|-------------|---------------|--------|"
} > "$OUT"

# Append full table rows (same order as mapping)
awk -F'|' '
/^\| [A-Z0-9_]+ \| [a-z0-9]+ \|/ {
  gsub(/^ +| +$/, "", $2)
  gsub(/^ +| +$/, "", $3)
  gsub(/^ +| +$/, "", $4)
  gsub(/^ +| +$/, "", $5)
  if ($2 == "" || $3 == "") next
  status = ($5 == "done" ? "Implemented" : "Unimplemented")
  print "| " $2 " | " $3 " | " $4 " | " status " |"
}
' "$MAPPING" >> "$OUT"

rm -f "${ROOT}/.features_summary.txt"
echo "Generated $OUT"
