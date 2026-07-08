#!/usr/bin/env bash
# Merge per-batch outputs into combined tables, one per set, each row prefixed
# with the genome name (joined from that batch's orgs.org via orgId). Also emits
# a combined orgs.tsv.
#   merge_tables.sh <out_dir> "<sets>" <keep_cand 0|1> <batch_dir> [batch_dir ...]
set -euo pipefail

out="${1:?}"; sets="${2:?}"; keep_cand="${3:-1}"; shift 3
dirs=("$@")
mkdir -p "$out"

# Combined orgs table (name map + proteome sizes).
orgs="$out/orgs.tsv"; : > "$orgs.tmp"; ow=0
for d in "${dirs[@]}"; do
  [ -s "$d/orgs.org" ] || continue
  if [ "$ow" = 0 ]; then cat "$d/orgs.org" >> "$orgs.tmp"; ow=1
  else tail -n +2 "$d/orgs.org" >> "$orgs.tmp"; fi
done
[ "$ow" = 1 ] && mv "$orgs.tmp" "$orgs" || rm -f "$orgs.tmp"

merge_one() {  # <set> <suffix>
  local s="$1" suf="$2" o="$out/$s.sum.$suf" tmp="$out/.$s.$suf.tmp" w=0 found=0
  : > "$tmp"
  for d in "${dirs[@]}"; do
    local f="$d/$s.sum.$suf" org="$d/orgs.org"
    [ -s "$f" ] && [ -s "$org" ] || continue
    found=1
    awk -v WH="$w" '
      BEGIN { FS = OFS = "\t" }
      FNR == NR { if (FNR > 1) name[$1] = $4; next }
      FNR == 1 { if (WH == 0) print "name", $0; next }
      { print (($1 in name) ? name[$1] : $1), $0 }
    ' "$org" "$f" >> "$tmp"
    w=1
  done
  [ "$found" = 1 ] && mv "$tmp" "$o" || rm -f "$tmp"
}

for s in $sets; do
  merge_one "$s" rules
  merge_one "$s" steps
  merge_one "$s" warn
  [ "$keep_cand" = 1 ] && merge_one "$s" cand
  [ "$s" = aa ] && merge_one "$s" knownsim
done
echo ">> merged tables in $out"
