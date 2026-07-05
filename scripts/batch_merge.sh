#!/usr/bin/env bash
# Merge per-batch GapMind outputs into combined tables, one per set, with each
# row labelled by the genome "name" from your manifest (joined from the batch's
# orgs.org via orgId). Reports which batches are complete vs missing. Safe to run
# repeatedly; merges whatever is done so far.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

[ -s "$BATCH_DIR/nbatches.txt" ] || { echo "ERROR: no $BATCH_DIR/nbatches.txt (run batch-prepare)" >&2; exit 1; }
nb="$(cat "$BATCH_DIR/nbatches.txt")"
mkdir -p "$BATCH_DIR/merged"

# ---- completeness report ----
done=0; : > "$BATCH_DIR/incomplete.txt"
for i in $(seq 1 "$nb"); do
  printf -v bid "%05d" "$i"
  if [ -f "$BATCH_DIR/results/batch_$bid/.done" ]; then done=$((done+1)); else echo "$i" >> "$BATCH_DIR/incomplete.txt"; fi
done
echo ">> $done / $nb batches complete"
if [ "$done" -ne "$nb" ]; then
  echo ">> $((nb-done)) incomplete (ids in $BATCH_DIR/incomplete.txt); merging the finished ones."
  echo "   Re-run 'make batch-submit' to finish them, then merge again."
fi

# ---- merge one table type across all completed batches ----
# Prepends a "name" column looked up from each batch's orgs.org (orgId -> genomeName).
merge_table() {
  local set="$1" suf="$2"
  local out="$BATCH_DIR/merged/$set.sum.$suf"
  local tmp="$out.tmp.$$" wrote=0 found=0
  : > "$tmp"
  for d in "$BATCH_DIR"/results/batch_*/; do
    [ -f "${d}.done" ] || continue
    local f="${d}${set}.sum.${suf}" org="${d}orgs.org"
    [ -s "$f" ] && [ -s "$org" ] || continue
    found=1
    awk -v WH="$wrote" '
      BEGIN { FS = OFS = "\t" }
      FNR == NR { if (FNR > 1) name[$1] = $4; next }      # orgs.org: orgId=$1 name=$4
      FNR == 1 { if (WH == 0) print "name", $0; next }    # original header once
      { print (($1 in name) ? name[$1] : $1), $0 }
    ' "$org" "$f" >> "$tmp"
    wrote=1
  done
  if [ "$found" = 1 ]; then mv "$tmp" "$out"; echo "   wrote $out"; else rm -f "$tmp"; fi
}

for set in $SETS; do
  echo ">> merging set $set"
  merge_table "$set" rules
  merge_table "$set" steps
  merge_table "$set" warn
  [ "$KEEP_CAND" = 1 ] && merge_table "$set" cand
  [ "$set" = "aa" ] && merge_table "$set" knownsim
done

# Combined orgs table (name map + proteome sizes) -- used by the plots and handy
# on its own. genomeName holds the manifest name; nProteins is the proteome size.
orgs_out="$BATCH_DIR/merged/orgs.tsv"; : > "$orgs_out.tmp"; owrote=0
for d in "$BATCH_DIR"/results/batch_*/; do
  [ -f "${d}.done" ] && [ -s "${d}orgs.org" ] || continue
  if [ "$owrote" = 0 ]; then cat "${d}orgs.org" >> "$orgs_out.tmp"; owrote=1
  else tail -n +2 "${d}orgs.org" >> "$orgs_out.tmp"; fi
done
[ "$owrote" = 1 ] && mv "$orgs_out.tmp" "$orgs_out" || rm -f "$orgs_out.tmp"

echo ">> Merged tables in $BATCH_DIR/merged/"
echo "   Start with $BATCH_DIR/merged/<set>.sum.rules (one row per genome x pathway;"
echo "   nHi/nMed/nLo columns), labelled by your manifest name."

# Summary / QC / biological figures (best-effort; needs matplotlib in the env).
if [ -s "$orgs_out" ] && python -c 'import matplotlib' >/dev/null 2>&1; then
  echo ">> generating plots"
  python "$here/scripts/make_plots.py" --tables "$BATCH_DIR/merged" --orgs "$orgs_out" \
      --sets "$SETS" --code-dir "$CODE_DIR" --out "$BATCH_DIR/plots" \
    || echo "   note: plotting failed; rerun with 'make plots'"
else
  echo ">> skipping plots (matplotlib not in env yet); run 'make env' then 'make plots'"
fi
