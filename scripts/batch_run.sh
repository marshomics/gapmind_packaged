#!/usr/bin/env bash
# Analyze ONE batch of genomes. Called by each SGE array task as
#   scripts/batch_run.sh $SGE_TASK_ID
# (also runnable directly for calibration). Idempotent: a completed batch has a
# .done marker and is skipped. Results per genome are keyed by GapMind's internal
# orgId; orgs.org in the batch dir maps orgId -> the name from your manifest.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

id="${1:?Usage: batch_run.sh <batch-number>}"
printf -v bid "%05d" "$id"
orgfile="$BATCH_DIR/batches/$bid.orgfile"
bdir="$BATCH_DIR/results/batch_$bid"
[ -s "$orgfile" ] || { echo "ERROR: no orgfile $orgfile" >&2; exit 1; }

if [ -f "$bdir/.done" ]; then
  echo ">> batch $bid already complete; skipping"
  exit 0
fi

# Threads: SGE sets NSLOTS to the slots granted by the PE; fall back to SGE_SLOTS.
NS="${NSLOTS:-$SGE_SLOTS}"
export MC_CORES="$NS"

cd "$CODE_DIR"                       # so bin/* and the default tmp/path.<set> resolve
rm -rf "$bdir"; mkdir -p "$bdir/cache"

echo ">> batch $bid: buildorgs ($(wc -l < "$orgfile") genomes, $NS threads)"
bin/buildorgs.pl -out "$bdir/orgs" -orgfile "$orgfile" -cache "$bdir/cache"

if [ "$SEARCH_TOOL" != "usearch" ]; then
  bin/diamond makedb --quiet --in "$bdir/orgs.faa" -d "$bdir/orgs.faa.dmnd"
fi

for set in $SETS; do
  echo ">> batch $bid: analyzing set $set"
  if [ "$SEARCH_TOOL" = "usearch" ]; then
    bin/gapsearch.pl    -orgs "$bdir/orgs" -set "$set" -out "$bdir/$set.hits" -nCPU "$NS"
    bin/gaprevsearch.pl -orgs "$bdir/orgs" -hits "$bdir/$set.hits" \
        -curated "$CODE_DIR/tmp/path.$set/curated.faa.udb" -out "$bdir/$set.revhits" -nCPU "$NS"
  else
    bin/gapsearch.pl    -diamond -orgs "$bdir/orgs" -set "$set" -out "$bdir/$set.hits" -nCPU "$NS"
    bin/gaprevsearch.pl -diamond -orgs "$bdir/orgs" -hits "$bdir/$set.hits" \
        -curated "$CODE_DIR/tmp/path.$set/curated.faa.dmnd" -out "$bdir/$set.revhits" -nCPU "$NS"
  fi

  bin/gapsummary.pl -orgs "$bdir/orgs" -set "$set" \
      -hits "$bdir/$set.hits" -revhits "$bdir/$set.revhits" -out "$bdir/$set.sum"

  # Dependency warnings (reads the shared steps.db via the results/path.<set> symlink).
  bin/checkGapRequirements.pl -results "$BATCH_DIR/results" -org "batch_$bid" \
      -set "$set" -out "$bdir/$set.sum.warn" || echo "   note: checkGapRequirements failed for $set"

  # Known-gap comparison (amino-acid set only, if the marker file exists). Uses usearch.
  if [ "$set" = "aa" ] && [ -s gaps/aa/aa.known.gaps.markers.faa ]; then
    bin/orgsVsMarkers.pl -orgs "$bdir/orgs" -vs gaps/aa/aa.known.gaps.markers.faa \
        -out "$bdir/$set.sum.knownsim" -nCPU "$NS" || echo "   note: orgsVsMarkers failed for $set"
  fi
done

# Drop the big intermediates; keep orgs.org (name map) and the *.sum.* tables.
rm -f "$bdir/orgs.faa" "$bdir/orgs.faa.dmnd" "$bdir"/*.hits "$bdir"/*.revhits
rm -rf "$bdir/cache"
: > "$bdir/.done"
echo ">> batch $bid: done ($bdir)"
