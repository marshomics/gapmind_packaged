#!/usr/bin/env bash
# Analyze one batch of genomes (the orgfile lists them) against the prepared
# databases. Runs inside the container; writes batch_<id>/ with the *.sum.*
# tables and orgs.org, cleaning the large intermediates.
#   run_batch.sh <id> <orgfile> <db_dir> "<sets>" <cpus> <search_tool> <knownsim0|1>
set -euo pipefail

id="${1:?}"; orgfile="${2:?}"; db="${3:?}"; sets="${4:?}"
ncpu="${5:-4}"; tool="${6:-diamond}"; knownsim="${7:-1}"
: "${GAPMIND_DIR:?GAPMIND_DIR not set (run inside the gapmind container)}"
export MC_CORES="$ncpu"

out="batch_${id}"
mkdir -p "$out"

# Build the orgs table. Cache in the (writable) work dir, not the read-only image.
buildorgs.pl -out "$out/orgs" -orgfile "$orgfile" -cache ./cache

if [ "$tool" != usearch ]; then
  diamond makedb --quiet --in "$out/orgs.faa" -d "$out/orgs.faa.dmnd"
fi

# checkGapRequirements.pl reads <results>/path.<set>/steps.db and ignores its own
# -stepsDb flag, so expose the DB via a symlink in the work dir and use -results .
for s in $sets; do ln -sfn "$db/path.$s" "path.$s"; done

for s in $sets; do
  echo ">> batch $id: set $s"
  if [ "$tool" = usearch ]; then
    gapsearch.pl    -orgs "$out/orgs" -set "$s" -dir "$db/path.$s" -out "$out/$s.hits" -nCPU "$ncpu"
    gaprevsearch.pl -orgs "$out/orgs" -hits "$out/$s.hits" \
        -curated "$db/path.$s/curated.faa.udb" -out "$out/$s.revhits" -nCPU "$ncpu"
  else
    gapsearch.pl    -diamond -orgs "$out/orgs" -set "$s" -dir "$db/path.$s" -out "$out/$s.hits" -nCPU "$ncpu"
    gaprevsearch.pl -diamond -orgs "$out/orgs" -hits "$out/$s.hits" \
        -curated "$db/path.$s/curated.faa.dmnd" -out "$out/$s.revhits" -nCPU "$ncpu"
  fi
  gapsummary.pl -orgs "$out/orgs" -set "$s" -dbDir "$db/path.$s" \
      -hits "$out/$s.hits" -revhits "$out/$s.revhits" -out "$out/$s.sum"
  checkGapRequirements.pl -results . -org "$out" -set "$s" -out "$out/$s.sum.warn" \
      || echo "   note: checkGapRequirements failed for $s"
  markers="$GAPMIND_DIR/gaps/aa/aa.known.gaps.markers.faa"
  if [ "$knownsim" = 1 ] && [ "$s" = aa ] && [ -s "$markers" ]; then
    orgsVsMarkers.pl -orgs "$out/orgs" -vs "$markers" -out "$out/$s.sum.knownsim" -nCPU "$ncpu" \
        || echo "   note: orgsVsMarkers failed for $s"
  fi
done

# Keep the tables + orgs.org; drop the big intermediates and temp symlinks.
rm -f "$out/orgs.faa" "$out/orgs.faa.dmnd" "$out"/*.hits "$out"/*.revhits
rm -rf ./cache
for s in $sets; do rm -f "path.$s"; done
echo ">> $out done"
