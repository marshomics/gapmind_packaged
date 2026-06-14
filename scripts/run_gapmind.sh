#!/usr/bin/env bash
# Phase 3: analyze one proteome (a protein FASTA) against the built databases.
# Runs the full command-line chain: buildorgs -> gapsearch -> gaprevsearch ->
# gapsummary -> checkGapRequirements -> orgsVsMarkers -> buildGapsDb.
#
#   run_gapmind.sh -f proteins.faa -n myorg [-s "aa carbon"] [-t 8]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

usage() {
  echo "Usage: $0 -f proteins.faa -n NAME [-s 'aa carbon'] [-t THREADS]" >&2
  exit 1
}

FAA=""; NAME=""; RUN_SETS="$SETS"; T="$THREADS"
while getopts "f:n:s:t:h" opt; do
  case "$opt" in
    f) FAA="$OPTARG" ;;
    n) NAME="$OPTARG" ;;
    s) RUN_SETS="$OPTARG" ;;
    t) T="$OPTARG" ;;
    *) usage ;;
  esac
done
[ -n "$FAA" ] && [ -n "$NAME" ] || usage
[ -s "$FAA" ] || { echo "ERROR: no such FASTA: $FAA" >&2; exit 1; }
case "$NAME" in
  *[!A-Za-z0-9_.-]*) echo "ERROR: NAME may contain only letters, digits, '_', '.', '-'" >&2; exit 1 ;;
esac

# Absolutize the FASTA path before we cd into the code base.
FAA="$(cd "$(dirname "$FAA")" && pwd)/$(basename "$FAA")"

cd "$CODE_DIR"
out="tmp/$NAME"
mkdir -p "$out"

# buildorgs.pl splits the org specifier on ":" and rejects whitespace in the file
# path, so hand it a clean, space-free relative path via a symlink.
ln -sf "$FAA" "$out/input.faa"
echo ">> buildorgs ($NAME)"
bin/buildorgs.pl -out "$out/orgs" -orgs "file:$out/input.faa:$NAME"

# gapsearch.pl -diamond requires <orgprefix>.faa.dmnd, which is set-independent,
# so build it once here rather than per set.
if [ "$SEARCH_TOOL" != "usearch" ]; then
  bin/diamond makedb --quiet --in "$out/orgs.faa" -d "$out/orgs.faa.dmnd"
fi

for set in $RUN_SETS; do
  echo "=================================================================="
  echo ">> Analyzing $NAME : set $set"
  echo "=================================================================="

  if [ "$SEARCH_TOOL" = "usearch" ]; then
    bin/gapsearch.pl    -orgs "$out/orgs" -set "$set" -out "$out/$set.hits" -nCPU "$T"
    bin/gaprevsearch.pl -orgs "$out/orgs" -hits "$out/$set.hits" \
        -curated "tmp/path.$set/curated.faa.udb" -out "$out/$set.revhits" -nCPU "$T"
  else
    bin/gapsearch.pl    -diamond -orgs "$out/orgs" -set "$set" -out "$out/$set.hits" -nCPU "$T"
    bin/gaprevsearch.pl -diamond -orgs "$out/orgs" -hits "$out/$set.hits" \
        -curated "tmp/path.$set/curated.faa.dmnd" -out "$out/$set.revhits" -nCPU "$T"
  fi

  bin/gapsummary.pl -orgs "$out/orgs" -set "$set" \
      -hits "$out/$set.hits" -revhits "$out/$set.revhits" -out "$out/$set.sum"

  # Dependency check between pathways. buildGapsDb.pl requires -requirements, so
  # this step is mandatory; it still writes a header-only file when a set has no
  # requirements, which buildGapsDb reads happily.
  bin/checkGapRequirements.pl -org "$NAME" -set "$set" -out "$out/$set.sum.warn"

  cmd=(bin/buildGapsDb.pl -gaps "$out/$set.sum"
       -requirements "$out/$set.sum.warn"
       -steps "tmp/path.$set/steps.db" -out "$out/$set.sum.db")

  # Optional comparison to organisms with known gaps (amino-acid set only, and
  # only if the marker file was produced when the steps DB was built).
  # orgsVsMarkers.pl uses usearch regardless of $SEARCH_TOOL.
  if [ "$set" = "aa" ] && [ -s gaps/aa/aa.known.gaps.markers.faa ]; then
    if bin/orgsVsMarkers.pl -orgs "$out/orgs" \
         -vs gaps/aa/aa.known.gaps.markers.faa -out "$out/$set.sum.knownsim"; then
      cmd+=(-markersim "$out/$set.sum.knownsim")
    else
      echo "   note: orgsVsMarkers failed; continuing without it"
    fi
  fi

  # Combine everything into one sqlite3 database.
  "${cmd[@]}"

  echo ">> $set results: $CODE_DIR/$out/$set.sum.{rules,steps,cand} and $set.sum.db"
done

echo ">> Analysis complete for $NAME (outputs under $CODE_DIR/$out/)"
