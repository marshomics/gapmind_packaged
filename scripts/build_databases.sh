#!/usr/bin/env bash
# Phase 1 (curated database) + Phase 2 (steps database) for every set in $SETS,
# then format the curated database for the chosen search tool.
#
# Heavy and network-bound: downloads the PaperBLAST data (~1.5 GB), Swiss-Prot
# (~0.6 GB), and runs hmmsearch over every Pfam model (Phase 1). Budget several
# hours and tens of GB of disk. Phase 2 also fetches some UniProt sequences,
# so an internet connection is required throughout.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

cd "$CODE_DIR"

# Thread count for the steps that honor it (e.g. the hmmsearch in runPfamHits).
export MC_CORES="$THREADS"

# ----------------------------------------------------------------------------
# DB_SOURCE=prebuilt (default): download the maintainer's consistent curated +
# steps databases per set and format them. No Phase 1/2, no PaperBLAST data, no
# Swiss-Prot, no Pfam download -- steps.db already carries its HMMs, which
# extractHmms.pl pulls out. This matches the official, validated GapMind.
# ----------------------------------------------------------------------------
if [ "$DB_SOURCE" = "prebuilt" ]; then
  for set in $SETS; do
    echo "=================================================================="
    echo ">> Prebuilt databases: $set"
    echo "=================================================================="
    mkdir -p "tmp/path.$set"
    if [ -f "tmp/path.$set/.prebuilt" ]; then
      echo ">> already downloaded for $set; skipping (rm tmp/path.$set/.prebuilt to refetch)"
    else
      # Overwrites any from-scratch files so the curated + steps DBs stay a
      # consistent set (the official ones reference each other's IDs).
      for f in curated.faa curated.db steps.db; do
        echo ">> Downloading path.$set/$f"
        curl -fSL "$PREBUILT_BASE/path.$set/$f" -o "tmp/path.$set/$f"
      done
      : > "tmp/path.$set/.prebuilt"
    fi
    bin/extractHmms.pl "tmp/path.$set/steps.db" "tmp/path.$set"
    bin/blast/formatdb -p T -o T -i "tmp/path.$set/curated.faa"
    if [ "$SEARCH_TOOL" = "usearch" ]; then
      bin/usearch -makeudb_ublast "tmp/path.$set/curated.faa" -output "tmp/path.$set/curated.faa.udb"
    else
      bin/diamond makedb --in "tmp/path.$set/curated.faa" -d "tmp/path.$set/curated.faa.dmnd"
    fi
    echo ">> Done (prebuilt): $set"
  done
  echo ">> Prebuilt databases ready under $CODE_DIR/tmp/path.*"
  exit 0
fi

# ============================================================================
# DB_SOURCE=build: rebuild the curated + steps databases from scratch (Phases
# 1-2). Everything below runs with the repo root as the working directory
# because the Swiss-Prot parsers use `use lib "SWISS/lib"` and the scripts use
# relative tmp/ and data/ paths.
# ============================================================================

# --- shared inputs (downloaded once) --------------------------------------
mkdir -p data ind

for f in litsearch.db uniq.faa; do
  if [ ! -s "data/$f" ]; then
    echo ">> Downloading data/$f"
    curl -fSL "$PAPERBLAST_DATA_BASE/$f" -o "data/$f"
  fi
done
# GapMind's curatedFaa.pl reads uniq.faa as a plain FASTA (not via fastacmd), so
# it does not need to be formatted with formatdb. The PaperBLAST/SitesBLAST web
# CGIs would need that, but this pipeline targets the GapMind command line.

if [ ! -s ind/uniprot_sprot.dat.gz ]; then
  echo ">> Downloading Swiss-Prot"
  curl -fSL "$SPROT_URL" -o ind/uniprot_sprot.dat.gz
fi

# Characterized Swiss-Prot entries, in curated_parsed format (Phase 1 input).
if [ ! -s sprot.curated_parsed ]; then
  echo ">> Parsing Swiss-Prot (sprotCharacterized.pl)"
  gzip -dc ind/uniprot_sprot.dat.gz | perl -I SWISS/lib bin/sprotCharacterized.pl > sprot.curated_parsed
fi

# --- per-set build ---------------------------------------------------------
for set in $SETS; do
  echo "=================================================================="
  echo ">> Building set: $set"
  echo "=================================================================="

  # Phase 1: curated.faa, curated.db, curated2.faa, hetero.tab, pfam.hits.tab ...
  # curated.db is written by the final setupGaps step, so treat it as the "done"
  # marker and skip the (slow) Phase 1 rebuild when re-running.
  if [ -s "tmp/path.$set/curated.db" ]; then
    echo ">> Phase 1 already complete for $set (curated.db present); skipping setupGaps"
  else
    perl bin/setupGaps.pl -ind ind -set "$set" -data data -sprot sprot.curated_parsed
  fi

  # Phase 2: -doquery runs gapquery.pl for every pathway to build the per-pathway
  # *.query files (this also fetches a few UniProt sequences over the network),
  # then assembles steps.db. steps.db is written last, so use it as the marker.
  if [ -s "tmp/path.$set/steps.db" ]; then
    echo ">> Phase 2 already complete for $set (steps.db present); skipping buildStepsDb"
  else
    perl bin/buildStepsDb.pl -set "$set" -doquery
  fi

  # Pull the HMM models referenced by the steps DB into the working dir.
  bin/extractHmms.pl "tmp/path.$set/steps.db" "tmp/path.$set"

  # Format the curated DB: BLAST (for the web viewer) + the chosen search tool.
  bin/blast/formatdb -p T -o T -i "tmp/path.$set/curated.faa"
  if [ "$SEARCH_TOOL" = "usearch" ]; then
    bin/usearch -makeudb_ublast "tmp/path.$set/curated.faa" -output "tmp/path.$set/curated.faa.udb"
  else
    bin/diamond makedb --in "tmp/path.$set/curated.faa" -d "tmp/path.$set/curated.faa.dmnd"
  fi

  echo ">> Done set: $set  (databases in tmp/path.$set/)"
done

echo ">> All curated + steps databases built under $CODE_DIR/tmp/path.*"
