#!/usr/bin/env bash
# Download the prebuilt curated + steps databases per set and format them for the
# analysis. Runs inside the container (extractHmms.pl, diamond, usearch on PATH).
#   prepare_db.sh "<sets>" <out_db_dir> <prebuilt_base_url> <search_tool>
set -euo pipefail

sets="${1:?}"; db="${2:?}"; base="${3:?}"; tool="${4:-diamond}"
: "${GAPMIND_DIR:?GAPMIND_DIR not set (run inside the gapmind container)}"

mkdir -p "$db"
for s in $sets; do
  d="$db/path.$s"; mkdir -p "$d"
  for f in curated.faa curated.db steps.db; do
    echo ">> downloading path.$s/$f"
    curl -fSL "$base/path.$s/$f" -o "$d/$f"
  done
  # Pull the HMM models that steps.db references, alongside it.
  "$GAPMIND_DIR/bin/extractHmms.pl" "$d/steps.db" "$d"
  # Format the curated DB for the chosen search tool.
  if [ "$tool" = usearch ]; then
    "$GAPMIND_DIR/bin/usearch" -makeudb_ublast "$d/curated.faa" -output "$d/curated.faa.udb"
  else
    diamond makedb --quiet --in "$d/curated.faa" -d "$d/curated.faa.dmnd"
  fi
  echo ">> path.$s ready"
done
