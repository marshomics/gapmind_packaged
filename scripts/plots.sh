#!/usr/bin/env bash
# Generate summary / QC / biological figures from GapMind output.
#   scripts/plots.sh [tables_dir] [orgs_table] [out_dir]
# Defaults to the merged batch output. For a single genome from `make analyze`:
#   scripts/plots.sh PaperBLAST/tmp/myorg PaperBLAST/tmp/myorg/orgs.org PaperBLAST/tmp/myorg/plots
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

tables="${1:-$BATCH_DIR/merged}"
orgs="${2:-$BATCH_DIR/merged/orgs.tsv}"
out="${3:-$BATCH_DIR/plots}"

[ -s "$orgs" ] || { echo "ERROR: orgs table not found: $orgs (run 'make batch-merge' first)" >&2; exit 1; }
python "$here/scripts/make_plots.py" \
  --tables "$tables" --orgs "$orgs" --sets "$SETS" --code-dir "$CODE_DIR" --out "$out"
echo ">> figures written to $out"
