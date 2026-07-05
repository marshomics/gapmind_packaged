#!/usr/bin/env bash
# Build genome x pathway presence/absence (0/1) + confidence (2/1/0) matrices.
#   scripts/presence.sh [tables_dir] [orgs_table] [out_dir]
# Defaults to the merged batch output, writing the matrices alongside it.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

tables="${1:-$BATCH_DIR/merged}"
orgs="${2:-$BATCH_DIR/merged/orgs.tsv}"
out="${3:-$BATCH_DIR/merged}"

[ -s "$orgs" ] || { echo "ERROR: orgs table not found: $orgs (run 'make batch-merge' first)" >&2; exit 1; }
python "$here/scripts/presence_absence.py" \
  --tables "$tables" --orgs "$orgs" --sets "$SETS" --code-dir "$CODE_DIR" \
  --mode "$PA_MODE" --out "$out"
echo ">> presence/confidence matrices written to $out"
