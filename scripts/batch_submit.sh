#!/usr/bin/env bash
# Submit the SGE array job (one task per batch) plus a merge job held until the
# array finishes. Run after scripts/batch_prepare.sh. Reads sizing from config.sh
# (SGE_PE, SGE_SLOTS, SGE_TC, SGE_QUEUE, SGE_H_VMEM, SGE_H_RT).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$here/config.sh"

command -v qsub >/dev/null || { echo "ERROR: qsub not found (not on an SGE submit host?)" >&2; exit 1; }
[ -s "$BATCH_DIR/nbatches.txt" ] || { echo "ERROR: run 'make batch-prepare' first (no $BATCH_DIR/nbatches.txt)" >&2; exit 1; }
nb="$(cat "$BATCH_DIR/nbatches.txt")"
[ "$nb" -ge 1 ] 2>/dev/null || { echo "ERROR: bad batch count '$nb'" >&2; exit 1; }

conda_base="$(conda info --base)"
envv="CONDA_BASE=$conda_base,CONDA_ENV=$CONDA_ENV,PIPELINE_DIR=$here"
sge="$here/scripts/gapmind_job.sge"
logs="$BATCH_DIR/logs"

# Options common to both jobs. -pe is added only for multi-slot tasks; with
# SGE_SLOTS=1 no parallel environment is requested (most robust).
common=(-o "$logs" -j y -cwd -v "$envv")
[ -n "$SGE_QUEUE" ] && common+=(-q "$SGE_QUEUE")
pe=()
if [ "${SGE_SLOTS:-1}" -gt 1 ] 2>/dev/null; then
  pe=(-pe "$SGE_PE" "$SGE_SLOTS")
  echo ">> $nb array tasks: -pe $SGE_PE $SGE_SLOTS, up to $SGE_TC at once (<= $((SGE_TC*SGE_SLOTS)) slots), h_vmem=$SGE_H_VMEM h_rt=$SGE_H_RT"
else
  echo ">> $nb array tasks: 1 slot each (no PE), up to $SGE_TC at once (<= $SGE_TC slots), h_vmem=$SGE_H_VMEM h_rt=$SGE_H_RT"
fi

# Smoke test: TEST_TASKS=N submits only tasks 1..N (scattered across whatever
# nodes SGE picks) and skips the merge -- use it to confirm cross-node execution
# before committing the full run.
range="1-$nb"
if [ -n "${TEST_TASKS:-}" ]; then
  range="1-$TEST_TASKS"
  echo ">> SMOKE TEST: submitting only tasks $range (no merge). Rerun without TEST_TASKS for the full run."
fi

arr_args=(-terse -N gapmind -t "$range" -tc "$SGE_TC")
[ "${#pe[@]}" -gt 0 ] && arr_args+=("${pe[@]}")
arr_args+=(-l h_vmem="$SGE_H_VMEM" -l h_rt="$SGE_H_RT" "${common[@]}" "$sge")
arr_jid="$(qsub "${arr_args[@]}" | head -1 | cut -d. -f1)"
echo ">> array job id: $arr_jid"

if [ -n "${TEST_TASKS:-}" ]; then
  echo ">> When it finishes, check $logs/ and $BATCH_DIR/results/batch_00001/ (should contain a .done and *.sum.* files)."
  exit 0
fi

# Merge job: single slot, held until every array task finishes.
merge_args=(-terse -N gapmind_merge -hold_jid "$arr_jid"
            -l h_vmem=4G -l h_rt="$SGE_H_RT" "${common[@]}" "$sge")
merge_jid="$(qsub "${merge_args[@]}" | head -1 | cut -d. -f1)"
echo ">> merge job id: $merge_jid (runs after the array completes)"
echo
echo "Track:   qstat -u \$USER        (or: make batch-status)"
echo "Logs:    $logs/"
echo "Results: $BATCH_DIR/merged/ once the merge job finishes"
echo "Re-run failed tasks: fix the cause, then 'make batch-submit' again -- finished"
echo "batches (with a .done marker) are skipped, so only the incomplete ones rerun."
