#!/usr/bin/env bash
# Pin containers/environment.yml to an exact, hash-verified lock so the image is
# reproducible build-to-build. This sandbox has no conda; run this on any host
# that does have it (a workstation, a CI runner, or inside the built image).
#
#   ./lock.sh                       # linux-64 by default
#   PLATFORMS="linux-64 linux-aarch64 osx-arm64" ./lock.sh
#
# Two methods, preferred first:
#   1) conda-lock  -> conda-lock.yml   (multi-platform, per-package hashes)
#   2) conda       -> conda-<plat>.lock via --explicit (this platform only)
#
# To consume the lock in the image, replace the env step in the Dockerfile with:
#   COPY environment.yml conda-lock.yml /tmp/
#   RUN conda install -n base -c conda-forge conda-lock \
#    && conda-lock install -n base /tmp/conda-lock.yml && conda clean -afy
set -euo pipefail
cd "$(dirname "$0")"

ENV_FILE=environment.yml
PLATFORMS="${PLATFORMS:-linux-64}"     # space-separated; add arches as needed

if [ ! -f "$ENV_FILE" ]; then
  echo "!! $ENV_FILE not found (run from containers/)" >&2
  exit 1
fi

if command -v conda-lock >/dev/null 2>&1; then
  echo ">> conda-lock for: $PLATFORMS"
  plat_args=()
  for p in $PLATFORMS; do plat_args+=(-p "$p"); done
  conda-lock lock -f "$ENV_FILE" "${plat_args[@]}" --lockfile conda-lock.yml
  # Also render human-diffable explicit files (one per platform).
  conda-lock render -k explicit conda-lock.yml || true
  echo ">> wrote conda-lock.yml"
elif command -v conda >/dev/null 2>&1; then
  echo ">> conda-lock not found; explicit lock for the CURRENT platform only"
  tmp="_gapmind_lock_$$"
  cleanup() { conda env remove -n "$tmp" -y >/dev/null 2>&1 || true; }
  trap cleanup EXIT
  conda env create -n "$tmp" -f "$ENV_FILE"
  out="conda-$(conda info --json | python -c 'import json,sys;print(json.load(sys.stdin)["platform"])').lock"
  conda list -n "$tmp" --explicit --md5 > "$out"
  echo ">> wrote $out (single platform)"
else
  echo "!! neither conda-lock nor conda is on PATH. Install one first:" >&2
  echo "     pipx install conda-lock        # or" >&2
  echo "     conda install -n base -c conda-forge conda-lock" >&2
  exit 1
fi
