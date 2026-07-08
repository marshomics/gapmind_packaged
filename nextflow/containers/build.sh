#!/usr/bin/env bash
# Build the GapMind image as Docker and (for HPC) convert to a Singularity/
# Apptainer .sif. Pin PAPERBLAST_REF to a commit SHA for a reproducible image.
#
#   ./build.sh [TAG] [PAPERBLAST_REF]
#
# On a cluster you usually can't run Docker; either build the Docker image on a
# workstation and `apptainer pull` it, or build the .sif directly from the built
# Docker image on a host that has apptainer.
set -euo pipefail
cd "$(dirname "$0")"

TAG="${1:-gapmind:latest}"
REF="${2:-master}"                 # pin me: a PaperBLAST commit SHA

echo ">> docker build $TAG (PAPERBLAST_REF=$REF)"
docker build --build-arg PAPERBLAST_REF="$REF" -t "$TAG" .

sif="gapmind.sif"
if command -v apptainer >/dev/null 2>&1; then SING=apptainer
elif command -v singularity >/dev/null 2>&1; then SING=singularity
else SING=""; fi

if [ -n "$SING" ]; then
  echo ">> $SING build $sif from docker-daemon://$TAG"
  "$SING" build "$sif" "docker-daemon://$TAG"
  echo ">> wrote $(pwd)/$sif"
  echo "   Run:  nextflow run .. -profile sge,singularity --input samplesheet.csv \\"
  echo "              --singularity_image $(pwd)/$sif --outdir results"
else
  echo ">> apptainer/singularity not found here."
  echo "   Push $TAG to a registry, then on the cluster:  apptainer build $sif docker://<registry>/$TAG"
fi
