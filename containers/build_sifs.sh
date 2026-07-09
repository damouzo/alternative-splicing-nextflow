#!/usr/bin/env bash
# Build ISAR, report, and sashimi containers as Apptainer SIF files.
#
# Usage:
#   sbatch containers/build_sifs.sh                         # default output dir
#   sbatch containers/build_sifs.sh /custom/path/containers # custom output dir
#
# This is only needed when running without internet access to a public registry,
# or when you want a local copy of the containers.
# Normal usage: images are pulled automatically from ghcr.io on first pipeline run.
#
# After building, add to your params.yaml:
#   isar_container:   <CONTAINERS_DIR>/isar-1.0.0.sif
#   report_container:  <CONTAINERS_DIR>/report-1.0.0.sif
#   sashimi_container: <CONTAINERS_DIR>/sashimi-1.0.0.sif

#SBATCH --job-name=build_sifs
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=containers/build_sifs_%j.log

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINERS_DIR="${1:-/data/BCI-KRP/containers}"

mkdir -p "${CONTAINERS_DIR}"

echo "=== Building ISAR SIF ==="
echo "Output: ${CONTAINERS_DIR}/isar-1.0.0.sif"
apptainer build "${CONTAINERS_DIR}/isar-1.0.0.sif" "${SCRIPT_DIR}/isar/isar.def"

echo ""
echo "=== Building report SIF ==="
echo "Output: ${CONTAINERS_DIR}/report-1.0.0.sif"
apptainer build "${CONTAINERS_DIR}/report-1.0.0.sif" "${SCRIPT_DIR}/report/report.def"

echo ""
echo "=== Building sashimi SIF ==="
echo "Output: ${CONTAINERS_DIR}/sashimi-1.0.0.sif"
apptainer build "${CONTAINERS_DIR}/sashimi-1.0.0.sif" "${SCRIPT_DIR}/sashimi/sashimi.def"

echo ""
echo "=== Verification ==="
apptainer exec "${CONTAINERS_DIR}/isar-1.0.0.sif" \
    Rscript -e "library(IsoformSwitchAnalyzeR); cat('ISAR OK\n')"
apptainer exec "${CONTAINERS_DIR}/report-1.0.0.sif" \
    Rscript -e "library(DT); library(plotly); cat('Report OK\n')"
apptainer exec "${CONTAINERS_DIR}/sashimi-1.0.0.sif" \
    /bin/bash -lc "command -v rmats2sashimiplot >/dev/null && command -v samtools >/dev/null && echo 'Sashimi OK'"

echo ""
echo "=== Done ==="
echo "Add to your params.yaml:"
echo "  isar_container:   ${CONTAINERS_DIR}/isar-1.0.0.sif"
echo "  report_container: ${CONTAINERS_DIR}/report-1.0.0.sif"
echo "  sashimi_container: ${CONTAINERS_DIR}/sashimi-1.0.0.sif"
