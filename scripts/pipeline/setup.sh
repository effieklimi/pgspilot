#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# scripts/pipeline/setup.sh
# Dual-mode wrapper that:
#   • On the HOST: mounts the project into the Docker image and re-executes
#     itself inside the container.
#   • In the CONTAINER: runs fit_pca_1kg.py with sensible MVP defaults.
###############################################################################

# ---------- Fixed paths inside the container ----------
VCF_PATTERN="/app/genome_data/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
LABELS_PATH="/app/scripts/helpers/integrated_call_samples_v3.20130502.ALL.panel"
OUT_DIR="/app/pca_model"
SITES_PATH="/app/pca_model/pca_sites.b38.tsv"
PRUNED_PATH="/app/pca_model/pca_sites.b38.tsv"   # marker file for pruned panel

# ---------- MVP defaults (can still be overridden via CLI) ----------
PCS=4             # only the first 2 principal components
RANDOM_SEED=42    # deterministic pilot

###############################################################################
### MODE 1: HOST (wrapper) ####################################################
###############################################################################
if [ -z "${INSIDE_DOCKER:-}" ]; then
  echo "==> [HOST] Preparing PCA model setup..."

  # Detect project root (assumes this script lives in scripts/)
  THIS_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
  PROJECT_ROOT=$(cd -- "${THIS_SCRIPT_DIR}/../../" &> /dev/null && pwd)

  # Forward any CLI overrides to the inner call
  EXTRA_ARGS=("$@")

  echo "==> [HOST] Launching Docker container for PCA setup..."
  docker run --rm \
    -m 48g --memory-swap -1 \
    -e INSIDE_DOCKER=1 \
    -v "${PROJECT_ROOT}:/app" \
    pgspilot \
    /app/scripts/pipeline/setup.sh "${EXTRA_ARGS[@]}"

  echo "✓ [HOST] PCA setup finished."
  exit 0
fi

###############################################################################
### MODE 2: CONTAINER (worker) ###############################################
###############################################################################
echo "==> [CONTAINER] Executing PCA setup…"
bash "/app/scripts/pipeline/pruned_panel.sh"

echo "==> [CONTAINER] Running per-ancestry MAF script..."
bash "/app/scripts/pipeline/per_ancestry_maf.sh"

echo "==> [CONTAINER] Running PCA setup..."
exec python3 -u /app/scripts/analyses/fit_pca_1kg.py \
  --vcf-pattern "${VCF_PATTERN}" \
  --labels      "${LABELS_PATH}" \
  --out         "${OUT_DIR}" \
  --pcs         "${PCS}" \
  --sites       "${SITES_PATH}" \
  --random-seed "${RANDOM_SEED}" \
  "$@"
