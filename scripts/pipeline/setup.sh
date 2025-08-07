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

# ---------- MVP defaults (can still be overridden via CLI) ----------
PCS=2            # only the first 2 principal components
MAX_SNPS=5000    # 10× smaller than full run
THIN_KB=100      # looser LD pruning
MAF=0.05
MAX_MISSING=0.05
RANDOM_SEED=1    # deterministic pilot

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

exec python3 /app/scripts/fit_pca_1kg.py \
  --vcf-pattern "${VCF_PATTERN}" \
  --labels      "${LABELS_PATH}" \
  --out         "${OUT_DIR}" \
  --pcs         "${PCS}" \
  --max-snps    "${MAX_SNPS}" \
  --thin-kb     "${THIN_KB}" \
  --maf         "${MAF}" \
  --max-missing "${MAX_MISSING}" \
  --random-seed "${RANDOM_SEED}" \
  "$@"
