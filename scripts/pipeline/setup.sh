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

  echo "==> [HOST] Launching Docker container for PCA setup..."
  docker run --rm \
    -m 48g --memory-swap -1 \
    -e INSIDE_DOCKER=1 \
    -v "${PROJECT_ROOT}:/app" \
    pgspilot \
    /app/scripts/pipeline/setup.sh "$@"

  echo "✓ [HOST] PCA setup finished."
  exit 0
fi

###############################################################################
### MODE 2: CONTAINER (worker) ###############################################
###############################################################################

# 0. Ensure required genome_data downloads and local reference builds
echo "==> [CONTAINER] Ensuring required genome_data downloads..."
bash "/app/scripts/helpers/download_data.sh" || exit 1

echo "==> [CONTAINER] Ensuring reference builds (BCF/CSI and bref3)..."
bash "/app/scripts/helpers/build_refs.sh" || exit 1

# 0.5 Generate/ensure ALT alleles lookup database
ALT_DB="/app/genome_data/alt_alleles.db"
REF_BCFS_DIR="/app/genome_data/ref_bcfs_b38"
if [ -s "$ALT_DB" ]; then
  echo "✓ [CONTAINER] ALT alleles DB already exists at $ALT_DB. Skipping."
else
  echo "==> [CONTAINER] Building ALT alleles DB at $ALT_DB..."
  mkdir -p "/app/genome_data"
  # build_alt_db.py writes to ./genome_data by default; ensure cwd=/app and pass explicit ref dir
  (
    cd /app && \
    python3 /app/scripts/helpers/build_alt_db.py --ref-dir "$REF_BCFS_DIR"
  ) || { echo "✗ [CONTAINER] Failed to build ALT alleles DB" >&2; exit 1; }
  if [ -s "$ALT_DB" ]; then
    echo "✓ [CONTAINER] ALT alleles DB ready: $ALT_DB"
  else
    echo "✗ [CONTAINER] ALT alleles DB missing or empty after build." >&2
    exit 1
  fi
fi

# 1. Generate LD-pruned panel 
if [ -f "$PRUNED_PATH" ]; then
  echo "✓ [CONTAINER] Pruned panel already exists at $PRUNED_PATH. Skipping."
else
  echo "==> [CONTAINER] Running pruned panel script..."
  bash "/app/scripts/helpers/pruned_panel.sh"
fi

# 2. Make per-ancestry maf files
MAF_FILES=(
  "/app/pca_model/maf_AFR.grch38.tsv.gz"
  "/app/pca_model/maf_AMR.grch38.tsv.gz"
  "/app/pca_model/maf_EAS.grch38.tsv.gz"
  "/app/pca_model/maf_EUR.grch38.tsv.gz"
  "/app/pca_model/maf_SAS.grch38.tsv.gz"
)

# Check if all files already exist
all_exist=true
for f in "${MAF_FILES[@]}"; do
  if [ ! -f "$f" ]; then
    all_exist=false
    break
  fi
done

if [ "$all_exist" = true ]; then
  echo "✓ [CONTAINER] All per-ancestry MAF files already exist. Skipping generation."
else
  echo "==> [CONTAINER] Running per-ancestry MAF script..."
  bash "/app/scripts/helpers/per_ancestry_maf.sh"
fi


# 3. Run ancestry PCA
echo "==> [CONTAINER] Running PCA setup..."
PCA_FILES=(
  "/app/pca_model/pca_sites.b38.tsv"
  "/app/pca_model/ref_means.npy"
  "/app/pca_model/ref_stds.npy"
  "/app/pca_model/loadings.npy"
  "/app/pca_model/ref_scores.csv"
  "/app/pca_model/classifier.pkl"
  "/app/pca_model/meta.json"
)

all_exist=true
for f in "${PCA_FILES[@]}"; do
  if [ ! -f "$f" ]; then
    all_exist=false
    break
  fi
done

if [ "$all_exist" = true ]; then
  echo "✓ [CONTAINER] All PCA files already exist. Skipping generation."
else
  echo "==> [CONTAINER] Running PCA script..."
  python3 /app/scripts/analyses/fit_pca_1kg.py fit-ref \
    --vcf-pattern "${VCF_PATTERN}" \
    --labels "${LABELS_PATH}" \
    --sites "${SITES_PATH}" \
    --out "${OUT_DIR}" \
    --pcs 12 \
    --exclude-mhc --exclude-long-range


  echo "==> [CONTAINER] Cleaning up non-essential PCA files..."
  find /app/pca_model -mindepth 1 -maxdepth 1 \
    ! -name "pca_sites.b38.tsv" \
    ! -name "ref_means.npy" \
    ! -name "ref_stds.npy" \
    ! -name "loadings.npy" \
    ! -name "ref_scores.csv" \
    ! -name "classifier.pkl" \
    ! -name "meta.json" \
    ! -name "maf_AFR.grch38.tsv.gz" \
    ! -name "maf_AMR.grch38.tsv.gz" \
    ! -name "maf_EAS.grch38.tsv.gz" \
    ! -name "maf_EUR.grch38.tsv.gz" \
    ! -name "maf_SAS.grch38.tsv.gz" \
    ! -name "labels.keep" \
    -exec rm -rf -- {} +
fi


