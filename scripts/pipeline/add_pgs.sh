#!/usr/bin/env bash
set -Eeuo pipefail

# scripts/pipeline/add_pgs.sh
# Dual-personality script (like user.sh):
#  - On HOST: wrapper that launches Docker with project mounted at /app
#  - In CONTAINER: worker that runs the harmonizer as before

# Host vs container switch
if [ -z "${INSIDE_DOCKER:-}" ]; then
  ##############################################################################
  ### MODE 1: HOST (WRAPPER)                                                 ###
  ##############################################################################

  if [[ $# -lt 1 ]]; then
    echo "Usage (on Host): $0 PGS_ID [PGS_ID ...]" >&2
    exit 1
  fi

  # Locate project root (repo root = two levels up from this script)
  THIS_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
  PROJECT_ROOT=$(cd -- "${THIS_SCRIPT_DIR}/../../" &> /dev/null && pwd)

  DOCKER_IMAGE="${DOCKER_IMAGE:-pgspilot}"

  # Re-invoke this script inside the container, forwarding all args
  docker run --rm \
    -e INSIDE_DOCKER=1 \
    -v "${PROJECT_ROOT}:/app" \
    "$DOCKER_IMAGE" \
    /app/scripts/pipeline/add_pgs.sh "$@"

  exit 0
fi

##############################################################################
### MODE 2: CONTAINER (WORKER)                                             ###
##############################################################################

# Resolve this script's directory and source config in-container
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
CONFIG="${SCRIPT_DIR}/config.sh"
if [[ -f "$CONFIG" ]]; then
  # shellcheck source=/dev/null
  source "$CONFIG"
else
  echo "WARN: No config.sh at $CONFIG — using built-in defaults." >&2
fi

if [[ $# -lt 1 ]]; then
  echo "Usage (in Container): $0 PGS_ID [PGS_ID ...]" >&2
  exit 1
fi

echo "==> [CONTAINER] Fetching harmonized weights + subpop MAF filtering..."
for PGSID in "$@"; do
  echo "  [+] ${PGSID}"
  python3 "${WEIGHTS_HARMONIZER}" \
    --pgs-id "${PGSID}" \
    --maf-dir "${PCA_DIR}" \
    --out-dir "${WEIGHTS_HM_DIR}" \
    --gzip -v
done

echo "✅ Done. Outputs in: ${WEIGHTS_HM_DIR}"



# echo "==> [2/3] Building frozen include lists…"
# SITE_MASK_ARGS=()
# [[ -n "$SITE_MASK" ]] && SITE_MASK_ARGS=(--site-mask "$SITE_MASK" --mask-min-dr2 "$MASK_MIN_DR2")
# $PYTHON "$INCLUDE_BUILDER" \
#   --weights "$WEIGHTS_HM_DIR" \
#   --scorable "$SCORABLE_SITES" \
#   --include-chroms "$INCLUDE_CHROMS" \
#   $DROP_MHC \
#   "${SITE_MASK_ARGS[@]}" \
#   --out "$TRAITS_DIR"

# echo "==> [3/3] Building calibration from 1000G…"
# $PYTHON "$CALIB_BUILDER" \
#   --weights_dir "$WEIGHTS_HM_DIR" \
#   --traits_dir "$TRAITS_DIR" \
#   --vcf_pattern "$ONEKG_VCF_PATTERN" \
#   --labels "$ONEKG_LABELS" \
#   --out "$CALIB_DIR" \
#   --quantiles "$QUANTILES" \
#   --min-variants "$MIN_VARIANTS"

# echo "✓ Prep complete."
# echo "   Harmonized weights: $WEIGHTS_HM_DIR"
# echo "   Include lists:      $TRAITS_DIR"
# echo "   Calibration:        $CALIB_DIR"
