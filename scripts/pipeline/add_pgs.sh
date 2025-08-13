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

REGISTRY_CSV="${WEIGHTS_HM_DIR}/registry.csv"
STD_MAIN="/app/pgs/weights/standardization/standardization.tsv"

pgs_already_harmonized(){
  local id="$1"
  # Prefer registry check (CSV, first column is pgs_id)
  if [[ -s "$REGISTRY_CSV" ]]; then
    if awk -F',' -v id="$id" 'NR>1 && $1==id {exit 0} END{exit 1}' "$REGISTRY_CSV"; then
      return 0
    fi
  fi
  # Fallback: any harmonized file present
  shopt -s nullglob
  local have=0 f
  for f in "${WEIGHTS_HM_DIR}/${id}."*.b38.tsv*; do
    if [[ -s "$f" ]]; then have=1; break; fi
  done
  shopt -u nullglob
  [[ "$have" -eq 1 ]]
}

pgs_standardized(){
  local id="$1"
  # Check standardization table for rows with this PGS ID (TSV: first column pgs_id)
  if [[ -s "$STD_MAIN" ]] && awk -F'\t' -v id="$id" 'NR>1 && $1==id {exit 0} END{exit 1}' "$STD_MAIN"; then
    return 0
  fi
  return 1
}

MISSING_STD=()
for PGSID in "$@"; do
  echo "  [+] ${PGSID}"
  if pgs_already_harmonized "$PGSID"; then
    echo "     ↳ Found in registry or files; skipping harmonization."
  else
    python3 "${WEIGHTS_HARMONIZER}" \
      --pgs-id "${PGSID}" \
      --maf-dir "${PCA_DIR}" \
      --out-dir "${WEIGHTS_HM_DIR}" \
      --gzip -v
  fi

  if ! pgs_standardized "$PGSID"; then
    MISSING_STD+=("$PGSID")
  fi
done

# Build/refresh standardization only if any requested PGS lacks entries
if (( ${#MISSING_STD[@]} > 0 )); then
  echo "==> [CONTAINER] Building/refreshing standardization (missing for: ${MISSING_STD[*]})"
  python3 "/app/scripts/analyses/build_standardization.py"
  # Standardization written to $STD_MAIN
else
  echo "✓ [CONTAINER] Standardization already present for all requested PGS IDs. Skipping."
fi

echo "✅ Done. Outputs in: ${WEIGHTS_HM_DIR}"
