#!/usr/bin/env bash
set -Eeuo pipefail

# user.sh
# Dual-personality script:
#  - On HOST: Acts as a Docker WRAPPER to launch the container.
#  - In CONTAINER: Acts as the WORKER to run the pipeline.

# The switch is the $INSIDE_DOCKER environment variable.
if [ -z "${INSIDE_DOCKER:-}" ]; then
  ##############################################################################
  ### MODE 1: HOST (WRAPPER)                                                 ###
  # This block runs when the script is executed on your local machine.       #
  ##############################################################################

  echo "==> [HOST] Detected script running on host. Preparing Docker container..."

  THIS_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
  PROJECT_ROOT=$(cd -- "${THIS_SCRIPT_DIR}/../../" &> /dev/null && pwd)
  # --- END OF THE FIX ---

  # 2. Argument and Sanity Checks
  if [[ $# -lt 1 ]]; then
    echo "Usage (on Host): $0 </path/to/your/file> [pipeline options...]"
    exit 1
  fi
  IN_FILE_PATH="$1"
  PIPELINE_ARGS=("${@:2}")

  # 3. Resolve absolute path to the input file and its directory.
  if [[ ! -f "$IN_FILE_PATH" ]]; then
    echo "Error: Input file not found at '$IN_FILE_PATH'" >&2
    exit 1
  fi
  INPUT_DIR_RAW=$(dirname "$IN_FILE_PATH")
  FILENAME=$(basename "$IN_FILE_PATH")
  ABS_INPUT_DIR=$(cd -- "$INPUT_DIR_RAW" &> /dev/null && pwd)

  # 4. Build Docker command and launch
  echo "==> [HOST] Mounting '$ABS_INPUT_DIR' to /app/input_data"
  DOCKER_IMAGE="pgspilot"

  CMD_TO_RUN=(/app/scripts/pipeline/user.sh "/app/input_data/$FILENAME")
  if [[ ${#PIPELINE_ARGS[@]} -gt 0 ]]; then
    CMD_TO_RUN+=("${PIPELINE_ARGS[@]}")
  fi

  # Now, use the reliable $PROJECT_ROOT variable for all mounts
  docker run --rm \
    -e INSIDE_DOCKER=1 \
    -v "$ABS_INPUT_DIR:/app/input_data" \
    -v "${PROJECT_ROOT}/genome_data:/app/genome_data" \
    -v "${PROJECT_ROOT}/pca_model:/app/pca_model" \
    -v "${PROJECT_ROOT}/traits:/app/traits" \
    -v "${PROJECT_ROOT}/weights_hm:/app/weights_hm" \
    -v "${PROJECT_ROOT}/weights_raw:/app/weights_raw" \
    -v "${PROJECT_ROOT}/users:/app/users" \
    -v "${PROJECT_ROOT}/calibration:/app/calibration" \
    -v "${PROJECT_ROOT}/scripts:/app/scripts:ro" \
    "$DOCKER_IMAGE" \
    "${CMD_TO_RUN[@]}"

  echo "✓ [HOST] Pipeline finished."

else
  ##############################################################################
  ### MODE 2: CONTAINER (WORKER)                                             ###
  # This block runs when the script is executed inside the Docker container.   #
  ##############################################################################

  echo "==> [CONTAINER] Script running inside Docker. Starting pipeline..."


  # --- This is your original pipeline logic ---
  source "/app/scripts/pipeline/config.sh"

  if [[ $# -lt 1 ]]; then
    echo "Usage (in Container): $0 <23andme.txt[.gz] | imputed.vcf.gz> [--traits \"Height,LDL\"] [--coverage 0.95]" >&2
    exit 1
  fi

  INPUT="$1"; shift || true
  ONLY_TRAITS=""
  COVERAGE_THRESH="${COVERAGE_THRESH:-0.95}"
  # parse optional flags
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --traits)
        ONLY_TRAITS="$2"; shift 2;;
      --coverage)
        COVERAGE_THRESH="$2"; shift 2;;
      *)
        echo "Unknown option: $1" >&2; exit 1;;
    esac
  done

  # Figure out STEM and prepare user folder
  BASE=$(basename "$INPUT")
  if [[ "$INPUT" == *".vcf.gz" ]]; then
    STEM=${BASE%%.vcf.gz}
    STEM=${STEM%%_imputed_all}        # normalize
  elif [[ "$INPUT" == *".txt" ]] || [[ "$INPUT" == *".txt.gz" ]]; then
    STEM=${BASE%%.*}
  else
    echo "Input must be a .txt(.gz) 23andMe file or an imputed .vcf.gz" >&2
    exit 1
  fi

  # All paths from config.sh are now correct inside the container (e.g., /app/users)
  USER_DIR="/app/users/$STEM"; mkdir -p "$USER_DIR"
  echo "==> Input: $INPUT"
  echo "==> STEM:  $STEM"
  echo "==> Out:   $USER_DIR"

  # 1) Impute if needed
  FINAL_VCF="$USER_DIR/${STEM}_imputed_all.vcf.gz"

  if [[ $INPUT == *.vcf.gz ]]; then
    echo "==> [CONTAINER] Using provided VCF (skipping imputation)"
    cp -v "$INPUT" "$FINAL_VCF"
    tabix -f -p vcf "$FINAL_VCF" || true
  else
    echo "==> [CONTAINER] Running imputation…"
    OUT_DIR="$USER_DIR" "$IMPUTE_SH" "$INPUT"

    # sanity: ensure final file exists and is indexed
    [[ -f "$FINAL_VCF" ]] || { echo "✗ Could not find expected output: $FINAL_VCF" >&2; exit 1; }
    [[ -f "${FINAL_VCF}.tbi" ]] || tabix -f -p vcf "$FINAL_VCF" || true
  fi



  # # 2) Project user PCs and assign ancestry
  # echo "==> Calling ancestry…"
  # USER_PCS="$USER_DIR/user_pcs.tsv"
  # python /app/scripts/analyses/call_ancestry.py \
  #   --user-pfile "$USER_PFILE" \
  #   --out-ancestry "$USER_DIR/ancestry.tsv" \
  #   --out-pcs "$USER_PCS" \
  #   --pcs 6

  # SUBPOP=$(awk -F'\t' 'NR==2{print $2}' "$USER_DIR/ancestry.tsv")
  # if [[ -z "$SUBPOP" || "$SUBPOP" == "Uncertain" ]]; then
  #   echo "✗ Ancestry uncertain; skipping scoring. (Consider fallback policy.)" >&2
  #   exit 2
  # fi
  # echo "==> Ancestry: $SUBPOP"


  # # 3) Score all PGS for that ancestry
  # SCORES_DIR="$USER_DIR/pgs_scores"; mkdir -p "$SCORES_DIR"
  # bash /app/scripts/analyses/score_all_pgs.sh "$USER_PFILE" "$STEM" "$SUBPOP" "$SCORES_DIR"

  # # 4) Collect + standardize → JSON
  # STD_TABLE="/app/weights_hm/pgs_standardization.tsv"
  # OUT_JSON="$USER_DIR/${STEM}_${SUBPOP}_pgs.json"
  # python /app/scripts/analyses/collect_scores.py \
  #   --results-dir "$SCORES_DIR" \
  #   --standardization "$STD_TABLE" \
  #   --user-iid "$STEM" \
  #   --subpop "$SUBPOP" \
  #   --out-json "$OUT_JSON"

  # echo "✓ Done. Results: $OUT_JSON"

fi