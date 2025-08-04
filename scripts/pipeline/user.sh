#!/usr/bin/env bash
set -Eeuo pipefail

# run_user.sh
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

  # 1. Argument and Sanity Checks
  if [[ $# -lt 1 ]]; then
    echo "Usage (on Host): $0 </path/to/your/file> [pipeline options...]"
    exit 1
  fi
  IN_FILE_PATH="$1"
  PIPELINE_ARGS=("${@:2}") # Capture all other arguments to pass them to the container

  # 2. Resolve absolute path to the input file and its directory.
  #    This dependency-free method works on both macOS and Linux.
  if [[ ! -f "$IN_FILE_PATH" ]]; then
    echo "Error: Input file not found at '$IN_FILE_PATH'" >&2
    exit 1
  fi

  ORIGINAL_PWD=$(pwd) # Save our starting directory
  INPUT_DIR_RAW=$(dirname "$IN_FILE_PATH")
  FILENAME=$(basename "$IN_FILE_PATH")
  
  cd "$INPUT_DIR_RAW" # Temporarily change to the input file's directory
  ABS_INPUT_DIR=$(pwd) # Get the absolute path of that directory
  cd "$ORIGINAL_PWD" # IMPORTANT: Go back to where we started

  # 3. Build Docker command and launch
  echo "==> [HOST] Mounting '$ABS_INPUT_DIR' to /app/input"
  DOCKER_IMAGE="pgspilot"

  # --- START of the fix ---
  # Define the command and its mandatory arguments to be run inside the container
  CMD_TO_RUN=(/app/scripts/pipeline/user.sh "/app/input/$FILENAME")

  # If the PIPELINE_ARGS array is not empty, append its contents
  if [[ ${#PIPELINE_ARGS[@]} -gt 0 ]]; then
    CMD_TO_RUN+=("${PIPELINE_ARGS[@]}")
  fi
  # --- END of the fix ---

  docker run --rm \
    -e INSIDE_DOCKER=true \
    -v "$ABS_INPUT_DIR:/app/input:ro" \
    -v "$(pwd)/genome_data:/app/genome_data:ro" \
    -v "$(pwd)/pca_model:/app/pca_model:ro" \
    -v "$(pwd)/traits:/app/traits:ro" \
    -v "$(pwd)/weights_hm:/app/weights_hm:ro" \
    -v "$(pwd)/weights_raw:/app/weights_raw:ro" \
    -v "$(pwd)/users:/app/users" \
    -v "$(pwd)/calibration:/app/calibration" \
    -v "$(pwd)/scripts:/app/scripts:ro" \
    "$DOCKER_IMAGE" \
    "${CMD_TO_RUN[@]}" # <-- Use the safely constructed command array here

  echo "✓ [HOST] Pipeline finished."

else
  ##############################################################################
  ### MODE 2: CONTAINER (WORKER)                                             ###
  # This block runs when the script is executed inside the Docker container.   #
  ##############################################################################

  echo "==> [CONTAINER] Script running inside Docker. Starting pipeline..."

  # --- This is your original pipeline logic ---

  # Source config from its fixed path inside the container
  source "/app/scripts/config.sh"

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
  USER_DIR="$USERS_DIR/$STEM"
  mkdir -p "$USER_DIR"

  echo "==> [CONTAINER] Input: $INPUT"
  echo "==> [CONTAINER] STEM:  $STEM"
  echo "==> [CONTAINER] Out:   $USER_DIR"

  # 1) Impute if needed
  FINAL_VCF="$USER_DIR/${STEM}_imputed_all.vcf.gz"
  if [[ "$INPUT" == *".vcf.gz" ]]; then
    echo "==> [CONTAINER] Using provided VCF (skipping imputation)"
    if [[ ! -f "$FINAL_VCF" ]]; then
      # Copy into the standard location (keep user dir self-contained)
      cp -v "$INPUT" "$FINAL_VCF"
      tabix -f -p vcf "$FINAL_VCF" || true
    fi
  else
    echo "==> [CONTAINER] Running imputation…"
    # This calls the imputation script, which is found at the path defined in config.sh
    # e.g., IMPUTE_SH="/app/scripts/impute.sh"
    "$IMPUTE_SH" "$INPUT"

    # The imputation script creates its output in a subdir of /app/users/
    IMPUTE_OUT_DIR="/app/users/${STEM}_results"
    CAND="${IMPUTE_OUT_DIR}/${STEM}_imputed_all.vcf.gz"

    [[ -f "$CAND" ]] || { echo "✗ Could not find imputation output: $CAND" >&2; exit 1; }

    echo "==> [CONTAINER] Moving imputation results to final user directory."
    mv -v "$CAND" "$FINAL_VCF"
    mv -v "$CAND.tbi" "$FINAL_VCF.tbi"
    # Optional: Move the entire imputed/phased dirs for inspection
    mv -v "${IMPUTE_OUT_DIR}/imputed_dir" "${IMPUTE_OUT_DIR}/phased_dir" "$USER_DIR/"
    rmdir "$IMPUTE_OUT_DIR" # Clean up the now-empty results folder
  fi

  # 2) Sanity: required assets
  # for req in "$FINAL_VCF" "$PCA_DIR/loadings.npy" "$TRAITS_DIR" "$WEIGHTS_HM_DIR" "$CALIB_DIR"; do
  #   [[ -e "$req" ]] || { echo "✗ Missing required input/asset: $req" >&2; exit 1; }
  # done

  # # 3) Score + report
  # echo "==> [CONTAINER] Scoring & reporting…"
  # SCORER_ARGS=()
  # [[ -n "$ONLY_TRAITS" ]] && SCORER_ARGS+=(--only-traits "$ONLY_TRAITS")
  # [[ -n "$COVERAGE_THRESH" ]] && SCORER_ARGS+=(--coverage-thresh "$COVERAGE_THRESH")

  # $PYTHON "$SCORER" \
  #   --vcf "$FINAL_VCF" \
  #   --pca "$PCA_DIR" \
  #   --traits_dir "$TRAITS_DIR" \
  #   --weights_dir "$WEIGHTS_HM_DIR" \
  #   --calib_dir "$CALIB_DIR" \
  #   --out_dir "$USER_DIR" \
  #   --report-md \
  #   "${SCORER_ARGS[@]}"

  echo "CONTAINER: Done."
  echo "CONTAINER: Imputation finished successfully."
  echo "Imputed VCF located at: $FINAL_VCF"
  # echo "  Ancestry JSON:    $USER_DIR/ancestry.json"
  # echo "  All scores JSON:  $USER_DIR/scores_all.json"
  # echo "  Per-trait JSONs:  $USER_DIR/scores_<trait>.json"
  # echo "  Reports:          $USER_DIR/report_<trait>.md"
fi