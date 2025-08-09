#!/usr/bin/env bash
set -Eeuo pipefail

# user.sh
# Dual-personality script:
#  - On HOST: Acts as a Docker WRAPPER to launch the container.
#  - In CONTAINER: Acts as the WORKER to run the pipeline.

# The switch is the $INSIDE_DOCKER environment variable.

# (Put near top of script)
set -Eeuo pipefail
IFS=$'\n\t'
log(){ printf '%s\n' "$*" >&2; }
die(){ log "✗ $*"; exit 1; }

cleanup_imputation_tmp() {
  shopt -s nullglob
  rm -rf "$USER_DIR/imputed_dir" \
         "$USER_DIR/phased_dir" \
         "$USER_DIR"/isec_chr*/
  rm -f "$USER_DIR"/*.build37.chr.vcf.gz* \
        "$USER_DIR"/*.build37.chr.alt.vcf.gz* \
        "$USER_DIR"/*.lift38.vcf \
        "$USER_DIR"/*.lift38.*.vcf.gz* \
        "$USER_DIR"/*_stranded_clean.vcf.gz* \
        "$USER_DIR"/*.reheadered.vcf.gz* \
        "$USER_DIR"/*.alt.vcf.gz* \
        "$USER_DIR/$STEM.txt.gz" \
        2>/dev/null || true
  shopt -u nullglob
}


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

  # 1) Impute if needed (but reuse if already exists)
  FINAL_VCF="$USER_DIR/${STEM}_imputed_all.vcf.gz"

  if [[ -f "$FINAL_VCF" ]]; then
    echo "==> [CONTAINER] Found existing FINAL_VCF: $FINAL_VCF (skipping imputation)"
    # ensure it's indexed
    [[ -f "${FINAL_VCF}.tbi" ]] || tabix -f -p vcf "$FINAL_VCF" || true
  # Replace your imputation block with this
  else
    echo "==> [CONTAINER] Running imputation…"

    # Run the flaky step without aborting the whole script, but record the status.
    set +e
    OUT_DIR="$USER_DIR" "$IMPUTE_SH" "$INPUT"
    imp_status=$?
    set -e

    # Cleanup should always run, success or failure.
    cleanup_imputation_tmp || true

    # Validate success by artifacts, not by exit code.
    [[ -s "$FINAL_VCF" ]] || die "Imputation exited $imp_status and produced no $FINAL_VCF"

    # Basic integrity checks
    tabix -f -p vcf "$FINAL_VCF" || die "Failed to index $FINAL_VCF"
    nvar=$( (zgrep -vc '^#' "$FINAL_VCF" || echo 0) )
    (( nvar > 0 )) || die "Imputed VCF has zero variants"

    # Optional: only treat certain non-zero codes as warnings
    if [[ "$imp_status" -ne 0 ]]; then
      log "⚠ Imputation returned $imp_status, but outputs validated (variants: $nvar); continuing."
    fi
  fi



  # 2) Preparing PLINK2 PFILE from imputed VCF…
  VCF_TO_PFILE_SH="/app/scripts/pipeline/vcf_to_pfile.sh"

  INFO_KEY="${IMP_INFO_KEY:-}"           # e.g., R2, INFO, or IMPINFO; leave empty to auto-detect
  INFO_MIN="${IMP_INFO_MIN:-0.8}"        # threshold for the chosen INFO key
  KEEP_UNFILTERED="${KEEP_FILTERED:-0}"  # 1 to skip bcftools filtering

  USER_PFILE="$USER_DIR/pfiles/user"
  PANEL_ID_FILE="/app/pca_model/panel.ids"

  have_pfiles() {
    local pref="$1"
    # All three exist and are non-empty?
    [[ -s "${pref}.pgen" && -s "${pref}.pvar" && -s "${pref}.psam" ]] || return 1
    # .pvar has at least 1 variant row
    local nvar
    nvar=$(grep -vc '^#' "${pref}.pvar" || echo 0)
    [[ "$nvar" -gt 0 ]] || return 1
    # If we know the source VCF, ensure pfiles are not older than it
    if [[ -n "${FINAL_VCF:-}" && -f "$FINAL_VCF" ]]; then
      [[ "${pref}.pgen" -nt "$FINAL_VCF" && "${pref}.pvar" -nt "$FINAL_VCF" && "${pref}.psam" -nt "$FINAL_VCF" ]] || return 1
    fi
    return 0
  }

  if have_pfiles "$USER_PFILE"; then
    echo "==> [CONTAINER] Found existing PFILE trio at ${USER_PFILE}.* (reusing)"
  else
    echo "==> [CONTAINER] Building PFILE trio from $FINAL_VCF"
    mkdir -p "$(dirname "$USER_PFILE")"

    args=( --vcf "$FINAL_VCF" --out-prefix "$USER_PFILE" )
    [[ -n "${INFO_KEY:-}" ]] && args+=( --info-key "$INFO_KEY" )
    [[ -n "${INFO_MIN:-}" ]] && args+=( --info-min "$INFO_MIN" )

    if [[ "$KEEP_UNFILTERED" -eq 1 ]]; then
      args+=( --extract-snps "$PANEL_ID_FILE" --keep-unfiltered )
    fi

    "$VCF_TO_PFILE_SH" "${args[@]}" >/dev/null

  fi

  # 3) Project user PCs and assign ancestry
  echo "==> Calling ancestry…"
  USER_PCS="$USER_DIR/user_pcs.tsv"
  python /app/scripts/analyses/assign_ancestry.py \
    --user-pfile "$USER_PFILE" \
    --out-ancestry "$USER_DIR/ancestry.tsv" \
    --out-pcs "$USER_PCS" \
    --pcs 6

  SUBPOP=$(awk -F'\t' 'NR==2{print $2}' "$USER_DIR/ancestry.tsv")
  if [[ -z "$SUBPOP" || "$SUBPOP" == "Uncertain" ]]; then
    echo "✗ Ancestry uncertain; skipping scoring. (Consider fallback policy.)" >&2
    exit 2
  fi
  echo "==> Ancestry: $SUBPOP"


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