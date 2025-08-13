#!/usr/bin/env bash
# /app/scripts/analyses/score_all_pgs.sh
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <USER_PFILE_PREFIX> <USER_IID> <SUBPOP> <RESULTS_DIR>" >&2
  exit 1
fi

USER_PFILE="$1"
USER_ID="$2"
SUBPOP="$3"   # EUR/AFR/EAS/AMR/SAS
RESULTS_DIR="$4"

# Use fixed harmonized weights directory inside the container
HM_DIR="/app/pgs/weights/harmonized"
if [[ ! -d "$HM_DIR" ]]; then
  echo "[ERROR] Harmonized weights directory not found: $HM_DIR" >&2
  exit 2
fi

work="$(mktemp -d -p /tmp scorepgs.XXXX)"
trap 'rm -rf "$work"' EXIT
mkdir -p "$RESULTS_DIR"

# Normalize SUBPOP (uppercase, strip CR/spaces)
SUBPOP=$(printf '%s' "$SUBPOP" | tr '\r' '\n' | head -n1 | tr -d '[:space:]')
SUBPOP=${SUBPOP^^}

# 1) Standardize IDs to match weights' 'chr' convention when needed
PVAR_IN="${USER_PFILE}.pvar"
if [[ ! -s "$PVAR_IN" ]]; then
  echo "[ERROR] Missing input PVAR: $PVAR_IN" >&2; exit 2
fi
CHR_FIELD=$(awk 'BEGIN{FS="\t"} !/^#/ {print $1; exit}' "$PVAR_IN" 2>/dev/null || true)
SET_FMT="@:#:\$r:\$a"
if [[ "$CHR_FIELD" != chr* ]]; then
  # Input uses numeric contigs; add 'chr' prefix so IDs match weights (e.g., chr1:pos:ref:alt)
  SET_FMT="chr@:#:\$r:\$a"
fi
plink2 --pfile "$USER_PFILE" \
  --set-all-var-ids "$SET_FMT" --new-id-max-allele-len 200 truncate \
  --make-pgen --out "$work/user_stdids"

shopt -s nullglob
# Collect matching files; tolerate both plain and gz
declare -a WEIGHT_FILES=()
for pat in \
  "$HM_DIR"/*."${SUBPOP}".b38.tsv \
  "$HM_DIR"/*."${SUBPOP}".b38.tsv.gz \
  "$HM_DIR"/*."${SUBPOP}".b38.tsv.gs; do
  for f in $pat; do
    [[ -f "$f" ]] && WEIGHT_FILES+=("$f")
  done
done

if (( ${#WEIGHT_FILES[@]} == 0 )); then
  echo "[INFO] No weights found in $HM_DIR for SUBPOP=$SUBPOP (patterns: *.${SUBPOP}.b38.tsv[.gz])." >&2
  echo "[OK] Nothing to score." >&2
  exit 0
fi
IFS=$'\n' WEIGHT_FILES=($(printf '%s\n' "${WEIGHT_FILES[@]}" | sort)); unset IFS
echo "[INFO] Using weights dir: $HM_DIR (${#WEIGHT_FILES[@]} files for $SUBPOP)" >&2

score_one(){
  local weights_path="$1"
  local base
  base=$(basename -- "$weights_path")
  local pgsid
  pgsid="${base%%.*}"

  local score_tsv="$work/${pgsid}.${SUBPOP}.score.tsv"
  # Build minimal score (ID A1 BETA) from harmonized weights
  if [[ "$weights_path" == *.gz || "$weights_path" == *.gs ]]; then
    gzip -cd -- "$weights_path" | \
      awk 'BEGIN{FS=OFS="\t"; printed=0} NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "ID","A1","BETA"; next}
           {id=$h["chr"] ":" $h["pos"] ":" $h["ref"] ":" $h["alt"];
            if(h["effect_allele"] && h["beta"]) print id, $h["effect_allele"], $h["beta"];}' \
      > "$score_tsv"
  else
    awk 'BEGIN{FS=OFS="\t"; printed=0} NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "ID","A1","BETA"; next}
         {id=$h["chr"] ":" $h["pos"] ":" $h["ref"] ":" $h["alt"];
          if(h["effect_allele"] && h["beta"]) print id, $h["effect_allele"], $h["beta"];}' \
      "$weights_path" > "$score_tsv"
  fi

  if [[ ! -s "$score_tsv" ]]; then
    echo "[WARN] Empty score table for $weights_path; skipping." >&2
    return 0
  fi

  echo "[+] Scoring $pgsid ($SUBPOP)"
  plink2 \
    --pfile "$work/user_stdids" \
    --score "$score_tsv" 1 2 3 header-read no-mean-imputation list-variants \
    --out "${RESULTS_DIR}/${pgsid}.${SUBPOP}"
}

for wf in "${WEIGHT_FILES[@]}"; do
  if [[ ! -f "$wf" ]]; then
    echo "[WARN] Missing weights path: $wf" >&2
    continue
  fi
  score_one "$wf"
done

echo "[OK] Scoring done -> $RESULTS_DIR"
