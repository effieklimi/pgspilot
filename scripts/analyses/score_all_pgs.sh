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

REGISTRY="/app/weights_hm/pgs_harmonization_registry.tsv"

work="$(mktemp -d -p /tmp scorepgs.XXXX)"
trap 'rm -rf "$work"' EXIT
mkdir -p "$RESULTS_DIR"

# 1) Standardize IDs
plink2 --pfile "$USER_PFILE" \
  --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 200 truncate \
  --make-pgen --out "$work/user_stdids"

# 2) For each matching (PGS, SUBPOP) in registry, score
# registry columns we wrote: pgs_id (col2), ancestry (col7), out_path (col12)
tail -n +2 "$REGISTRY" | awk -v sp="$SUBPOP" -F'\t' '$7==sp {print $2"\t"$12}' \
| while IFS=$'\t' read -r PGSID OUTPATH; do
  [[ -f "$OUTPATH" ]] || { echo "[WARN] Missing weights $OUTPATH"; continue; }

  awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
       {id=$h["chr"] ":" $h["pos"] ":" $h["ref"] ":" $h["alt"];
        print id, $h["effect_allele"], $h["beta"]}' \
      "$OUTPATH" > "$work/${PGSID}.${SUBPOP}.score.tsv"

  plink2 \
    --pfile "$work/user_stdids" \
    --score "$work/${PGSID}.${SUBPOP}.score.tsv" 1 2 3 header-read no-mean-imputation list-variants \
    --out "${RESULTS_DIR}/${PGSID}.${SUBPOP}"
done

echo "[OK] Scoring done -> $RESULTS_DIR"
