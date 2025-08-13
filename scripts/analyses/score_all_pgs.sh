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

# 1) Standardize IDs; later we'll align score ID chr-prefix to match these IDs
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

# Build a fast lookup from (chr:pos:ref:alt) → dataset VARID from user_stdids.pvar
# We normalize chr to numeric (strip leading 'chr') for robust joins.
PVAR_LOOKUP="$work/pvar_lookup.tsv"
awk 'BEGIN{FS=OFS="\t"} 
     /^#/ {next}
     {
       ch=$1; sub(/^chr/, "", ch);
       ref=toupper($4);
       alts=$5; gsub(/ /, "", alts);
       n=split(alts, arr, ",");
       for(i=1;i<=n;i++){
         alt=toupper(arr[i]);
         key1=ch ":" $2 ":" ref ":" alt;
         key2=ch ":" $2 ":" alt ":" ref;
         print key1, $3;
         print key2, $3;
       }
     }' \
    "$work/user_stdids.pvar" > "$PVAR_LOOKUP"
if [[ ! -s "$PVAR_LOOKUP" ]]; then
  echo "[ERROR] Failed to build PVAR lookup from $work/user_stdids.pvar" >&2
  exit 2
fi

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

# Deduplicate by logical PGS identifier to avoid overwriting outputs when both
# compressed and uncompressed weight files are present for the same PGS and SUBPOP.
# Preference order: .tsv.gz > .tsv > .tsv.gs
if (( ${#WEIGHT_FILES[@]} > 1 )); then
  declare -A _chosen_path
  declare -A _chosen_rank
  for _wf in "${WEIGHT_FILES[@]}"; do
    _base=$(basename -- "$_wf")
    if [[ "$_base" =~ ^(PGS[0-9]+)\. ]]; then
      _pgsid="${BASH_REMATCH[1]}"
    else
      _pgsid="${_base%%.*}"
    fi
    _key="${_pgsid}.${SUBPOP}"
    _rank=0
    if [[ "$_base" == *.tsv.gz ]]; then
      _rank=3
    elif [[ "$_base" == *.tsv ]]; then
      _rank=2
    elif [[ "$_base" == *.tsv.gs ]]; then
      _rank=1
    fi
    _prev_rank=${_chosen_rank["$_key"]:- -1}
    if [[ -z "$_prev_rank" ]]; then _prev_rank=-1; fi
    if (( _rank > _prev_rank )); then
      _chosen_path["$_key"]="$_wf"
      _chosen_rank["$_key"]=$_rank
    fi
  done
  # Rebuild WEIGHT_FILES from chosen unique keys
  WEIGHT_FILES=()
  for _k in "${!_chosen_path[@]}"; do
    WEIGHT_FILES+=( "${_chosen_path[$_k]}" )
  done
  IFS=$'\n' WEIGHT_FILES=($(printf '%s\n' "${WEIGHT_FILES[@]}" | sort)); unset IFS
  echo "[INFO] Selected ${#WEIGHT_FILES[@]} unique PGS IDs for $SUBPOP after de-duplication" >&2
fi

score_one(){
  local weights_path="$1"
  local base
  base=$(basename -- "$weights_path")
  local pgsid
  pgsid="${base%%.*}"

  local score_tsv="$work/${pgsid}.${SUBPOP}.score.tsv"
  local stats_path="$work/${pgsid}.${SUBPOP}.stats"
  # Build minimal score (ID A1 BETA) from harmonized weights
  if [[ "$weights_path" == *.gz || "$weights_path" == *.gs ]]; then
    gzip -cd -- "$weights_path" | awk -v map_path="$PVAR_LOOKUP" -v label="$weights_path" -v stats_path="$stats_path" '
           BEGIN{FS="[\t ]+"; OFS="\t"; miss_printed=0; while ((getline < map_path) > 0) { if($1!="") m[$1]=$2; nmap++ } close(map_path)}
           NR==1{
             # Remove potential UTF-8 BOM on first header token
             sub(/^\xef\xbb\xbf/, "", $1)
             # Build a lowercase header-name -> column-index map (alnum + underscore)
             for(i=1;i<=NF;i++) { k=tolower($i); gsub(/[^a-z0-9_]/, "", k); h[k]=i }
            # Resolve synonyms (prefer harmonized coordinates/alleles when available)
            c_chr = (h["hm_chr"] ? h["hm_chr"] : (h["hm_chrom"] ? h["hm_chrom"] : (h["chr"] ? h["chr"] : (h["chrom"] ? h["chrom"] : (h["chromosome"] ? h["chromosome"] : (h["chr_name"] ? h["chr_name"] : 0))))))
            c_pos = (h["hm_pos"] ? h["hm_pos"] : (h["hm_position"] ? h["hm_position"] : (h["pos_b38"] ? h["pos_b38"] : (h["bp38"] ? h["bp38"] : (h["base_pair_location"] ? h["base_pair_location"] : (h["pos"] ? h["pos"] : (h["position"] ? h["position"] : (h["chr_position"] ? h["chr_position"] : (h["bp"] ? h["bp"] : 0)))))))))
            c_ref = (h["hm_ref"] ? h["hm_ref"] : (h["hm_other_allele"] ? h["hm_other_allele"] : (h["ref"] ? h["ref"] : (h["a0"] ? h["a0"] : (h["other_allele"] ? h["other_allele"] : (h["non_effect_allele"] ? h["non_effect_allele"] : (h["noneffect"] ? h["noneffect"] : 0)))))))
            c_alt = (h["hm_alt"] ? h["hm_alt"] : (h["hm_effect_allele"] ? h["hm_effect_allele"] : (h["alt"] ? h["alt"] : (h["a1"] ? h["a1"] : (h["effect_allele"] ? h["effect_allele"] : (h["effect"] ? h["effect"] : 0))))))
            c_a1  = (h["hm_effect_allele"] ? h["hm_effect_allele"] : (h["effect_allele"] ? h["effect_allele"] : (h["a1"] ? h["a1"] : (h["alt"] ? h["alt"] : 0))))
             c_bta = (h["beta"] ? h["beta"] : (h["effect_weight"] ? h["effect_weight"] : (h["weight"] ? h["weight"] : (h["log_odds"] ? h["log_odds"] : (h["logodds"] ? h["logodds"] : 0)))))
             print "SNP","A1","BETA"; next
           }
           {
             if (!c_a1 || !c_bta) next
             if (!c_chr || !c_pos || !c_ref || !c_alt) next
             chrval = $(c_chr); pos=$(c_pos); ref=$(c_ref); alt=$(c_alt);
             # Normalize chr prefix from weights to numeric string
             sub(/^chr/, "", chrval)
             gsub(/[ \r]/, "", chrval); gsub(/[ \r]/, "", pos); gsub(/[ \r]/, "", ref); gsub(/[ \r]/, "", alt)
             if (chrval=="" || pos=="" || ref=="" || alt=="" || chrval=="." || pos=="." || ref=="." || alt==".") next
             key = chrval ":" pos ":" toupper(ref) ":" toupper(alt)
             id = m[key]
             a1 = $(c_a1); b = $(c_bta)
             gsub(/[ \r]/, "", a1); gsub(/[ \r]/, "", b)
             if (a1=="" || b=="" || a1=="." || toupper(a1)=="NA" || toupper(b)=="NA") next
             tot++
             if (id!="") { matched++; print id, a1, b }
             else {
               if (miss_printed < 5) { miss_printed++; printf("[DEBUG] key_miss: %s\n", key) > "/dev/stderr" }
             }
           }
           END{
            printf("[DEBUG] %s: scanned=%d matched=%d\n", label, tot+0, matched+0) > "/dev/stderr"
            if (stats_path!="") { printf("scanned\t%d\nmatched\t%d\n", tot+0, matched+0) > stats_path }
           }' > "$score_tsv"
  else
    awk -v map_path="$PVAR_LOOKUP" -v label="$weights_path" -v stats_path="$stats_path" '
         BEGIN{FS="[\t ]+"; OFS="\t"; miss_printed=0; while ((getline < map_path) > 0) { if($1!="") m[$1]=$2; nmap++ } close(map_path)}
         NR==1{
           sub(/^\xef\xbb\xbf/, "", $1)
          for(i=1;i<=NF;i++) { k=tolower($i); gsub(/[^a-z0-9_]/, "", k); h[k]=i }
          c_chr = (h["hm_chr"] ? h["hm_chr"] : (h["hm_chrom"] ? h["hm_chrom"] : (h["chr"] ? h["chr"] : (h["chrom"] ? h["chrom"] : (h["chromosome"] ? h["chromosome"] : (h["chr_name"] ? h["chr_name"] : 0))))))
          c_pos = (h["hm_pos"] ? h["hm_pos"] : (h["hm_position"] ? h["hm_position"] : (h["pos_b38"] ? h["pos_b38"] : (h["bp38"] ? h["bp38"] : (h["base_pair_location"] ? h["base_pair_location"] : (h["pos"] ? h["pos"] : (h["position"] ? h["position"] : (h["chr_position"] ? h["chr_position"] : (h["bp"] ? h["bp"] : 0)))))))))
          c_ref = (h["hm_ref"] ? h["hm_ref"] : (h["hm_other_allele"] ? h["hm_other_allele"] : (h["ref"] ? h["ref"] : (h["a0"] ? h["a0"] : (h["other_allele"] ? h["other_allele"] : (h["non_effect_allele"] ? h["non_effect_allele"] : (h["noneffect"] ? h["noneffect"] : 0)))))))
          c_alt = (h["hm_alt"] ? h["hm_alt"] : (h["hm_effect_allele"] ? h["hm_effect_allele"] : (h["alt"] ? h["alt"] : (h["a1"] ? h["a1"] : (h["effect_allele"] ? h["effect_allele"] : (h["effect"] ? h["effect"] : 0))))))
           c_a1  = (h["hm_effect_allele"] ? h["hm_effect_allele"] : (h["effect_allele"] ? h["effect_allele"] : (h["a1"] ? h["a1"] : (h["alt"] ? h["alt"] : 0))))
           c_bta = (h["beta"] ? h["beta"] : (h["effect_weight"] ? h["effect_weight"] : (h["weight"] ? h["weight"] : (h["log_odds"] ? h["log_odds"] : (h["logodds"] ? h["logodds"] : 0)))))
           print "SNP","A1","BETA"; next
         }
         {
           if (!c_a1 || !c_bta) next
           if (!c_chr || !c_pos || !c_ref || !c_alt) next
           chrval = $(c_chr); pos=$(c_pos); ref=$(c_ref); alt=$(c_alt);
           sub(/^chr/, "", chrval)
           gsub(/[ \r]/, "", chrval); gsub(/[ \r]/, "", pos); gsub(/[ \r]/, "", ref); gsub(/[ \r]/, "", alt)
           if (chrval=="" || pos=="" || ref=="" || alt=="" || chrval=="." || pos=="." || ref=="." || alt==".") next
           key = chrval ":" pos ":" toupper(ref) ":" toupper(alt)
           id = m[key]
           a1 = $(c_a1); b = $(c_bta)
           gsub(/[ \r]/, "", a1); gsub(/[ \r]/, "", b)
           if (a1=="" || b=="" || a1=="." || toupper(a1)=="NA" || toupper(b)=="NA") next
           tot++
           if (id!="") { matched++; print id, a1, b }
           else {
             if (miss_printed < 5) { miss_printed++; printf("[DEBUG] key_miss: %s\n", key) > "/dev/stderr" }
           }
         }
         END{
          printf("[DEBUG] %s: scanned=%d matched=%d\n", label, tot+0, matched+0) > "/dev/stderr"
          if (stats_path!="") { printf("scanned\t%d\nmatched\t%d\n", tot+0, matched+0) > stats_path }
         }' "$weights_path" > "$score_tsv"
  fi

  # Debug + emptiness handling based on data rows (exclude header)
  local n_rows
  n_rows=$(awk 'NR>1{c++} END{print c+0}' "$score_tsv" 2>/dev/null || echo 0)
  local scanned=0
  local matched=0
  if [[ -s "$stats_path" ]]; then
    scanned=$(awk -F'\t' '$1=="scanned"{print $2}' "$stats_path" 2>/dev/null || echo 0)
    matched=$(awk -F'\t' '$1=="matched"{print $2}' "$stats_path" 2>/dev/null || echo 0)
  fi
  awk -F'\t' 'BEGIN{empty=0; total=0}
       NR>1{total++; if($1=="" || $1=="." || toupper($1)=="NA") empty++}
       END{printf("[DEBUG] %s: rows=%d, empty_id=%d\n", ARGV[1], total, empty) > "/dev/stderr"}' "$score_tsv" || true
  awk 'NR==1{print "[DEBUG] header:", $0 > "/dev/stderr"; next} NR<=6{print "[DEBUG] row:", $0 > "/dev/stderr"}' OFS='\t' "$score_tsv" || true
  if [[ "${n_rows}" -eq 0 ]]; then
    echo "[WARN] Empty score table for $weights_path; skipping." >&2
    # Emit a clear skip marker so downstream steps and users can see it was attempted
    {
      echo "pgs_id	${pgsid}"
      echo "subpop	${SUBPOP}"
      echo "reason	no_overlap"
      echo "weights_path	${weights_path}"
      echo "variants_in_weights	${scanned}"
      echo "variants_matched	${matched}"
    } > "${RESULTS_DIR}/${pgsid}.${SUBPOP}.skipped.tsv"

    # Also emit empty-shaped PLINK outputs so downstream can consume uniformly
    local out_pref="${RESULTS_DIR}/${pgsid}.${SUBPOP}"
    {
      printf "FID\tIID\tSCORE1_SUM\tSCORE1_AVG\n"
      printf "%s\t%s\t0\t0\n" "${USER_ID}" "${USER_ID}"
    } > "${out_pref}.sscore"
    printf "ID\n" > "${out_pref}.sscore.vars"
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
