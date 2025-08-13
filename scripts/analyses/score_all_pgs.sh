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

# Standardization table (PGS x SUBPOP -> mean, sd) for z-scoring
# Priority: $STANDARDIZATION_TSV -> /app/... -> ./pgs/...
STANDARDIZATION_TSV_DEFAULT="/app/pgs/weights/standardization/standardization.tsv"
STANDARDIZATION_TSV_ALT="/pgs/weights/standardization/standardization.tsv"
STANDARDIZATION_TSV_LOCAL="pgs/weights/standardization/standardization.tsv"
STANDARDIZATION_TSV="${STANDARDIZATION_TSV:-}"
if [[ -z "$STANDARDIZATION_TSV" ]]; then
  if [[ -s "$STANDARDIZATION_TSV_DEFAULT" ]]; then
    STANDARDIZATION_TSV="$STANDARDIZATION_TSV_DEFAULT"
  elif [[ -s "$STANDARDIZATION_TSV_ALT" ]]; then
    STANDARDIZATION_TSV="$STANDARDIZATION_TSV_ALT"
  elif [[ -s "$STANDARDIZATION_TSV_LOCAL" ]]; then
    STANDARDIZATION_TSV="$STANDARDIZATION_TSV_LOCAL"
  else
    STANDARDIZATION_TSV=""
  fi
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

  # Ensure per-PGS output directory and prefix
  local out_dir="${RESULTS_DIR}/${pgsid}"
  mkdir -p "$out_dir"
  local out_pref="${out_dir}/${pgsid}.${SUBPOP}"

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
    } > "${out_pref}.skipped.tsv"

    # Also emit empty-shaped PLINK outputs so downstream can consume uniformly
    {
      printf "FID\tIID\tSCORE1_SUM\tSCORE1_AVG\n"
      printf "%s\t%s\t0\t0\n" "${USER_ID}" "${USER_ID}"
    } > "${out_pref}.sscore"
    printf "ID\n" > "${out_pref}.sscore.vars"

    # Emit normalized summary as well, using standardization mean/sd when available
    local mean_val="" sd_val="" trait_val="" zscore="" err=""
    if [[ -n "$STANDARDIZATION_TSV" && -s "$STANDARDIZATION_TSV" ]]; then
      local _std_line
      _std_line=$(awk -v pgs="$pgsid" -v sp="$SUBPOP" 'BEGIN{FS="\t"; OFS="\t"}
        NR==1{
          for(i=1;i<=NF;i++){k=$i; gsub(/^[[:space:]]+|[[:space:]]+$/, "", k); k=tolower(k); m[k]=i}
          c_pgs=(m["pgs_id"]?m["pgs_id"]:(m["pgs"]?m["pgs"]:0));
          c_sp=(m["subpop"]?m["subpop"]:(m["population"]?m["population"]:0));
          c_mean=(m["mean"]?m["mean"]:0);
          c_sd=(m["sd"]?m["sd"]:0);
          c_trait=(m["trait"]?m["trait"]:0);
          next
        }
        c_pgs && c_sp && c_mean && c_sd {
          if($c_pgs==pgs && toupper($c_sp)==sp){
            printf("%s\t%s\t%s\n", $c_mean, $c_sd, (c_trait?$c_trait:"") );
            exit
          }
        }
      ' "$STANDARDIZATION_TSV" 2>/dev/null || true)
      mean_val=$(printf '%s' "$_std_line" | cut -f1)
      sd_val=$(printf '%s' "$_std_line" | cut -f2)
      trait_val=$(printf '%s' "$_std_line" | cut -f3)
    fi
    if [[ -n "$mean_val" && -n "$sd_val" ]] && awk -v sd="$sd_val" 'BEGIN{exit !(sd+0>0)}'; then
      zscore=$(awk -v r=0 -v mu="$mean_val" -v sd="$sd_val" 'BEGIN{printf("%.10g", (r-mu)/sd)}')
    else
      err="std_row_missing"
      if [[ -z "$STANDARDIZATION_TSV" || ! -s "$STANDARDIZATION_TSV" ]]; then err="no_standardization_table"; fi
    fi
    {
      printf "pgs_id\tsubpop\tiid\traw_score\tzscore\tmean\tsd\tmatched_variant_count\tweights_path\ttrait\terror\n"
      printf "%s\t%s\t%s\t0\t%s\t%s\t%s\t0\t%s\t%s\t%s\n" \
        "$pgsid" "$SUBPOP" "$USER_ID" "${zscore:-}" "${mean_val:-}" "${sd_val:-}" \
        "$weights_path" "${trait_val:-}" "${err:-}"
    } > "${out_pref}.normalized.tsv"
    return 0
  fi

  echo "[+] Scoring $pgsid ($SUBPOP)"
  plink2 \
    --pfile "$work/user_stdids" \
    --score "$score_tsv" 1 2 3 header-read no-mean-imputation list-variants \
    --out "${out_pref}"

  # 2) Derive user-level normalized score (z-score) using standardization.tsv if available
  #    Output one-row TSV: <pgs_id, subpop, iid, raw_score, zscore, mean, sd, matched_variant_count, weights_path, trait>
  local sscore_path="${out_pref}.sscore"
  local svars_path="${out_pref}.sscore.vars"
  local norm_out="${out_pref}.normalized.tsv"

  # Count matched variants (exclude header)
  local matched_variants=0
  if [[ -s "$svars_path" ]]; then
    matched_variants=$(awk 'NR>1{c++} END{print c+0}' "$svars_path" 2>/dev/null || echo 0)
  fi

  # Extract user's raw score from .sscore (prefer SCORE1_SUM, else SCORE1_AVG, else SCORE)
  local raw_score=""
  if [[ -s "$sscore_path" ]]; then
    raw_score=$(awk -v uid="$USER_ID" '
      BEGIN{FS=OFS="\t"}
      NR==1{
        for(i=1;i<=NF;i++){h[$i]=i}
        sc=-1; if("SCORE1_SUM" in h) sc=h["SCORE1_SUM"]; else if("SCORE1_AVG" in h) sc=h["SCORE1_AVG"]; else if("SCORE" in h) sc=h["SCORE"]; next
      }
      sc>0 && $2==uid { print $sc; exit }
    ' "$sscore_path" 2>/dev/null || true)
    # If empty, try whitespace-sep read
    if [[ -z "$raw_score" ]]; then
      raw_score=$(awk -v uid="$USER_ID" '
        NR==1{for(i=1;i<=NF;i++){h[$i]=i}; sc=-1; if("SCORE1_SUM" in h) sc=h["SCORE1_SUM"]; else if("SCORE1_AVG" in h) sc=h["SCORE1_AVG"]; else if("SCORE" in h) sc=h["SCORE"]; next}
        sc>0 && $2==uid { print $sc; exit }
      ' "$sscore_path" 2>/dev/null || true)
    fi
  fi

  # Lookup mean/sd/trait from standardization table
  local mean_val="" sd_val="" trait_val=""
  if [[ -n "$STANDARDIZATION_TSV" && -s "$STANDARDIZATION_TSV" ]]; then
    # Capture into a single line to avoid set -e failures on empty reads
    local _std_line
    _std_line=$(awk -v pgs="$pgsid" -v sp="$SUBPOP" 'BEGIN{FS="\t"; OFS="\t"}
      NR==1{
        for(i=1;i<=NF;i++){k=$i; gsub(/^[[:space:]]+|[[:space:]]+$/, "", k); k=tolower(k); m[k]=i}
        c_pgs=(m["pgs_id"]?m["pgs_id"]:(m["pgs"]?m["pgs"]:0));
        c_sp=(m["subpop"]?m["subpop"]:(m["population"]?m["population"]:0));
        c_mean=(m["mean"]?m["mean"]:0);
        c_sd=(m["sd"]?m["sd"]:0);
        c_trait=(m["trait"]?m["trait"]:0);
        next
      }
      c_pgs && c_sp && c_mean && c_sd {
        if($c_pgs==pgs && toupper($c_sp)==sp){
          printf("%s\t%s\t%s\n", $c_mean, $c_sd, (c_trait?$c_trait:"") );
          exit
        }
      }
    ' "$STANDARDIZATION_TSV" 2>/dev/null || true)
    mean_val=$(printf '%s' "$_std_line" | cut -f1)
    sd_val=$(printf '%s' "$_std_line" | cut -f2)
    trait_val=$(printf '%s' "$_std_line" | cut -f3)
  fi

  # Compute z-score if possible
  local zscore=""; local err=""
  if [[ -n "$raw_score" && -n "$mean_val" && -n "$sd_val" ]]; then
    # Validate numeric sd
    if awk -v sd="$sd_val" 'BEGIN{exit !(sd+0>0)}'; then
      zscore=$(awk -v r="$raw_score" -v mu="$mean_val" -v sd="$sd_val" 'BEGIN{printf("%.10g", (r-mu)/sd)}')
    else
      err="non_positive_sd"
    fi
  else
    if [[ -z "$STANDARDIZATION_TSV" || ! -s "$STANDARDIZATION_TSV" ]]; then
      err="no_standardization_table"
    elif [[ -z "$raw_score" ]]; then
      err="no_user_score"
    else
      err="std_row_missing"
    fi
  fi

  {
    printf "pgs_id\tsubpop\tiid\traw_score\tzscore\tmean\tsd\tmatched_variant_count\tweights_path\ttrait\terror\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n" \
      "$pgsid" "$SUBPOP" "$USER_ID" "${raw_score:-}" "${zscore:-}" "${mean_val:-}" "${sd_val:-}" \
      "$matched_variants" "$weights_path" "${trait_val:-}" "${err:-}"
  } > "$norm_out"
}

idx=0
total=${#WEIGHT_FILES[@]}
for wf in "${WEIGHT_FILES[@]}"; do
  idx=$((idx+1))
  echo "[INFO] (${idx}/${total}) Processing weights: ${wf}" >&2
  if [[ ! -f "$wf" ]]; then
    echo "[WARN] Missing weights path: $wf" >&2
    continue
  fi
  score_one "$wf"
done

echo "[OK] Scoring done -> $RESULTS_DIR"
