#!/usr/bin/env bash
###############################################################################
# impute.sh – 23andMe text → GRCh38 phased & imputed VCF (chr1–22)
#
# USAGE:
#   ./impute_23andme.sh  <path/to/23andme_genome>.txt[.gz]
#
# OUTPUTS (per sample STEM from input filename):
#   ${STEM}_results/
#     ├─ phased_dir/            # Eagle output per chr
#     ├─ imputed_dir/           # Beagle output per chr (with DS & GP)
#     ├─ imputation_qc_reports/            # QC analysis reports
#     └─ ${STEM}_imputed_all.vcf.gz (+ .tbi)
#
# NOTES:
#   • Input build auto-detected (b37|b38); b36/hg18 not supported here.
#   • Uses 1000G high-coverage GRCh38 refs for Eagle and Beagle.
#   • DS (dosage of ALT) is included for scoring; GP also kept (fallback).
#   • Ancestry-specific PRS standardization is NOT done here.
###############################################################################
set -Eeuo pipefail
[[ "${DEBUG:-}" == 1 ]] && set -x   # run pipeline with DEBUG=1 for verbose log
exec 2>&1  

### --- START: DEBUG HELPER FUNCTION ---
# This function will count the number of data records (lines without '#') in a VCF file.
count_vcf_records() {
  local vcf_file="$1"
  local description="$2"
  local count
  if [[ ! -f "$vcf_file" ]]; then
    echo "  [DEBUG] ✗ File not found for counting: $vcf_file"
    return
  fi
  # Use bcftools view -H to get only data lines, then count them.
  count=$(bcftools view -H "$vcf_file" | wc -l)
  echo "  [DEBUG] ✓ Records in $description ($vcf_file): $count"
}

### --- QC ANALYSIS FUNCTION ---
run_qc_analysis() {
  local vcf_file="$1"
  local qc_type="$2"  # "pre_imputation" or "post_imputation"
  local qc_dir="$3"
  
  echo "==> Running $qc_type QC analysis on: $vcf_file"
  
  if [[ ! -f "$vcf_file" ]]; then
    echo "  [QC] ✗ VCF file not found: $vcf_file"
    return 1
  fi
  
  mkdir -p "$qc_dir"
  local report_file="$qc_dir/${qc_type}_qc_report.txt"
  
  echo "=== $qc_type QC REPORT ===" > "$report_file"
  echo "VCF: $vcf_file" >> "$report_file"
  echo "Timestamp: $(date)" >> "$report_file"
  echo "" >> "$report_file"
  
  # Basic stats
  echo "=== BASIC STATS ===" >> "$report_file"
  bcftools stats "$vcf_file" >> "$report_file" 2>/dev/null || echo "Could not generate basic stats" >> "$report_file"
  
  # Sample count
  echo "=== SAMPLE INFORMATION ===" >> "$report_file"
  local sample_count=$(bcftools query -l "$vcf_file" | wc -l)
  echo "Number of samples: $sample_count" >> "$report_file"
  bcftools query -l "$vcf_file" >> "$report_file"
  
  # Variant count per chromosome
  echo "=== VARIANTS PER CHROMOSOME ===" >> "$report_file"
  bcftools query -f '%CHROM\n' "$vcf_file" | sort | uniq -c >> "$report_file"
  
  # Missing rate analysis
  echo "=== MISSING RATE ANALYSIS ===" >> "$report_file"
  local total_gts=$(bcftools query -f '[%GT\t]\n' "$vcf_file" | wc -w)
  local missing_gts=$(bcftools query -f '[%GT\t]\n' "$vcf_file" | grep -o '\./\.\|\./' | wc -l)
  local missing_rate=$(echo "scale=6; $missing_gts / $total_gts" | bc -l 2>/dev/null || echo "0")
  echo "Total genotypes: $total_gts" >> "$report_file"
  echo "Missing genotypes: $missing_gts" >> "$report_file"
  echo "Overall missing rate: $missing_rate" >> "$report_file"
  
  # Check phasing status
  echo "=== PHASING STATUS ===" >> "$report_file"
  local phased_count=$(bcftools query -f '[%GT\t]\n' "$vcf_file" | head -1000 | grep -o '|' | wc -l)
  local unphased_count=$(bcftools query -f '[%GT\t]\n' "$vcf_file" | head -1000 | grep -o '/' | wc -l)
  echo "Phased genotypes (sample of 1000 records): $phased_count" >> "$report_file"
  echo "Unphased genotypes (sample of 1000 records): $unphased_count" >> "$report_file"
  
  if [[ "$qc_type" == "post_imputation" ]]; then
    # DR2 analysis for post-imputation
    echo "=== IMPUTATION QUALITY (DR2) ANALYSIS ===" >> "$report_file"
    
    # Check if DR2 is declared in the header
     if bcftools view -h "$vcf_file" | grep -q 'ID=DR2,'; then
      # DR2 statistics
      bcftools query -f '%INFO/DR2\n' "$vcf_file" | awk '
      BEGIN { sum=0; count=0; high_qual=0; med_qual=0; low_qual=0 }
      {
        if ($1 != ".") {
          sum += $1
          count++
          if ($1 > 0.8) high_qual++
          else if ($1 > 0.3) med_qual++
          else low_qual++
        }
      }
      END {
        if (count > 0) {
          print "Total variants with DR2:", count
          print "Mean DR2:", sum/count
          print "Variants with DR2 > 0.8 (high quality):", high_qual, "(" int(high_qual/count*100) "%)"
          print "Variants with DR2 0.3-0.8 (medium quality):", med_qual, "(" int(med_qual/count*100) "%)" 
          print "Variants with DR2 < 0.3 (low quality):", low_qual, "(" int(low_qual/count*100) "%)"
        } else {
          print "No DR2 values found"
        }
      }' >> "$report_file"
      
      # DR2 distribution
      echo "=== DR2 DISTRIBUTION ===" >> "$report_file"
      bcftools query -f '%INFO/DR2\n' "$vcf_file" | awk '
      BEGIN { bin1=0; bin2=0; bin3=0; bin4=0; bin5=0; total=0 }
      {
        if ($1 != ".") {
          total++
          if ($1 <= 0.1) bin1++
          else if ($1 <= 0.3) bin2++
          else if ($1 <= 0.5) bin3++
          else if ($1 <= 0.8) bin4++
          else bin5++
        }
      }
      END {
        if (total > 0) {
          print "DR2 0.0-0.1:", bin1, "(" int(bin1/total*100) "%)"
          print "DR2 0.1-0.3:", bin2, "(" int(bin2/total*100) "%)"
          print "DR2 0.3-0.5:", bin3, "(" int(bin3/total*100) "%)"
          print "DR2 0.5-0.8:", bin4, "(" int(bin4/total*100) "%)"
          print "DR2 0.8-1.0:", bin5, "(" int(bin5/total*100) "%)"
        }
      }' >> "$report_file"
    else
      echo "No DR2 field found in VCF" >> "$report_file"
    fi
  fi
  
  echo "  [QC] ✓ $qc_type QC report saved: $report_file"
  
  # Display key findings to console
  echo "  [QC] Key findings:"
  echo "    Sample count: $sample_count"
  echo "    Missing rate: $missing_rate"
  if [[ $phased_count -gt $unphased_count ]]; then
    echo "    Phasing: Data appears phased ✓"
  else
    echo "    Phasing: Data appears unphased ⚠️"
  fi
  
  if [[ "$qc_type" == "post_imputation" ]] && bcftools query -f '%INFO/DR2\n' "$vcf_file" 2>/dev/null | head -1 | grep -q .; then
    local mean_dr2=$(bcftools query -f '%INFO/DR2\n' "$vcf_file" | awk '{if($1!="."){s+=$1;c++}} END{if(c>0) print s/c; else print 0}')
    echo "    Mean DR2: $mean_dr2"
    if (( $(echo "$mean_dr2 < 0.1" | bc -l 2>/dev/null || echo "0") )); then
      echo "    ⚠️  WARNING: Very low DR2 - poor imputation quality!"
    elif (( $(echo "$mean_dr2 < 0.3" | bc -l 2>/dev/null || echo "0") )); then
      echo "    ⚠️  CAUTION: Low DR2 - review imputation parameters"
    else
      echo "    ✓ DR2 looks reasonable"
    fi
  fi
}

### ───────────────────────────── 0. sanity ─────────────────────────────────###
[[ $# -eq 1 ]] || { echo "Usage: $0 <path/to/genome_file>.txt[.gz]"; exit 1; }
IN_TXT=$1
[[ -f $IN_TXT ]] || { echo "Input file not found: $IN_TXT"; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1"; exit 1; }; }
need bcftools; need tabix; need python; need java
need bgzip || true

debug_vcf () {
  local vcf=$1
  local idx
  if [[ $vcf == *.bcf ]]; then idx="${vcf}.csi"; else idx="${vcf}.tbi"; fi

  # ensure index exists & up-to-date
  if [[ ! -f $idx || $(stat -c %Y "$vcf") -gt $(stat -c %Y "$idx") ]]; then
    if [[ $vcf == *.bcf ]]; then
      bcftools index -f "$vcf"
    else
      tabix -f -p vcf "$vcf"
    fi
  fi

  echo "[DEBUG] $(basename "$vcf") – first 3 records & contig list"
  bcftools index -s "$vcf" | head
  bcftools view -H -r 1:1-5000 "$vcf" | head -n 3
}

### ─────────────────────────── 1. paths/vars ───────────────────────────────###
THREADS=${THREADS:-4}

ROOT_DIR="/app"

BASE=$(basename "$IN_TXT")
STEM=${BASE%%.*}

INPUT_BUILD="${INPUT_BUILD:-auto}" # override with INPUT_BUILD=38 ./impute_23andme.sh file.txt

FASTA_37="${ROOT_DIR}/genome_data/fasta/Homo_sapiens_assembly19.fasta"
FASTA_38="${ROOT_DIR}/genome_data/fasta/Homo_sapiens_assembly38.fasta"
CHAIN="${ROOT_DIR}/genome_data/chain/hg19ToHg38.over.chain.gz"

JAR="${ROOT_DIR}/genome_data/jars/beagle.27Feb25.75f.jar"
EAGLE="eagle"

REF_DIR_EAGLE="${ROOT_DIR}/genome_data/ref_bcfs_b38"
MAP_DIR_EAGLE="${MAP_DIR_EAGLE:-${ROOT_DIR}/genome_data/eagle_maps_b38}"

REF_DIR_BEAGLE="${ROOT_DIR}/genome_data/ref_brefs_b38"
MAP_DIR_BEAGLE="${ROOT_DIR}/genome_data/beagle_maps_b38"

# Debug: Check if key files exist
echo "==> Checking key file paths:"
echo "  CHAIN: $CHAIN $([ -f "$CHAIN" ] && echo "✓" || echo "✗")"
echo "  MAP_DIR_EAGLE: $MAP_DIR_EAGLE $([ -d "$MAP_DIR_EAGLE" ] && echo "✓" || echo "✗")"
if [ -d "$MAP_DIR_EAGLE" ]; then
  echo "  Files in MAP_DIR_EAGLE:"
  ls -la "$MAP_DIR_EAGLE/" | head -5
fi

# ---- per-sample folders & files ------------------------------------------- #
OUT_DIR="${ROOT_DIR}/users/${STEM}"
PHASED_DIR="${OUT_DIR}/phased_dir"
IMPUTED_DIR="${OUT_DIR}/imputed_dir"
QC_DIR="${OUT_DIR}/imputation_qc_reports"
mkdir -p "$OUT_DIR" "$PHASED_DIR" "$IMPUTED_DIR" "$QC_DIR"

RAW_GZ="${OUT_DIR}/${STEM}.txt.gz"

VCF_CHR_GZ="${OUT_DIR}/${STEM}.build37.chr.vcf.gz"
RAW_VCF="${OUT_DIR}/${STEM}.build37.chr.alt.vcf.gz"

LIFT_VCF="${OUT_DIR}/${STEM}.lift38.vcf"
PRIM_VCF="${OUT_DIR}/${STEM}.lift38.primary.vcf.gz"
CHR_VCF="${OUT_DIR}/${STEM}.lift38.primary.chr.vcf.gz"

FINAL_VCF="${OUT_DIR}/${STEM}_imputed_all.vcf.gz"

### ────────────────── 2. gzip copy + detect build ──────────────────────────###
echo "==> copy + gzip raw file"
if [[ $IN_TXT == *.gz ]]; then
  cp "$IN_TXT" "$RAW_GZ"
else
  gzip -c "$IN_TXT" > "$RAW_GZ"
fi

detect_build() {
  local hdr
  hdr=$(zgrep -m1 -Ei 'build|grch|hg[0-9]+' "$RAW_GZ" || true)
  if echo "$hdr" | grep -Eiq 'build[^0-9]?37|grch37|hg19|b37'; then echo 37
  elif echo "$hdr" | grep -Eiq 'build[^0-9]?38|grch38|hg38|b38'; then echo 38
  elif echo "$hdr" | grep -Eiq 'build[^0-9]?36|hg18|b36'; then echo 36
  else echo "unknown"
  fi
}

if [[ "$INPUT_BUILD" == "auto" ]]; then
  GUESS=$(detect_build)
  case "$GUESS" in
  37) INPUT_BUILD=37 ;;
  38) INPUT_BUILD=38 ;;
  36)
    echo "✗ Detected hg18/build36; not supported."
    exit 1
    ;;
  *)
    echo "⚠️  Could not detect build; defaulting to GRCh37. Override with INPUT_BUILD=38"
    INPUT_BUILD=37
    ;;
  esac
fi
echo "==> Detected/selected input build: GRCh${INPUT_BUILD}"

### ─────────────────── 3. TSV→VCF and (conditional) liftover ─────────────###
if [[ "$INPUT_BUILD" == 37 ]]; then

  # 1) TSV → chr-prefixed VCF (build37)
  echo "==> Pre-processing TSV to add 'chr' prefix..."
  TEMP_TSV_GZ="${OUT_DIR}/temp_with_chr.tsv.gz"
  zcat "$RAW_GZ" \
    | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$2="chr"$2; print}' \
    | bgzip -c > "$TEMP_TSV_GZ"
  
  echo "==> Converting TSV → VCF (build37)"
  bcftools convert --tsv2vcf "$TEMP_TSV_GZ" -f "$FASTA_37" -s "$STEM" -Oz -o "$VCF_CHR_GZ"

  rm "$TEMP_TSV_GZ"
  tabix -f -p vcf "$VCF_CHR_GZ"
  count_vcf_records "$VCF_CHR_GZ" "Initial VCF (chr-prefix)"

  # 2) QC raw conversion
  echo "==> QC: initial_conversion"
  run_qc_analysis "$VCF_CHR_GZ" "initial_conversion" "$QC_DIR"

  # 3) Patch missing ALT alleles
  echo "==> Patching missing ALT alleles"
  ALT_DB="${ROOT_DIR}/genome_data/alt_alleles.db"
  [[ -f "$ALT_DB" ]] || { echo "✗ Missing $ALT_DB"; exit 1; }
  python scripts/helpers/alt_fix.py --db "$ALT_DB" "$VCF_CHR_GZ"
  [[ -f "$RAW_VCF" ]] || { echo "✗ Expected $RAW_VCF"; exit 1; }
  count_vcf_records "$RAW_VCF" "VCF after alt_fix"

  # 4) Liftover to GRCh38 (still chr-prefixed)
  echo "==> CrossMap liftover (chr-prefix contigs)"
  CrossMap vcf "$CHAIN" "$RAW_VCF" "$FASTA_38" "$LIFT_VCF"
  lift_count=$(grep -vc '^#' "$LIFT_VCF")
  echo "  [DEBUG] Lifted records: $lift_count"

  # 5) Sort & index lifted VCF
  echo "==> Sorting lifted VCF (chr-prefix)"
  bcftools sort -Oz -o "$CHR_VCF" "$LIFT_VCF"
  tabix -f -p vcf "$CHR_VCF"
  count_vcf_records "$CHR_VCF" "Post-liftover VCF (chr-prefixed)"

  debug_vcf "$CHR_VCF"
  
  # 6) QC post-lift
  echo "==> QC: post_liftover"
  run_qc_analysis "$CHR_VCF" "post_liftover" "$QC_DIR"

  # --- 7) Replace contig header with hg38 lengths ---------------------------
  echo "==> Re-headering $CHR_VCF with hg38 contig lengths"
  REHEADERED="${CHR_VCF%.vcf.gz}.reheadered.vcf.gz"

  # -f pulls contig names & lengths from the .fai; output goes to stdout
  bcftools reheader -f "${FASTA_38}.fai" "$CHR_VCF" \
    | bcftools view -Oz -o "$REHEADERED" -
  
  debug_vcf "$REHEADERED"

  tabix -f -p vcf "$REHEADERED"
  mv "$REHEADERED"        "$CHR_VCF"
  mv "${REHEADERED}.tbi"  "$CHR_VCF.tbi"
  echo "  ✓ Header replaced and indexed: $CHR_VCF"

  # 8) Strip ambiguous SNPs before phasing
  echo "==> Removing strand-ambiguous (A/T, C/G) SNPs"
  FILTERED_VCF="${OUT_DIR}/${STEM}_stranded_clean.vcf.gz"
  bcftools view -e '( REF="A" && ALT="T" ) || ( REF="T" && ALT="A" ) || ( REF="C" && ALT="G" ) || ( REF="G" && ALT="C" )' \
    -Oz -o "$FILTERED_VCF" "$CHR_VCF"

  debug_vcf "$FILTERED_VCF"

  tabix -f -p vcf "$FILTERED_VCF"
  mv "$FILTERED_VCF" "$CHR_VCF"
  mv "${FILTERED_VCF}.tbi"  "$CHR_VCF.tbi"
  count_vcf_records "$CHR_VCF" "Strand-cleaned VCF"

  # 9) Final pre-imputation QC
  echo "==> QC: pre_imputation"
  run_qc_analysis "$CHR_VCF" "pre_imputation" "$QC_DIR"


  elif [[ "$INPUT_BUILD" == 38 ]]; then
    echo "==> Pre-processing TSV to add 'chr' prefix (build38)…"
    TEMP_TSV_GZ="${OUT_DIR}/temp_with_chr.tsv.gz"
    zcat "$RAW_GZ" | \
      awk 'BEGIN{OFS="\t"} /^#/ {print; next} !/^#/ {$2="chr"$2; print}' | \
      bgzip -c > "$TEMP_TSV_GZ"

    echo "==> TSV → VCF (build38) using pre-processed TSV"
    bcftools convert \
      --tsv2vcf "$TEMP_TSV_GZ" \
      -f "$FASTA_38" \
      -s "$STEM" \
      -Oz -o "$CHR_VCF" \
      2>"${OUT_DIR}/${STEM}.tsv2vcf.b38.log"
    rm "$TEMP_TSV_GZ"
    tabix -f -p vcf "$CHR_VCF"
    count_vcf_records "$CHR_VCF" "Initial VCF (build38)"

    echo "==> Running QC on initial build38 VCF…"
    run_qc_analysis "$CHR_VCF" "initial_conversion" "$QC_DIR"

    ALT_DB="${ROOT_DIR}/static_files/alt_alleles.db"
    echo "==> patch missing ALT alleles"
    python "${ROOT_DIR}/scripts/helpers/alt_fix.py" --db "$ALT_DB" "$CHR_VCF"
    count_vcf_records "${STEM}.alt.vcf.gz" "VCF after alt_fix (build38)"

    echo "==> Running QC on ALT-patched VCF…"
    run_qc_analysis "${STEM}.alt.vcf.gz" "post_alt_fix" "$QC_DIR"

    echo "==> Normalizing and selecting primary contigs"
    bcftools norm "$CHR_VCF" -Oz -o "$PRIM_VCF"
    tabix -f -p vcf "$PRIM_VCF"
    count_vcf_records "$PRIM_VCF" "Normalized primary-contig VCF"

    echo "==> Running QC on final pre-imputation VCF (build38)…"
    run_qc_analysis "$PRIM_VCF" "pre_imputation" "$QC_DIR"

  else
    echo "✗ Unexpected INPUT_BUILD=$INPUT_BUILD"
    exit 1
  fi


  ### ─────────────────────────── 4. Eagle phasing ────────────────────────────###
  for CHR in {1..22}; do
    REF_BCF="${REF_DIR_EAGLE}/1000GP_chr${CHR}.bcf"
    [[ -f "$REF_BCF" ]] || { echo "✗ Missing Eagle ref: $REF_BCF"; exit 1; }

    # 1) Prefix Eagle map with chr
    RAW_MAP="${MAP_DIR_EAGLE}/eagle_chr${CHR}_b38.map.gz"
    TMP_MAP="${OUT_DIR}/tmp_eagle_chr${CHR}.map.gz"

    # debug_vcf "$CHR_VCF"
    # debug_vcf "$REF_BCF"
    # zcat "$TMP_MAP" | head -5 | cat -vte
    # bcftools isec -n=2 -w1 -p "${OUT_DIR}/isec_chr${CHR}" "$CHR_VCF" "$REF_BCF"
    (
      # ---- header -------------------------------------------------------------
      echo "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)"

      # ---- 1. normalise & dedupe ---------------------------------------------
      zcat "$RAW_MAP" |
        awk -v CHR="$CHR" '
          # skip any header-ish lines (first field not all digits)
          $1 !~ /^[0-9]+$/ { next }

          # parse 3- or 4-column maps
          NF==3 { pos=$1; rate=$2; map=$3 }
          NF>=4 { pos=$2; rate=$3; map=$4 }

          # keep the first occurrence of each physical position
          !seen[pos]++ { print CHR, pos, rate, map }
        ' OFS=" " |

      # ---- 2. sort by physical position --------------------------------------
        LC_ALL=C sort -t' ' -k2,2n |

      # ---- 3. drop rows whose cM falls below previous -------------------------
        awk '
          NR==1 { prev_cM=$4; print; next }
          $4 >= prev_cM { print; prev_cM=$4 }
        '
    ) | bgzip -c > "$TMP_MAP"





    echo "[DEBUG] validating map order:"
    zcat "$TMP_MAP" |
      awk 'NR>1 {
            if ($2 < prev_pos || $4 < prev_cM)
                print "OUT-OF-ORDER:", $0
            prev_pos = $2
            prev_cM  = $4
          }'



    # 2) Run Eagle
    echo "==> Eagle phasing chr${CHR}"
    echo "==> DEBUG: CHR_VCF filepath ${CHR_VCF}"
    echo "=== DEBUG: contig label parity on chr${CHR} ==="
    echo "-- target (${CHR_VCF})"
    bcftools index -s "$CHR_VCF" | head -n 5
    echo "-- reference (${REF_BCF})"
    bcftools index -s "$REF_BCF" | head -n 5
    echo "-- variant intersection count"
    # bcftools isec -n=2 -w1 -p "${OUT_DIR}/isec_chr${CHR}" "$CHR_VCF" "$REF_BCF"
    # echo "  $(grep -vc '^#' "${OUT_DIR}"/isec_chr"${CHR}"/0003.vcf) shared sites"

    "$EAGLE" \
      --vcfRef         "$REF_BCF" \
      --vcfTarget      "$CHR_VCF" \
      --geneticMapFile "$TMP_MAP" \
      --chrom          "$CHR" \
      --outPrefix      "${PHASED_DIR}/${STEM}_phased_chr${CHR}" \
      --numThreads     "$THREADS" \
      --allowRefAltSwap

    rm "$TMP_MAP"
  done

    # 3) Post-phasing QC (sample chr21 & chr22)
   echo "==> QC: post_phasing"
   for CHR in 21 22; do
      run_qc_analysis \
        "${PHASED_DIR}/${STEM}_phased_chr${CHR}.vcf.gz" \
        "post_phasing_chr${CHR}" \
        "$QC_DIR"
  done

### ─────────────────────────── 5. Beagle imputation ────────────────────────###
for CHR in {1..22}; do 
  GT="${PHASED_DIR}/${STEM}_phased_chr${CHR}.vcf.gz"
  REF="${REF_DIR_BEAGLE}/1000GP_chr${CHR}.bref3"
  MAP="${MAP_DIR_BEAGLE}/beagle_chr${CHR}_b38.map"
  TEMP_MAP=$(mktemp) 

  trap 'rm -f "$TEMP_MAP"' EXIT HUP INT QUIT TERM

  if [[ ! -f "$GT" ]] || [[ ! -f "$REF" ]] || [[ ! -f "$MAP" ]]; then
    echo "==> skip chr${CHR} (missing one or more input files: VCF, REF, or MAP)"
    rm -f "$TEMP_MAP" 
    continue 
  fi

  echo "==> Preparing chr${CHR} for Beagle: adding 'chr' prefix to map file"
  awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' "$MAP" > "$TEMP_MAP"

  echo -e "\n==> Beagle imputation chr${CHR}"
  java -Xmx8g -jar "$JAR" \
    gt="$GT" \
    ref="$REF" \
    map="$TEMP_MAP" \
    window=40 overlap=2 ne=20000 \
    impute=true \
    gp=true \
    nthreads="$THREADS" \
    out="${IMPUTED_DIR}/${STEM}_imputed_chr${CHR}"

  rm -f "$TEMP_MAP"
done

trap 'rm -f "$LIST"' EXIT

echo "==> indexing and concatenating imputed chromosomes"
LIST=$(mktemp)
missing_flag=0

for CHR in {1..22}; do
  VC="${IMPUTED_DIR}/${STEM}_imputed_chr${CHR}.vcf.gz"
  if [[ ! -f "$VC" ]]; then
    echo "✗ Missing imputed VCF for chr${CHR}: $VC" >&2
    missing_flag=1
  else
    tabix -f -p vcf "$VC"
    echo "$VC" >> "$LIST"
  fi
done

if [[ $missing_flag -ne 0 ]]; then
  echo "⚠️  Aborting: one or more per-chr imputed VCFs are missing." >&2
  rm -f "$LIST"
  exit 1
fi

bcftools concat -f "$LIST" -Oz -o "$FINAL_VCF"
tabix -f -p vcf "$FINAL_VCF"
rm -f "$LIST"


# FINAL POST-IMPUTATION QC CHECK: Complete imputed VCF
echo "==> Running comprehensive QC on final imputed VCF..."
run_qc_analysis "$FINAL_VCF" "post_imputation" "$QC_DIR"

echo -e "\n✓ Imputation complete."
echo "  Final VCF: $FINAL_VCF"
echo "  QC Reports: $QC_DIR/"
echo ""
echo "  Final VCF: $FINAL_VCF"
echo "  QC Reports: $QC_DIR/"

echo "==> VERIFY: INFO lines in final VCF header"
bcftools view -h "$FINAL_VCF" | grep '^##INFO='

echo "==> VERIFY: first few DR2 values"
bcftools query -f '%INFO/DR2\n' "$FINAL_VCF" | head -n 10

bcftools +fill-tags "$FINAL_VCF" -- -t AF | \
  bcftools query -f '%INFO/AF\t%INFO/DR2\n' | \
  awk '
    {maf=$1<0.5?$1:1-$1; dr2=$2}
    maf<0.005 {rare+=1; rareDR+=dr2}
    maf>=0.05 {common+=1; commonDR+=dr2}
    END{
      print "Mean DR² common (MAF>=5%):", commonDR/common
      print "Mean DR² rare  (MAF<0.5%):", rareDR/rare
    }'


echo "==> QC SUMMARY:"
echo "  Check these files for detailed analysis:"
echo "    - Initial conversion QC: $QC_DIR/initial_conversion_qc_report.txt"
echo "    - Pre-imputation QC: $QC_DIR/pre_imputation_qc_report.txt"
echo "    - Post-phasing QC: $QC_DIR/post_phasing_chr1_qc_report.txt"
echo "    - Final imputation QC: $QC_DIR/post_imputation_qc_report.txt"
echo ""
echo "  Key things to check in the reports:"
echo "    1. Missing rates should be <5% at each step"
echo "    2. Data should be properly phased before imputation"
echo "    3. Mean DR2 should be >0.3 (ideally >0.5) in final output"
echo "    4. >30% of variants should have DR2 >0.8 for good quality"