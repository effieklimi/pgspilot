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
#     ├─ qc_reports/            # QC analysis reports
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
echo "  CHAIN: $CHAIN $([ -f "$CHAIN" ] && echo "✓" || echo "✗")"
echo "  MAP_DIR_EAGLE: $MAP_DIR_EAGLE $([ -d "$MAP_DIR_EAGLE" ] && echo "✓" || echo "✗")"

# ---- per-sample folders & files ------------------------------------------- #
OUT_DIR="${ROOT_DIR}/users/${STEM}"
PHASED_DIR="${OUT_DIR}/phased_dir"
IMPUTED_DIR="${OUT_DIR}/imputed_dir"
QC_DIR="${OUT_DIR}/qc_reports"
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

  # 1. TSV → chr-prefixed VCF (build37)
  echo "Pre-processing your genome"
  TEMP_TSV_GZ="${OUT_DIR}/temp_with_chr.tsv.gz"
  zcat "$RAW_GZ" \
    | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$2="chr"$2; print}' \
    | bgzip -c > "$TEMP_TSV_GZ"
  
  bcftools convert --tsv2vcf "$TEMP_TSV_GZ" -f "$FASTA_37" -s "$STEM" -Oz -o "$VCF_CHR_GZ"
  rm "$TEMP_TSV_GZ"
  tabix -f -p vcf "$VCF_CHR_GZ"

  # 2. Patch missing ALT alleles
  echo "Patching missing alternate alleles"
  ALT_DB="${ROOT_DIR}/genome_data/alt_alleles.db"
  [[ -f "$ALT_DB" ]] || { echo "✗ Missing $ALT_DB"; exit 1; }
  python "${ROOT_DIR}/scripts/helpers/alt_fix.py" --db "$ALT_DB" "$VCF_CHR_GZ"
  [[ -f "$RAW_VCF" ]] || { echo "✗ Expected $RAW_VCF"; exit 1; }

  # 3. Liftover to GRCh38 (still chr-prefixed)
  echo "Switching to the latest human genome assembly"
  CrossMap vcf "$CHAIN" "$RAW_VCF" "$FASTA_38" "$LIFT_VCF"
  # Sorting lifted VCF (chr-prefix):
  bcftools sort -Oz -o "$CHR_VCF" "$LIFT_VCF"
  tabix -f -p vcf "$CHR_VCF"
  debug_vcf "$CHR_VCF"
  run_qc_analysis "$CHR_VCF" "post_liftover" "$QC_DIR"
  # Replace contig header with hg38 lengths 
  # Re-headering $CHR_VCF with hg38 contig lengths
  REHEADERED="${CHR_VCF%.vcf.gz}.reheadered.vcf.gz"
  # -f pulls contig names & lengths from the .fai; output goes to stdout
  bcftools reheader -f "${FASTA_38}.fai" "$CHR_VCF" \
    | bcftools view -Oz -o "$REHEADERED" -
  debug_vcf "$REHEADERED"
  tabix -f -p vcf "$REHEADERED"
  mv "$REHEADERED"        "$CHR_VCF"
  mv "${REHEADERED}.tbi"  "$CHR_VCF.tbi"

  # 4. Strip ambiguous SNPs before phasing
  echo "Removing strand-ambiguous SNPs" # (A/T, C/G) 
  FILTERED_VCF="${OUT_DIR}/${STEM}_stranded_clean.vcf.gz"
  bcftools view -e '( REF="A" && ALT="T" ) || ( REF="T" && ALT="A" ) || ( REF="C" && ALT="G" ) || ( REF="G" && ALT="C" )' \
    -Oz -o "$FILTERED_VCF" "$CHR_VCF"
  debug_vcf "$FILTERED_VCF"
  tabix -f -p vcf "$FILTERED_VCF"
  mv "$FILTERED_VCF" "$CHR_VCF"
  mv "${FILTERED_VCF}.tbi"  "$CHR_VCF.tbi"
  run_qc_analysis "$CHR_VCF" "pre_imputation" "$QC_DIR"

elif [[ "$INPUT_BUILD" == 38 ]]; then
  echo "Pre-processing your genome"
  TEMP_TSV_GZ="${OUT_DIR}/temp_with_chr.tsv.gz"
  zcat "$RAW_GZ" | \
    awk 'BEGIN{OFS="\t"} /^#/ {print; next} !/^#/ {$2="chr"$2; print}' | \
    bgzip -c > "$TEMP_TSV_GZ"
  bcftools convert \
    --tsv2vcf "$TEMP_TSV_GZ" \
    -f "$FASTA_38" \
    -s "$STEM" \
    -Oz -o "$CHR_VCF" \
    2>"${OUT_DIR}/${STEM}.tsv2vcf.b38.log"
  rm "$TEMP_TSV_GZ"
  tabix -f -p vcf "$CHR_VCF"
  run_qc_analysis "$CHR_VCF" "initial_conversion" "$QC_DIR"

  # 2. Patch missing ALT alleles
  echo "Patching missing alternate alleles"
  ALT_DB="${ROOT_DIR}/genome_data/alt_alleles.db"
  [[ -f "$ALT_DB" ]] || { echo "Missing $ALT_DB"; exit 1; }
  python "${ROOT_DIR}/scripts/helpers/alt_fix.py" --db "$ALT_DB" "$CHR_VCF"
  [[ -f "$RAW_VCF" ]] || { echo "Expected $RAW_VCF"; exit 1; }

  # 3. Normalizing and selecting primary contigs
  echo "Preparation and quality control before phasing"
  bcftools norm "$CHR_VCF" -Oz -o "$PRIM_VCF"
  tabix -f -p vcf "$PRIM_VCF"
  run_qc_analysis "$PRIM_VCF" "pre_imputation" "$QC_DIR"

else
  echo "This genome build is not supported: $INPUT_BUILD"
  exit 1
fi


### ─────────────────────────── 4. Eagle phasing ────────────────────────────###
for CHR in {21..22}; do
  REF_BCF="${REF_DIR_EAGLE}/1000GP_chr${CHR}.bcf"
  [[ -f "$REF_BCF" ]] || { echo "✗ Missing Eagle ref: $REF_BCF"; exit 1; }

  # 1) Prefix Eagle map with chr
  RAW_MAP="${MAP_DIR_EAGLE}/eagle_chr${CHR}_b38.map.gz"
  TMP_MAP="${OUT_DIR}/tmp_eagle_chr${CHR}.map.gz"

  (
    # 1. normalise & dedupe 
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

    # 2. sort by physical position 
      LC_ALL=C sort -t' ' -k2,2n |

    # 3. drop rows whose cM falls below previous 
      awk '
        NR==1 { prev_cM=$4; print; next }
        $4 >= prev_cM { print; prev_cM=$4 }
      '
  ) | bgzip -c > "$TMP_MAP"

  # 2) Run Eagle
  echo "Phasing for chromosome ${CHR}"
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
echo "Post-phasing quality control"
for CHR in 21 22; do
  run_qc_analysis \
    "${PHASED_DIR}/${STEM}_phased_chr${CHR}.vcf.gz" \
    "post_phasing_chr${CHR}" \
    "$QC_DIR"
done

### ─────────────────────────── 5. Beagle imputation ────────────────────────###
for CHR in {21..22}; do 
  GT="${PHASED_DIR}/${STEM}_phased_chr${CHR}.vcf.gz"
  REF="${REF_DIR_BEAGLE}/1000GP_chr${CHR}.bref3"
  MAP="${MAP_DIR_BEAGLE}/beagle_chr${CHR}_b38.map"
  TEMP_MAP=$(mktemp) 

  trap 'rm -f "$TEMP_MAP"' EXIT HUP INT QUIT TERM

  if [[ ! -f "$GT" ]] || [[ ! -f "$REF" ]] || [[ ! -f "$MAP" ]]; then
    # skip chr${CHR} (missing one or more input files: VCF, REF, or MAP)
    rm -f "$TEMP_MAP" 
    continue 
  fi

  awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' "$MAP" > "$TEMP_MAP"

  echo -e "Imputation for chromosome ${CHR}"
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

echo "Indexing and concatenating imputed chromosomes"
LIST=$(mktemp)
missing_flag=0

for CHR in {21..22}; do
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


echo "Imputation complete"