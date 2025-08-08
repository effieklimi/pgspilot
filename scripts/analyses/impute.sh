#!/usr/bin/env bash
###############################################################################
# impute_23andme.sh – 23andMe text → GRCh38 phased & imputed VCF (chr1–22)
#
# USAGE:
#   ./impute_23andme.sh  <path/to/23andme_genome>.txt[.gz]
#
# OUTPUTS (per sample STEM from input filename):
#   ${ROOT_DIR}/users/${STEM}_results/
#     ├─ phased_dir/                      # Eagle output per chr
#     ├─ imputed_dir/                     # Beagle output per chr (with DS & GP)
#     ├─ qc_reports/                      # QC analysis reports
#     ├─ work/                            # transient intermediates (cleaner cwd)
#     └─ ${STEM}_imputed_all.vcf.gz (+ .tbi)
#
# NOTES:
#   • Input build auto-detected (b37|b38); b36/hg18 not supported.
#   • Uses 1000G high-coverage GRCh38 refs for Eagle and Beagle.
#   • DS (ALT dosage) + GP included.
###############################################################################
set -Eeuo pipefail
[[ "${DEBUG:-}" == 1 ]] && set -x
exec 2>&1

log()   { printf '%s %s\n' "$(date +'%F %T')" "$*"; }
debug() { [[ "${DEBUG:-0}" == 1 ]] && echo "[DEBUG] $*"; }

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1"; exit 1; }; }
need bcftools; need tabix; need bgzip; need python; need java; need CrossMap
# awk/grep/sort/wc/date come with coreutils; bc used for float compares if DEBUG
command -v bc >/dev/null 2>&1 || true

count_vcf_records() {
  [[ "${DEBUG:-0}" == 1 ]] || return 0
  local vcf="$1" desc="$2"
  [[ -f "$vcf" ]] || { echo "  [DEBUG] ✗ not found: $vcf"; return 0; }
  local n; n=$(bcftools view -H "$vcf" | wc -l)
  echo "  [DEBUG] $desc: $n variants"
}

run_qc_analysis() {
  # Writes a human-readable report; keeps console minimal.
  local vcf="$1" qc_type="$2" qc_dir="$3"
  [[ -f "$vcf" ]] || { echo "[QC] ✗ VCF not found: $vcf"; return 1; }
  mkdir -p "$qc_dir"
  local report="$qc_dir/${qc_type}_qc_report.txt"

  {
    echo "=== $qc_type QC REPORT ==="
    echo "VCF: $vcf"
    echo "Timestamp: $(date)"
    echo
    echo "=== BASIC STATS ==="
    bcftools stats "$vcf" 2>/dev/null || echo "Could not generate basic stats"

    echo
    echo "=== SAMPLE INFORMATION ==="
    local samples sample_count
    samples=$(bcftools query -l "$vcf")
    sample_count=$(printf '%s\n' "$samples" | wc -l | tr -d ' ')
    echo "Number of samples: $sample_count"
    printf '%s\n' "$samples"

    echo
    echo "=== VARIANTS PER CHROMOSOME ==="
    bcftools query -f '%CHROM\n' "$vcf" | sort | uniq -c

    echo
    echo "=== MISSING RATE ANALYSIS ==="
    local total_gts missing_gts missing_rate
    total_gts=$(bcftools query -f '[%GT\t]\n' "$vcf" | wc -w | tr -d ' ')
    missing_gts=$(bcftools query -f '[%GT\t]\n' "$vcf" | grep -oE '\./\.|\.\/' | wc -l | tr -d ' ')
    if command -v bc >/dev/null 2>&1; then
      missing_rate=$(echo "scale=6; $missing_gts / ($total_gts+0.000001)" | bc -l)
    else
      missing_rate="n/a (bc not present)"
    fi
    echo "Total genotypes: $total_gts"
    echo "Missing genotypes: $missing_gts"
    echo "Overall missing rate: $missing_rate"

    echo
    echo "=== PHASING STATUS (sample of 1000 records) ==="
    local phased_count unphased_count
    phased_count=$(bcftools query -f '[%GT\t]\n' "$vcf" | head -1000 | grep -o '|' | wc -l | tr -d ' ')
    unphased_count=$(bcftools query -f '[%GT\t]\n' "$vcf" | head -1000 | grep -o '/' | wc -l | tr -d ' ')
    echo "Phased genotypes: $phased_count"
    echo "Unphased genotypes: $unphased_count"

    if [[ "$qc_type" == "post_imputation" ]]; then
      echo
      echo "=== IMPUTATION QUALITY (DR2) ==="
      if bcftools view -h "$vcf" | grep -q 'ID=DR2,'; then
        bcftools query -f '%INFO/DR2\n' "$vcf" | awk '
          BEGIN { sum=0; c=0; h=0; m=0; l=0 }
          $1!="." { d=$1; sum+=d; c++; if(d>0.8)h++; else if(d>0.3)m++; else l++; }
          END {
            if(c>0){
              printf("Total variants with DR2: %d\nMean DR2: %.6f\n", c, sum/c);
              printf(">0.8: %d (%.0f%%)\n0.3–0.8: %d (%.0f%%)\n<0.3: %d (%.0f%%)\n",
                     h, 100*h/c, m, 100*m/c, l, 100*l/c);
            } else { print "No DR2 values found" }
          }'
        echo
        echo "=== DR2 DISTRIBUTION ==="
        bcftools query -f '%INFO/DR2\n' "$vcf" | awk '
          BEGIN{b1=b2=b3=b4=b5=0; t=0}
          $1!="." {d=$1; t++; if(d<=0.1)b1++; else if(d<=0.3)b2++; else if(d<=0.5)b3++; else if(d<=0.8)b4++; else b5++}
          END{
            if(t>0){
              printf("0.0–0.1: %d (%.0f%%)\n0.1–0.3: %d (%.0f%%)\n0.3–0.5: %d (%.0f%%)\n0.5–0.8: %d (%.0f%%)\n0.8–1.0: %d (%.0f%%)\n",
                     b1,100*b1/t,b2,100*b2/t,b3,100*b3/t,b4,100*b4/t,b5,100*b5/t)
            }
          }'
      else
        echo "No DR2 field found in VCF"
      fi
    fi
  } >"$report"

  log "[QC] ${qc_type} → $report"
}

### ─────────────────────────── 0) Args / paths ──────────────────────────────
[[ $# -eq 1 ]] || { echo "Usage: $0 <path/to/genome_file>.txt[.gz]"; exit 1; }
IN_TXT=$1
[[ -f $IN_TXT ]] || { echo "Input file not found: $IN_TXT"; exit 1; }

THREADS=${THREADS:-4}
ROOT_DIR="/app"

BASE=$(basename "$IN_TXT")
STEM=${BASE%%.*}

INPUT_BUILD="${INPUT_BUILD:-auto}"      # override: INPUT_BUILD=38 ./impute_23andme.sh file.txt

FASTA_37="${ROOT_DIR}/genome_data/fasta/Homo_sapiens_assembly19.fasta"
FASTA_38="${ROOT_DIR}/genome_data/fasta/Homo_sapiens_assembly38.fasta"
CHAIN="${ROOT_DIR}/genome_data/chain/hg19ToHg38.over.chain.gz"
ALT_DB="${ROOT_DIR}/genome_data/alt_alleles.db"

JAR="${ROOT_DIR}/genome_data/jars/beagle.27Feb25.75f.jar"
EAGLE="eagle"

REF_DIR_EAGLE="${ROOT_DIR}/genome_data/ref_bcfs_b38"
MAP_DIR_EAGLE="${MAP_DIR_EAGLE:-${ROOT_DIR}/genome_data/eagle_maps_b38}"

REF_DIR_BEAGLE="${ROOT_DIR}/genome_data/ref_brefs_b38"
MAP_DIR_BEAGLE="${ROOT_DIR}/genome_data/beagle_maps_b38"

# per-sample folders
OUT_DIR="${ROOT_DIR}/users/${STEM}_results"
PHASED_DIR="${OUT_DIR}/phased_dir"
IMPUTED_DIR="${OUT_DIR}/imputed_dir"
QC_DIR="${OUT_DIR}/qc_reports"
WORK_DIR="${OUT_DIR}/work"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$PHASED_DIR" "$IMPUTED_DIR" "$QC_DIR" "$WORK_DIR" "$LOG_DIR"

FINAL_VCF="${OUT_DIR}/${STEM}_imputed_all.vcf.gz"

# working/intermediate filenames live in WORK_DIR to keep things tidy
RAW_GZ="${WORK_DIR}/${STEM}.txt.gz"
VCF_CHR_GZ="${WORK_DIR}/${STEM}.build37.chr.vcf.gz"
ALT37_VCF="${WORK_DIR}/${STEM}.build37.chr.alt.vcf.gz"
LIFT_VCF="${WORK_DIR}/${STEM}.lift38.vcf"
PRIM_VCF="${WORK_DIR}/${STEM}.lift38.primary.vcf.gz"
CHR_VCF="${WORK_DIR}/${STEM}.lift38.primary.chr.vcf.gz"   # unified pre-imputation VCF (GRCh38, chr-prefixed)

### ───────────────────────── 1) gzip + build detect ─────────────────────────
log "Copying & gzipping input"
if [[ $IN_TXT == *.gz ]]; then cp "$IN_TXT" "$RAW_GZ"; else gzip -c "$IN_TXT" > "$RAW_GZ"; fi

detect_build() {
  local hdr
  hdr=$(zgrep -m1 -Ei 'build|grch|hg[0-9]+' "$RAW_GZ" || true)
  if   echo "$hdr" | grep -Eiq 'build[^0-9]?37|grch37|hg19|b37'; then echo 37
  elif echo "$hdr" | grep -Eiq 'build[^0-9]?38|grch38|hg38|b38'; then echo 38
  elif echo "$hdr" | grep -Eiq 'build[^0-9]?36|hg18|b36'; then echo 36
  else echo "unknown"; fi
}

if [[ "$INPUT_BUILD" == "auto" ]]; then
  case "$(detect_build)" in
    37) INPUT_BUILD=37 ;;
    38) INPUT_BUILD=38 ;;
    36) echo "Detected hg18/build36; not supported."; exit 1 ;;
    *)  echo "Build not detected; defaulting to GRCh37 (override with INPUT_BUILD=38)"; INPUT_BUILD=37 ;;
  esac
fi
log "Selected input build: GRCh${INPUT_BUILD}"

### ─────────────────────── 2) TSV→VCF & liftover/normalize ──────────────────
if [[ "$INPUT_BUILD" == 37 ]]; then
  log "TSV → VCF on GRCh37 (with chr prefix)"
  TMP_TSV="${WORK_DIR}/${STEM}.tmp_with_chr.tsv.gz"
  zcat "$RAW_GZ" | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$2="chr"$2; print}' | bgzip -c > "$TMP_TSV"
  bcftools convert --tsv2vcf "$TMP_TSV" -f "$FASTA_37" -s "$STEM" -Oz -o "$VCF_CHR_GZ"
  rm -f "$TMP_TSV"
  tabix -f -p vcf "$VCF_CHR_GZ"
  count_vcf_records "$VCF_CHR_GZ" "Initial (b37 chr-prefixed)"
  run_qc_analysis "$VCF_CHR_GZ" "initial_conversion" "$QC_DIR"

  log "Patching missing ALT alleles"
  [[ -f "$ALT_DB" ]] || { echo "Missing ALT DB: $ALT_DB"; exit 1; }
  python scripts/helpers/alt_fix.py --db "$ALT_DB" "$VCF_CHR_GZ"
  [[ -f "$ALT37_VCF" ]] || { echo "Expected output not found: $ALT37_VCF"; exit 1; }
  count_vcf_records "$ALT37_VCF" "After alt_fix"

  log "Liftover to GRCh38 (CrossMap)"
  CrossMap vcf "$CHAIN" "$ALT37_VCF" "$FASTA_38" "$LIFT_VCF"

  log "Sorting & indexing lifted VCF"
  bcftools sort -Oz -o "$CHR_VCF" "$LIFT_VCF"
  tabix -f -p vcf "$CHR_VCF"
  count_vcf_records "$CHR_VCF" "Post-liftover"

  run_qc_analysis "$CHR_VCF" "post_liftover" "$QC_DIR"

  log "Reheader with GRCh38 contig lengths"
  local_reheader="${WORK_DIR}/$(basename "${CHR_VCF%.vcf.gz}").reheadered.vcf.gz"
  bcftools reheader -f "${FASTA_38}.fai" "$CHR_VCF" | bcftools view -Oz -o "$local_reheader" -
  tabix -f -p vcf "$local_reheader"
  mv "$local_reheader" "$CHR_VCF"; mv "${local_reheader}.tbi" "${CHR_VCF}.tbi"

  log "Removing strand-ambiguous SNPs (A/T, C/G)"
  filtered="${WORK_DIR}/${STEM}.b38.primary.chr.stranded_clean.vcf.gz"
  bcftools view -e '( REF="A" && ALT="T" ) || ( REF="T" && ALT="A" ) || ( REF="C" && ALT="G" ) || ( REF="G" && ALT="C" )' -Oz -o "$filtered" "$CHR_VCF"
  tabix -f -p vcf "$filtered"
  mv "$filtered" "$CHR_VCF"; mv "${filtered}.tbi" "${CHR_VCF}.tbi"
  count_vcf_records "$CHR_VCF" "Strand-cleaned (pre-imputation)"
  run_qc_analysis "$CHR_VCF" "pre_imputation" "$QC_DIR"

elif [[ "$INPUT_BUILD" == 38 ]]; then
  log "TSV → VCF on GRCh38 (with chr prefix)"
  TMP_TSV="${WORK_DIR}/${STEM}.tmp_with_chr.tsv.gz"
  zcat "$RAW_GZ" | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$2="chr"$2; print}' | bgzip -c > "$TMP_TSV"
  CHR_VCF="${WORK_DIR}/${STEM}.build38.primary.chr.vcf.gz"
  bcftools convert --tsv2vcf "$TMP_TSV" -f "$FASTA_38" -s "$STEM" -Oz -o "$CHR_VCF" 2>"${LOG_DIR}/${STEM}.tsv2vcf.b38.log"
  rm -f "$TMP_TSV"
  tabix -f -p vcf "$CHR_VCF"
  count_vcf_records "$CHR_VCF" "Initial (b38 chr-prefixed)"
  run_qc_analysis "$CHR_VCF" "initial_conversion" "$QC_DIR"

  log "Patching missing ALT alleles"
  [[ -f "$ALT_DB" ]] || { echo "Missing ALT DB: $ALT_DB"; exit 1; }
  python "${ROOT_DIR}/scripts/helpers/alt_fix.py" --db "$ALT_DB" "$CHR_VCF"
  ALT38_VCF="${CHR_VCF%.vcf.gz}.alt.vcf.gz"
  [[ -f "$ALT38_VCF" ]] && CHR_VCF="$ALT38_VCF"
  count_vcf_records "$CHR_VCF" "After alt_fix (b38)"
  run_qc_analysis "$CHR_VCF" "post_alt_fix" "$QC_DIR"

  log "Normalize and select primary contigs"
  PRIM_VCF="${WORK_DIR}/${STEM}.build38.primary.vcf.gz"
  bcftools norm "$CHR_VCF" -Oz -o "$PRIM_VCF"
  tabix -f -p vcf "$PRIM_VCF"
  count_vcf_records "$PRIM_VCF" "Normalized primary-contig"
  run_qc_analysis "$PRIM_VCF" "pre_imputation" "$QC_DIR"
  CHR_VCF="$PRIM_VCF"
else
  echo "Unexpected INPUT_BUILD=$INPUT_BUILD"; exit 1
fi

### ─────────────────────────── 3) Eagle phasing ─────────────────────────────
for CHR in {1..22}; do
  REF_BCF="${REF_DIR_EAGLE}/1000GP_chr${CHR}.bcf"
  [[ -f "$REF_BCF" ]] || { echo "Missing Eagle ref: $REF_BCF"; exit 1; }

  RAW_MAP="${MAP_DIR_EAGLE}/eagle_chr${CHR}_b38.map.gz"
  TMP_MAP="${WORK_DIR}/eagle_chr${CHR}.map.gz"

  # sanitize & prefix map (monotonic cM, dedupe)
  (
    echo "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)"
    zcat "$RAW_MAP" | awk -v CHR="$CHR" '
      $1 !~ /^[0-9]+$/ { next }
      NF==3 { pos=$1; rate=$2; map=$3 }
      NF>=4 { pos=$2; rate=$3; map=$4 }
      !seen[pos]++ { print CHR, pos, rate, map }
    ' OFS=" " | LC_ALL=C sort -t' ' -k2,2n | awk '
      NR==1 { prev=$4; print; next }
      $4>=prev { print; prev=$4 }
    '
  ) | bgzip -c > "$TMP_MAP"

  log "Eagle phasing chr${CHR}"
  "$EAGLE" \
    --vcfRef         "$REF_BCF" \
    --vcfTarget      "$CHR_VCF" \
    --geneticMapFile "$TMP_MAP" \
    --chrom          "$CHR" \
    --outPrefix      "${PHASED_DIR}/${STEM}_phased_chr${CHR}" \
    --numThreads     "$THREADS" \
    --allowRefAltSwap

  rm -f "$TMP_MAP"
done

# quick post-phasing QC on two small chromosomes
for CHR in 21 22; do
  run_qc_analysis "${PHASED_DIR}/${STEM}_phased_chr${CHR}.vcf.gz" "post_phasing_chr${CHR}" "$QC_DIR"
done

### ───────────────────────── 4) Beagle imputation ───────────────────────────
for CHR in {1..22}; do
  GT="${PHASED_DIR}/${STEM}_phased_chr${CHR}.vcf.gz"
  REF="${REF_DIR_BEAGLE}/1000GP_chr${CHR}.bref3"
  MAP="${MAP_DIR_BEAGLE}/beagle_chr${CHR}_b38.map"
  [[ -f "$GT" && -f "$REF" && -f "$MAP" ]] || { log "Skipping chr${CHR} (missing GT/REF/MAP)"; continue; }

  # add chr prefix to MAP into a temp file
  TMPMAP="$(mktemp "${WORK_DIR}/beagle_map_chr${CHR}.XXXX")"
  awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' "$MAP" > "$TMPMAP"

  log "Beagle imputation chr${CHR}"
  java -Xmx8g -jar "$JAR" \
    gt="$GT" ref="$REF" map="$TMPMAP" \
    window=40 overlap=2 ne=20000 \
    impute=true gp=true nthreads="$THREADS" \
    out="${IMPUTED_DIR}/${STEM}_imputed_chr${CHR}"

  rm -f "$TMPMAP"
done

### ────────────────────── 5) Concatenate & final QC ─────────────────────────
log "Indexing & concatenating imputed chromosomes"
LIST="$(mktemp "${WORK_DIR}/imputed_list.XXXX")"; trap 'rm -f "$LIST"' EXIT
missing=0
for CHR in {1..22}; do
  VC="${IMPUTED_DIR}/${STEM}_imputed_chr${CHR}.vcf.gz"
  if [[ -f "$VC" ]]; then
    tabix -f -p vcf "$VC"
    echo "$VC" >> "$LIST"
  else
    echo "Missing imputed VCF for chr${CHR}: $VC" >&2
    missing=1
  fi
done
[[ $missing -eq 0 ]] || { echo "Aborting: one or more imputed VCFs missing."; exit 1; }

bcftools concat -f "$LIST" -Oz -o "$FINAL_VCF"
tabix -f -p vcf "$FINAL_VCF"
rm -f "$LIST"

log "Comprehensive QC on final imputed VCF"
run_qc_analysis "$FINAL_VCF" "post_imputation" "$QC_DIR"

# Quick DR2 sanity by MAF bins (console)
bcftools +fill-tags "$FINAL_VCF" -- -t AF | \
  bcftools query -f '%INFO/AF\t%INFO/DR2\n' | \
  awk '
    {maf=$1<0.5?$1:1-$1; dr2=$2}
    maf<0.005 {rare+=1; rareDR+=dr2}
    maf>=0.05 {common+=1; commonDR+=dr2}
    END{
      if(common>0) printf("Mean DR² common (MAF≥5%%): %.4f\n", commonDR/common);
      if(rare>0)   printf("Mean DR² rare   (MAF<0.5%%): %.4f\n", rareDR/rare);
    }'

echo
echo "✓ Imputation complete."
echo "  Final VCF : $FINAL_VCF"
echo "  QC reports: $QC_DIR/"
