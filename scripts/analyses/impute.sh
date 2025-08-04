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
#     └─ ${STEM}_imputed_all.vcf.gz (+ .tbi)
#
# NOTES:
#   • Input build auto-detected (b37|b38); b36/hg18 not supported here.
#   • Uses 1000G high-coverage GRCh38 refs for Eagle and Beagle.
#   • DS (dosage of ALT) is included for scoring; GP also kept (fallback).
#   • Ancestry-specific PRS standardization is NOT done here.
###############################################################################
set -Eeuo pipefail

### ───────────────────────────── 0. sanity ─────────────────────────────────###
[[ $# -eq 1 ]] || { echo "Usage: $0 <path/to/genome_file>.txt[.gz]"; exit 1; }
IN_TXT=$1
[[ -f $IN_TXT ]] || { echo "Input file not found: $IN_TXT"; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1"; exit 1; }; }
need bcftools; need tabix; need python; need java
need bgzip || true

### ─────────────────────────── 1. paths/vars ───────────────────────────────###
THREADS=${THREADS:-4}
ROOT_DIR=$(pwd)

BASE=$(basename "$IN_TXT")
STEM=${BASE%%.*}

INPUT_BUILD="${INPUT_BUILD:-auto}"   # override with INPUT_BUILD=38 ./impute_23andme.sh file.txt
TARGET_BUILD=38

# ---- static resources (set these to your actual paths) -------------------- #
FASTA_37="${ROOT_DIR}/genome_data/fasta/Homo_sapiens_assembly19.fasta"
FASTA_38="${ROOT_DIR}/genome_data/fasta/Homo_sapiens_assembly38.fasta"
CHAIN="${ROOT_DIR}/genome_data/chain/hg19ToHg38.over.chain.gz"

JAR="${ROOT_DIR}/genome_data/jars/beagle.27Feb25.75f.jar"
EAGLE="eagle"

REF_DIR_EAGLE="${ROOT_DIR}/genome_data/ref_bcfs_b38"
MAP_DIR_EAGLE="${ROOT_DIR}/genome_data/eagle_maps_b38"

REF_DIR_BEAGLE="${ROOT_DIR}/genome_data/ref_brefs_b38"
MAP_DIR_BEAGLE="${ROOT_DIR}/genome_data/beagle_maps_b38"

ADDCHR_MAP="${ROOT_DIR}/scripts/helpers/addchr.txt"   # e.g., "1  chr1" ..."MT chrM"
DROPCHR_MAP="${ROOT_DIR}/scripts/helpers/dropchr.txt" # e.g., "chr1 1" ..."chrM MT"

# ---- per-sample folders & files ------------------------------------------- #
OUT_DIR="${ROOT_DIR}/users/${STEM}_results"
PHASED_DIR="${OUT_DIR}/phased_dir"
IMPUTED_DIR="${OUT_DIR}/imputed_dir"
mkdir -p "$OUT_DIR" "$PHASED_DIR" "$IMPUTED_DIR"

RAW_GZ="${OUT_DIR}/${STEM}.txt.gz"

VCF_GZ="${OUT_DIR}/${STEM}.build37.vcf.gz"
VCF_CHR_GZ="${OUT_DIR}/${STEM}.build37.chr.vcf.gz"
RAW_VCF="${OUT_DIR}/${STEM}.build37.chr.alt.vcf.gz"

NOCHR_VCF="${OUT_DIR}/${STEM}.nochr.b37.vcf.gz"
LIFT_VCF="${OUT_DIR}/${STEM}.lift38.vcf"
SORT_VCF="${OUT_DIR}/${STEM}.lift38.sorted.vcf.gz"
NORM_VCF="${OUT_DIR}/${STEM}.lift38.norm.vcf.gz"
FIXREF_VCF="${OUT_DIR}/${STEM}.lift38.norm.fixref.vcf.gz"
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
  if   echo "$hdr" | grep -Eiq 'build[^0-9]?37|grch37|hg19|b37'; then echo 37
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
    36) echo "✗ Detected hg18/build36; not supported."; exit 1 ;;
    *)  echo "⚠️  Could not detect build; defaulting to GRCh37. Override with INPUT_BUILD=38"; INPUT_BUILD=37 ;;
  esac
fi
echo "==> Detected/selected input build: GRCh${INPUT_BUILD}"

### ─────────────────── 3. TSV→VCF and (conditional) liftover ─────────────###
if [[ "$INPUT_BUILD" == 37 ]]; then
  echo "==> TSV → VCF (build37)"
  # Use simple bcftools convert without column specification - let it auto-detect 23andMe format
  bcftools convert --tsv2vcf "$RAW_GZ" -f "$FASTA_37" -s "$STEM" -Oz -o "$VCF_GZ" \
    2> "${OUT_DIR}/${STEM}.tsv2vcf.b37.log"
  tabix -f -p vcf "$VCF_GZ"

  echo "==> add chr-prefixes"
  bcftools annotate --rename-chrs "$ADDCHR_MAP" -Oz -o "$VCF_CHR_GZ" "$VCF_GZ"
  tabix -f -p vcf "$VCF_CHR_GZ"

  echo "==> patch missing ALT alleles"
  python scripts/helpers/alt_fix.py "$VCF_CHR_GZ"
  [[ -f "$RAW_VCF" ]] || { echo "Expected $RAW_VCF after alt_fix.py, not found"; exit 1; }

  echo "==> chr→nochr rename (for CrossMap)"
  bcftools annotate --rename-chrs "$DROPCHR_MAP" -Oz -o "$NOCHR_VCF" "$RAW_VCF"
  tabix -f -p vcf "$NOCHR_VCF"

  echo "==> CrossMap liftover to GRCh38"
  crossmap vcf "$CHAIN" "$NOCHR_VCF" "$FASTA_38" "$LIFT_VCF"

  echo "==> sort + bgzip"
  bcftools sort "$LIFT_VCF" -Oz -o "$SORT_VCF"
  tabix -f -p vcf "$SORT_VCF"

  echo "==> keep primary contigs"
  CONTIGS=$(printf '%s,' {1..22} X Y MT); CONTIGS=${CONTIGS%,}
  bcftools view -r "$CONTIGS" -Oz -o "$PRIM_VCF" "$SORT_VCF"
  tabix -f -p vcf "$PRIM_VCF"

  echo "==> add chr-prefix again"
  bcftools annotate --rename-chrs "$ADDCHR_MAP" -Oz -o "$CHR_VCF" "$PRIM_VCF"
  tabix -f -p vcf "$CHR_VCF"

elif [[ "$INPUT_BUILD" == 38 ]]; then
  echo "==> TSV → VCF (build38)"
  VCF38_GZ="${OUT_DIR}/${STEM}.build38.vcf.gz"
  VCF38_CHR_GZ="${OUT_DIR}/${STEM}.build38.chr.vcf.gz"
  RAW38_VCF="${OUT_DIR}/${STEM}.build38.chr.alt.vcf.gz"

  # Use simple bcftools convert without column specification
  bcftools convert --tsv2vcf "$RAW_GZ" -f "$FASTA_38" -s "$STEM" -Oz -o "$VCF38_GZ" \
    2> "${OUT_DIR}/${STEM}.tsv2vcf.b38.log"
  tabix -f -p vcf "$VCF38_GZ"

  echo "==> add chr-prefixes (b38)"
  bcftools annotate --rename-chrs "$ADDCHR_MAP" -Oz -o "$VCF38_CHR_GZ" "$VCF38_GZ"
  tabix -f -p vcf "$VCF38_CHR_GZ"

  echo "==> patch missing ALT alleles"
  python scripts/helpers/alt_fix.py "$VCF38_CHR_GZ"
  [[ -f "$RAW38_VCF" ]] || { echo "Expected $RAW38_VCF after alt_fix.py, not found"; exit 1; }

  echo "==> normalize + fixref + keep primary contigs"
  bcftools norm -f "$FASTA_38" -m -both "$RAW38_VCF" -Oz -o "$NORM_VCF" && tabix -f -p vcf "$NORM_VCF"
  bcftools +fixref "$NORM_VCF" -- -f "$FASTA_38" -m flip -Oz -o "$FIXREF_VCF" && tabix -f -p vcf "$FIXREF_VCF"
  CONTIGS=$(printf 'chr%s,' {1..22} X Y M); CONTIGS=${CONTIGS%,}
  bcftools view -r "$CONTIGS" -Oz -o "$PRIM_VCF" "$FIXREF_VCF" && tabix -f -p vcf "$PRIM_VCF"
  CHR_VCF="$PRIM_VCF"
else
  echo "✗ Unexpected INPUT_BUILD=$INPUT_BUILD"; exit 1
fi

### ─────────────────────────── 4. Eagle phasing ────────────────────────────###
for CHR in {1..2}; do # for CHR in {1..22}; do
  REF_BCF="${REF_DIR_EAGLE}/1000GP_chr${CHR}.bcf"
  MAP="${MAP_DIR_EAGLE}/eagle_chr${CHR}_b38.map.gz"
  [[ -f $REF_BCF ]] || { echo "Missing Eagle ref: $REF_BCF"; exit 1; }
  [[ -f $MAP ]] || { echo "Missing Eagle map: $MAP"; exit 1; }

  echo -e "\n==> Eagle phasing chr${CHR}"
  "$EAGLE" \
    --vcfRef         "$REF_BCF" \
    --vcfTarget      "$CHR_VCF" \
    --geneticMapFile "$MAP" \
    --chrom          "$CHR" \
    --outPrefix      "${PHASED_DIR}/${STEM}_phased_chr${CHR}" \
    --numThreads     "$THREADS" \
    --allowRefAltSwap
done

### ─────────────────────────── 5. Beagle imputation ────────────────────────###
for CHR in {1..2}; do # for CHR in {1..22}; do
  GT="${PHASED_DIR}/${STEM}_phased_chr${CHR}.vcf.gz"
  REF="${REF_DIR_BEAGLE}/1000GP_chr${CHR}.bref3"
  MAP="${MAP_DIR_BEAGLE}/beagle_chr${CHR}_b38.map"

  for f in "$GT" "$REF" "$MAP"; do
    [[ -f $f ]] || { echo "==> skip chr${CHR} (missing $f)"; continue 2; }
  done

  echo -e "\n==> Beagle imputation chr${CHR}"
  java -Xmx8g -jar "$JAR" \
    gt="$GT" \
    ref="$REF" \
    map="$MAP" \
    impute=true \
    gp=true \
    dose=true \
    nthreads="$THREADS" \
    out="${IMPUTED_DIR}/${STEM}_imputed_chr${CHR}"
done

### ─────────────────────── 6. concatenate chr1-22 ─────────────────────────###
echo '==> indexing and concatenating imputed chromosomes'
LIST=$(mktemp); trap 'rm -f "$LIST"' EXIT
for CHR in {1..22}; do
  VC="${IMPUTED_DIR}/${STEM}_imputed_chr${CHR}.vcf.gz"
  [[ -f $VC ]] || { echo "Missing $VC; was this chromosome imputed?"; continue; }
  tabix -f -p vcf "$VC"
  echo "$VC" >> "$LIST"
done

bcftools concat -f "$LIST" -a -Oz -o "$FINAL_VCF"
tabix -f -p vcf "$FINAL_VCF"

echo -e "\n✓ Imputation complete."
echo "  Final VCF: $FINAL_VCF"