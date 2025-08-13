#!/usr/bin/env bash
set -Eeuo pipefail

# ---------- Fixed paths ----------
VCF_PATTERN="/app/genome_data/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
LABELS_PATH="/app/scripts/helpers/integrated_call_samples_v3.20130502.ALL.panel"

# Work dirs
PCA_DIR="/app/pca_model"
TMP_DIR="${PCA_DIR}/tmp"
mkdir -p "${TMP_DIR}"

# ---------- Prepare keep list ----------
# Make a two-column keep file with FID=0 to match PLINK's default FID
# Skip header; strip CR; tab-delimited
awk 'BEGIN{FS=OFS="\t"} NR>1 { gsub("\r",""); print 0, $1 }' \
  "${LABELS_PATH}" > "${PCA_DIR}/labels.keep"

# Sanity: compute overlap between VCF IIDs and keep IIDs (2nd col)
VCF_SAMPLE_LIST="${TMP_DIR}/vcf_samples.txt"
bcftools query -l "${VCF_PATTERN/\{chr\}/1}" > "$VCF_SAMPLE_LIST"

overlap=$(grep -Fxf "$VCF_SAMPLE_LIST" <(cut -f2 "${PCA_DIR}/labels.keep") | wc -l)
echo "[INFO] Overlap between VCF and labels.keep: $overlap samples"
if [ "$overlap" -eq 0 ]; then
  echo "[ERROR] No matching IDs found between VCF and labels.keep (check whitespace/IDs)."
  exit 1
fi

# Exclude xMHC region (support both '6' and 'chr6')
MHC_FILE="${TMP_DIR}/mhc_exclude.bed1"
{
  printf "6\t25000000\t34000000\n"
  printf "chr6\t25000000\t34000000\n"
} > "$MHC_FILE"

# Per-chromosome processing
for c in {1..22}; do
  VCF_PATH="${VCF_PATTERN/\{chr\}/$c}"

  # 2) Convert VCF -> PGEN with filters + stable IDs; drop xMHC
  plink2 \
    --vcf "${VCF_PATH}" \
    --keep "${PCA_DIR}/labels.keep" \
    --autosome --snps-only just-acgt --max-alleles 2 \
    --maf 0.05 --geno 0.02 \
    --exclude range "${MHC_FILE}" \
    --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 200 truncate \
    --make-pgen --out "${TMP_DIR}/1kgp38_chr${c}"

  # 3) LD prune per chromosome
  plink2 \
    --pfile "${TMP_DIR}/1kgp38_chr${c}" \
    --indep-pairwise 50kb 1 0.2 \
    --out "${TMP_DIR}/chr${c}"
done

# 4) Combine pruned IDs across chromosomes, de-dup
cat ${TMP_DIR}/chr*.prune.in | sort -u > "${PCA_DIR}/panel.ids"

# 5) Build chr pos ref alt table 
{
  echo -e "chr\tpos\tref\talt"
  awk 'NR==FNR{keep[$1]=1;next} ($3 in keep){print $1"\t"$2"\t"$4"\t"$5}' \
    "${PCA_DIR}/panel.ids" ${TMP_DIR}/1kgp38_chr*.pvar \
  | awk 'BEGIN{OFS="\t"} $4!~/,/ {print $1,$2,$3,$4}'
} > "${PCA_DIR}/pca_sites.b38.tsv"


echo "[OK] Wrote ${PCA_DIR}/pca_sites.b38.tsv with $(wc -l < "${PCA_DIR}/pca_sites.b38.tsv") sites."
