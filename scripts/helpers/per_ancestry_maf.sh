#!/usr/bin/env bash
set -Eeuo pipefail

# Fixed paths
VCF_PATTERN="/app/genome_data/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
LABELS_PATH="/app/scripts/helpers/integrated_call_samples_v3.20130502.ALL.panel"

# Work dirs
PCA_DIR="/app/pca_model"
TMP_DIR="${PCA_DIR}/tmp"
mkdir -p "${TMP_DIR}"

# Prepare keep list
# Build per-superpop keep lists (EUR/AFR/EAS/AMR/SAS) only if missing
for SP in AFR AMR EAS EUR SAS; do
  out="${PCA_DIR}/labels.${SP}.keep"
  if [ ! -s "$out" ]; then
    awk -v sp="$SP" 'BEGIN{FS=OFS="\t"} NR>1 {gsub("\r",""); if ($3==sp) print 0,$1}' \
      "$LABELS_PATH" > "$out"
  fi
done

# Init output tmp files with header
for SP in AFR AMR EAS EUR SAS; do
  : > "${TMP_DIR}/${SP}.maf.tmp"        # clear if exists
  printf "%s\n" $'chr\tpos\tref\talt\tmaf' > "${TMP_DIR}/${SP}.maf.tmp"
done

# Per-chromosome: make PGEN (no PCA filters!) and compute freq
for c in {1..22}; do
  VCF_PATH="${VCF_PATTERN/\{chr\}/$c}"

  plink2 \
    --vcf "${VCF_PATH}" \
    --snps-only just-acgt --max-alleles 2 \
    --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 200 truncate \
    --make-pgen --out "${TMP_DIR}/maf_chr${c}"

  # For each superpop, compute allele freqs on that chr and append as chr/pos/ref/alt/maf
  for SP in AFR AMR EAS EUR SAS; do
    plink2 \
      --pfile "${TMP_DIR}/maf_chr${c}" \
      --keep "${PCA_DIR}/labels.${SP}.keep" \
      --freq \
      --out "${TMP_DIR}/freq_${SP}_chr${c}"

    # Convert .afreq to the 5-column format (robust to column order)
    awk -v OFS='\t' -f \
      "/app/scripts/helpers/awk/afreq_to_maf.awk" \
      "${TMP_DIR}/freq_${SP}_chr${c}.afreq" >> "${TMP_DIR}/${SP}.maf.tmp"
  done
done

for SP in AFR AMR EAS EUR SAS; do
  sort -u "${TMP_DIR}/${SP}.maf.tmp" > "${PCA_DIR}/maf_${SP}.grch38.tsv"
  gzip -f "${PCA_DIR}/maf_${SP}.grch38.tsv"
  echo "[OK] Wrote ${PCA_DIR}/maf_${SP}.grch38.tsv.gz"
done
