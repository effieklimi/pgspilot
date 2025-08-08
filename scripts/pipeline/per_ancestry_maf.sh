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
# --------- Build per-superpop keep lists (EUR/AFR/EAS/AMR/SAS) ----------
for SP in AFR AMR EAS EUR SAS; do
  awk -v sp="$SP" 'BEGIN{FS=OFS="\t"} NR>1 {gsub("\r",""); if ($3==sp) print 0,$1}' \
    "$LABELS_PATH" > "${PCA_DIR}/labels.${SP}.keep"
done

# --------- Init output tmp files with header ----------
for SP in AFR AMR EAS EUR SAS; do
  : > "${TMP_DIR}/${SP}.maf.tmp"        # clear if exists
  echo -e "chr\tpos\tref\talt\tmaf" > "${TMP_DIR}/${SP}.maf.tmp"
done

# --------- Per-chromosome: make PGEN (no PCA filters!) and compute freq ----------
for c in {1..22}; do
  VCF_PATH="${VCF_PATTERN/\{chr\}/$c}"

  # Create inclusive per-chr PGEN for frequency calc
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

    # Convert .afreq to the 5-column format, robust to column order
    awk -v OFS="\t" '
      NR==1{
        for(i=1;i<=NF;i++){
          if($i=="#CHROM"||$i=="CHROM") C=i;
          else if($i=="POS") P=i;
          else if($i=="REF") R=i;
          else if($i=="ALT") A=i;
          else if($i=="MAF") M=i;
          else if($i=="A1_FREQ") F=i
        }
        next
      }
      {
        chr=$C; if (chr!~/^chr/) chr="chr"chr;
        maf=(M? $M : (($F>0.5)?1-$F:$F));
        print chr, $P, $R, $A, maf
      }' "${TMP_DIR}/freq_${SP}_chr${c}.afreq" >> "${TMP_DIR}/${SP}.maf.tmp"
  done
done

# --------- Finalise (sort/uniq, gzip) ----------
for SP in AFR AMR EAS EUR SAS; do
  sort -u "${TMP_DIR}/${SP}.maf.tmp" > "${PCA_DIR}/maf_${SP}.grch38.tsv"
  gzip -f "${PCA_DIR}/maf_${SP}.grch38.tsv"
  echo "[OK] Wrote ${PCA_DIR}/maf_${SP}.grch38.tsv.gz"
done
