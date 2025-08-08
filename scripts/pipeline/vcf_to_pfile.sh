#!/usr/bin/env bash
# Convert an (imputed) VCF.gz to a PLINK2 PFILE trio (.pgen/.pvar/.psam),
# with optional filtering to well-imputed, biallelic SNPs.
#
# Usage:
#   vcf_to_pfile.sh --vcf FILE.vcf.gz --out-prefix /path/to/user \
#       [--info-key R2|INFO|IMPINFO] [--info-min 0.8] [--keep-unfiltered]
#
# Produces: <out-prefix>.pgen/.pvar/.psam
# Prints:   the out-prefix on stdout (so callers can capture it if desired)
set -euo pipefail

# Defaults
INFO_KEY=""
INFO_MIN="0.8"
KEEP_UNFILTERED="0"

VCF=""
OUT_PREFIX=""

die() { echo "✗ $*" >&2; exit 1; }

# --- arg parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="${2:-}"; shift 2;;
    --out-prefix) OUT_PREFIX="${2:-}"; shift 2;;
    --info-key) INFO_KEY="${2:-}"; shift 2;;
    --info-min) INFO_MIN="${2:-}"; shift 2;;
    --keep-unfiltered) KEEP_UNFILTERED="1"; shift;;
    -h|--help) grep -E '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

[[ -n "$VCF" && -n "$OUT_PREFIX" ]] || die "Required: --vcf FILE.vcf.gz and --out-prefix PREFIX"
[[ -f "$VCF" ]] || die "VCF not found: $VCF"
command -v plink2 >/dev/null || die "plink2 not found in PATH"
command -v bcftools >/dev/null || die "bcftools not found in PATH"
command -v tabix >/dev/null || die "tabix not found in PATH"

# Ensure indexed
[[ -f "${VCF}.tbi" ]] || tabix -f -p vcf "$VCF" || true

echo "==> Preparing PLINK2 PFILE from imputed VCF…"
echo "   - Input VCF: $VCF"
echo "   - Output prefix: $OUT_PREFIX"

# Detect dosage field (prefer DS, fallback GP)
VCF_HEADER="$(bcftools view -h "$VCF")"
if echo "$VCF_HEADER" | grep -q 'ID=DS,'; then
  DOSAGE_FIELD="DS"
elif echo "$VCF_HEADER" | grep -q 'ID=GP,'; then
  DOSAGE_FIELD="GP"
else
  die "Neither DS nor GP found in VCF FORMAT fields; cannot read dosages."
fi
echo "   - Using dosage field: $DOSAGE_FIELD"

# Choose imputation quality key if not supplied
if [[ -z "$INFO_KEY" ]]; then
  if echo "$VCF_HEADER" | grep -q 'INFO=<ID=R2'; then
    INFO_KEY="R2"
  elif echo "$VCF_HEADER" | grep -q 'INFO=<ID=INFO'; then
    INFO_KEY="INFO"
  elif echo "$VCF_HEADER" | grep -q 'INFO=<ID=IMPINFO'; then
    INFO_KEY="IMPINFO"
  else
    echo "   ! No obvious INFO tag for imputation quality; proceeding without filter unless --keep-unfiltered=0 and INFO_KEY provided."
    INFO_KEY=""
  fi
fi

# Filter (unless keep-unfiltered) to well-imputed, biallelic SNPs
VCF_FOR_PLINK="$VCF"
TMP_FILT=""
if [[ "$KEEP_UNFILTERED" -ne 1 ]]; then
  if [[ -n "$INFO_KEY" ]]; then
    echo "   - Filtering: biallelic SNPs with ${INFO_KEY} >= ${INFO_MIN}"
    TMP_FILT="$(dirname "$OUT_PREFIX")/$(basename "$OUT_PREFIX").filt.vcf.gz"
    bcftools view -i "INFO/${INFO_KEY}>=${INFO_MIN} && TYPE=\"snp\" && N_ALT=1" "$VCF" -Oz -o "$TMP_FILT"
    tabix -f -p vcf "$TMP_FILT" || true
    VCF_FOR_PLINK="$TMP_FILT"
  else
    echo "   ! INFO key unknown; applying SNP-only & biallelic restriction inside plink2, but skipping INFO-based filter."
  fi
else
  echo "   - Skipping filtering (--keep-unfiltered)"
fi

# Convert to PFILE
plink2 \
  --vcf "$VCF_FOR_PLINK" dosage="$DOSAGE_FIELD" \
  --double-id \
  --snps-only just-acgt \
  --max-alleles 2 \
  --new-id-max-allele-len 200 truncate \
  --make-pgen \
  --out "$OUT_PREFIX"

# Cleanup temp filtered file
if [[ -n "$TMP_FILT" && -f "$TMP_FILT" ]]; then
  rm -f "$TMP_FILT" "${TMP_FILT}.tbi" || true
fi

echo "==> Wrote: ${OUT_PREFIX}.pgen .pvar .psam"
# Print the prefix so caller can capture it if desired
echo "$OUT_PREFIX"
