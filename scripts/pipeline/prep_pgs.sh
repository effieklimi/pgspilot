#!/usr/bin/env bash
set -Eeuo pipefail

# prep_pgs.sh
# Prepare new/updated PGS weights end-to-end: harmonize → include lists → calibration.
# Usage:
#   ./prep_pgs.sh weights_manifest.csv
# Manifest CSV columns (header required):
#   path,trait,ancestry
# Example:
#   weights_raw/HeartRate/PGS000300.txt,HeartRate,EUR
#   weights_raw/Height/PGS002123.txt,Height,EUR
#   weights_raw/Height/PGS002456.txt,Height,AFR

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 weights_manifest.csv" >&2
  exit 1
fi

MANIFEST="$1"
[[ -f "$MANIFEST" ]] || { echo "Manifest not found: $MANIFEST" >&2; exit 1; }

# Load config
source "./config.sh"

mkdir -p "$WEIGHTS_HM_DIR" "$TRAITS_DIR" "$CALIB_DIR"

echo "==> [1/3] Harmonizing weights to GRCh38 (per manifest)…"
# Skip header (NR>1). CSV: path,trait,ancestry
# handle commas inside paths by using awk with FPAT; simpler assumption: no commas in path
awk -F',' 'NR>1 {print $1"\t"$2"\t"$3}' "$MANIFEST" | while IFS=$'\t' read -r SRC TRAIT ANC; do
  [[ -f "$SRC" ]] || { echo "  ✗ Missing source: $SRC" >&2; exit 1; }
  echo "  [+] $TRAIT / $ANC  <=  $SRC"
  $PYTHON "$WEIGHTS_HARMONIZER" \
    --in "$SRC" \
    --out "$WEIGHTS_HM_DIR" \
    --fasta38 "$FASTA38" \
    --chain37to38 "$CHAIN_37_TO_38" \
    --ancestry "$ANC" \
    $DROP_MHC
done

echo "==> [2/3] Building frozen include lists…"
SITE_MASK_ARGS=()
[[ -n "$SITE_MASK" ]] && SITE_MASK_ARGS=(--site-mask "$SITE_MASK" --mask-min-dr2 "$MASK_MIN_DR2")
$PYTHON "$INCLUDE_BUILDER" \
  --weights "$WEIGHTS_HM_DIR" \
  --scorable "$SCORABLE_SITES" \
  --include-chroms "$INCLUDE_CHROMS" \
  $DROP_MHC \
  "${SITE_MASK_ARGS[@]}" \
  --out "$TRAITS_DIR"

echo "==> [3/3] Building calibration from 1000G…"
$PYTHON "$CALIB_BUILDER" \
  --weights_dir "$WEIGHTS_HM_DIR" \
  --traits_dir "$TRAITS_DIR" \
  --vcf_pattern "$ONEKG_VCF_PATTERN" \
  --labels "$ONEKG_LABELS" \
  --out "$CALIB_DIR" \
  --quantiles "$QUANTILES" \
  --min-variants "$MIN_VARIANTS"

echo "✓ Prep complete."
echo "   Harmonized weights: $WEIGHTS_HM_DIR"
echo "   Include lists:      $TRAITS_DIR"
echo "   Calibration:        $CALIB_DIR"
