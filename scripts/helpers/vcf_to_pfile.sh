#!/usr/bin/env bash
# Convert an imputed VCF.gz to a PLINK2 PFILE trio (.pgen/.pvar/.psam),

set -euo pipefail
set -E 

# -------- defaults --------
INFO_KEY=""
INFO_MIN="0.8"
THREADS="1"
KEEP_UNFILTERED="0"
EXTRACT_SNPS=""
NO_CLEANUP="0"

VCF=""
OUT_PREFIX=""

die() { echo "✗ $*" >&2; exit 1; }
log() { echo "==> $*"; }
warn() { echo "   ! $*" >&2; }

hash_file() {
  local f="$1"
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$f" | awk '{print $1}'
  elif command -v md5 >/dev/null 2>&1; then
    md5 -q "$f"
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$f" | awk '{print $1}'
  else
    echo "NA"
  fi
}

# Ensure VCF header declares INFO/END to avoid htslib warnings
add_end_header_if_missing() {
  local vcf="$1"
  [[ -f "$vcf" ]] || return 0
  if ! bcftools view -h "$vcf" | grep -q 'ID=END,'; then
    local tmp
    tmp="${vcf%.vcf.gz}.withEND.vcf.gz"
    HTS_LOG_LEVEL=ERROR bcftools annotate \
      -h <(printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n') \
      -Oz -o "$tmp" "$vcf"
    tabix -f -p vcf "$tmp" || true
    mv "$tmp" "$vcf" && mv "${tmp}.tbi" "${vcf}.tbi"
  fi
}

# -------- arg parsing --------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="${2:-}"; shift 2;;
    --out-prefix) OUT_PREFIX="${2:-}"; shift 2;;
    --info-key) INFO_KEY="${2:-}"; shift 2;;
    --info-min) INFO_MIN="${2:-}"; shift 2;;
    --threads) THREADS="${2:-}"; shift 2;;
    --extract-snps) EXTRACT_SNPS="${2:-}"; shift 2;;
    --keep-unfiltered) KEEP_UNFILTERED="1"; shift;;
    --no-cleanup) NO_CLEANUP="1"; shift;;
    -h|--help) grep -E '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

[[ -n "$VCF" && -n "$OUT_PREFIX" ]] || die "Required: --vcf FILE.vcf.gz and --out-prefix PREFIX"
[[ -f "$VCF" ]] || die "VCF not found: $VCF"
command -v plink2 >/dev/null || die "plink2 not found in PATH"
command -v bcftools >/dev/null || die "bcftools not found in PATH"
command -v tabix >/dev/null || die "tabix not found in PATH"
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be an integer"

# --- status handling for frontend ---
STATUS_JSON="${OUT_PREFIX}.status.json"
STAGE="startup"

json_string() {
  if command -v python3 >/dev/null 2>&1; then
    python3 - "$1" <<'PY'
import json,sys
print(json.dumps(sys.argv[1]))
PY
  else
    local s="${1//\\/\\\\}"; s="${s//\"/\\\"}"; s="${s//$'\n'/\\n}"
    printf '"%s"' "$s"
  fi
}

write_status() {
  local status="$1"; shift
  local message="${1:-}"
  {
    printf '{\n'
    printf '  "status": "%s",\n' "$status"
    printf '  "stage": "%s",\n' "$STAGE"
    printf '  "message": %s,\n' "$(json_string "$message")"
    printf '  "out_prefix": "%s",\n' "$OUT_PREFIX"
    printf '  "input_vcf": "%s"\n' "$(basename "$VCF")"
    printf '}\n'
  } > "$STATUS_JSON"
  echo "STATUS: ${status} ${STATUS_JSON}"
}

on_err() {
  local ec=$?
  write_status "failed" "cmd=${BASH_COMMAND}; line=${BASH_LINENO[0]}; exit=${ec}"
  exit "$ec"
}
trap on_err ERR
die() { write_status "failed" "$*"; echo "✗ $*" >&2; exit 1; }

# Temp workspace
TMPDIR_ROOT="${TMPDIR:-/tmp}"
WORKDIR="$(mktemp -d "${TMPDIR_ROOT%/}/vcf_to_pfile.XXXXXX")"
cleanup() {
  [[ "$NO_CLEANUP" = "1" ]] && { warn "Skipping cleanup (--no-cleanup). Workdir: $WORKDIR"; return 0; }
  rm -rf "$WORKDIR" || true
}
trap cleanup EXIT

# Ensure indexed
STAGE="index"
[[ -f "${VCF}.tbi" ]] || tabix -f -p vcf "$VCF" || true

# Add END header if missing before any bcftools/plink2 reads
add_end_header_if_missing "$VCF"

log "Preparing PLINK2 PFILE from imputed VCF…"
echo "   - Input VCF: $VCF"
echo "   - Output prefix: $OUT_PREFIX"
echo "   - Threads: $THREADS"

# Tool versions
BCFTOOLS_VER="$(bcftools --version | head -n1 || true)"
PLINK2_VER="$(plink2 --version 2>/dev/null | head -n1 || true)"
TABIX_VER="$(tabix --version 2>/dev/null | head -n1 || true)"

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

# ---- Detect imputation quality key (prefer DR2, then R2/AR2/RSQ/INFO/IMPINFO) ----
STAGE="detect_info"
if [[ -z "$INFO_KEY" ]]; then
  CANDIDATES=("DR2" "R2" "AR2" "RSQ" "Rsq" "INFO" "IMPINFO")
  for tag in "${CANDIDATES[@]}"; do
    if echo "$VCF_HEADER" | grep -q "INFO=<ID=${tag}\b"; then INFO_KEY="$tag"; break; fi
  done
  if [[ -z "$INFO_KEY" ]]; then
    for tag in "${CANDIDATES[@]}"; do
      if bcftools query -f "%INFO/${tag}\n" "$VCF" 2>/dev/null | head -n 10000 | grep -qv '^\.$'; then
        INFO_KEY="$tag"; break
      fi
    done
  fi
fi
[[ -n "$INFO_KEY" ]] && echo "   - INFO/imputation quality key: $INFO_KEY (min ${INFO_MIN})" || \
  warn "No obvious INFO tag for imputation quality; proceeding without INFO-based filter unless --info-key is provided."

# Count samples quickly
SAMPLE_CT="$(bcftools query -l "$VCF" | wc -l | awk '{print $1}')"

# ------------- filtering (NO MAF) -------------
VCF_FOR_PLINK="$VCF"
TMP_FILT=""

if [[ "$KEEP_UNFILTERED" -ne 1 ]]; then
  if [[ -n "$INFO_KEY" ]]; then
    STAGE="filter"
    log "Filtering: biallelic SNPs with ${INFO_KEY} >= ${INFO_MIN}"
    TMP_FILT="${WORKDIR}/$(basename "$OUT_PREFIX").filt.vcf.gz"
    bcftools view --threads "$THREADS" \
      -i "INFO/${INFO_KEY}>=${INFO_MIN} && TYPE=\"snp\" && N_ALT=1" \
      "$VCF" -Oz -o "$TMP_FILT"
    tabix -f -p vcf "$TMP_FILT" || true
    VCF_FOR_PLINK="$TMP_FILT"
  else
    warn "INFO key unknown; applying SNP-only & biallelic restriction inside plink2, skipping INFO filter."
  fi
else
  log "Skipping filtering (--keep-unfiltered)"
fi

# Ensure END header present on the exact file fed to plink2
add_end_header_if_missing "$VCF_FOR_PLINK"

# Double-check; in rare cases on some filesystems a rename may race with readers.
# If still missing, force-annotate to a new file and use that for plink2.
if ! bcftools view -h "$VCF_FOR_PLINK" | grep -q 'ID=END,'; then
  warn "'END' header still missing; forcing header injection for plink2 input."
  FOR_PLINK="${WORKDIR}/for_plink.vcf.gz"
  HTS_LOG_LEVEL=ERROR bcftools annotate \
    -h <(printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n') \
    -Oz -o "$FOR_PLINK" "$VCF_FOR_PLINK"
  tabix -f -p vcf "$FOR_PLINK" || true
  VCF_FOR_PLINK="$FOR_PLINK"
fi

# ------------- convert to PFILE -------------
STAGE="convert"
HTS_LOG_LEVEL=ERROR plink2 \
  --threads "$THREADS" \
  --vcf "$VCF_FOR_PLINK" dosage="$DOSAGE_FIELD" \
  --double-id \
  --snps-only just-acgt \
  --max-alleles 2 \
  --new-id-max-allele-len 200 truncate \
  --make-pgen \
  --out "$OUT_PREFIX"

# Optional: subset to a fixed SNP list (IDs)
SUBSET_PREFIX=""
SUBSET_SNPS_KEPT=0
if [[ -n "$EXTRACT_SNPS" ]]; then
  [[ -f "$EXTRACT_SNPS" ]] || die "--extract-snps file not found: $EXTRACT_SNPS"
  SUBSET_PREFIX="${OUT_PREFIX}.subset"
  STAGE="subset"
  log "Subsetting to SNP list: $EXTRACT_SNPS"
  plink2 \
    --threads "$THREADS" \
    --pfile "$OUT_PREFIX" \
    --extract "$EXTRACT_SNPS" \
    --make-pgen \
    --out "$SUBSET_PREFIX"
  SUBSET_SNPS_KEPT="$(awk 'BEGIN{c=0} !/^#/ {c++} END{print c}' "${SUBSET_PREFIX}.pvar")"
fi

# Panel size + fraction (only meaningful when --extract-snps used)
PANEL_SIZE=0
SUBSET_FRAC=""
if [[ -n "$EXTRACT_SNPS" ]]; then
  PANEL_SIZE=$(grep -cv '^[[:space:]]*$' "$EXTRACT_SNPS" || echo 0)
  if [[ "$PANEL_SIZE" -gt 0 ]]; then
    SUBSET_FRAC=$(awk -v k="$SUBSET_SNPS_KEPT" -v m="$PANEL_SIZE" 'BEGIN{printf "%.6f", (m>0? k/m : 0)}')
  else
    SUBSET_FRAC="0.0"
  fi
fi

# ------------- QC COUNTS + JSON REPORT -------------
STAGE="qc"
qc_json="${OUT_PREFIX}.qc.json"

vcount() {
  local vcf="$1"
  local expr="$2"
  bcftools view --threads "$THREADS" -H ${expr:+-i "$expr"} "$vcf" | wc -l | awk '{print $1}'
}

EXPR_SNP_BI='TYPE="snp" && N_ALT=1'
TOTAL_RECORDS=$(vcount "$VCF" "")
SNP_BIALLELIC=$(vcount "$VCF" "$EXPR_SNP_BI")

PASS_INFO_ONLY=0
if [[ -n "$INFO_KEY" ]]; then
  EXPR_INFO="INFO/${INFO_KEY}>=${INFO_MIN}"
  PASS_INFO_ONLY=$(vcount "$VCF" "$EXPR_SNP_BI && $EXPR_INFO")
fi

if [[ "$KEEP_UNFILTERED" -eq 0 && -n "$TMP_FILT" ]]; then
  KEPT_IN_FILTERED=$(vcount "$VCF_FOR_PLINK" "")
else
  KEPT_IN_FILTERED="$SNP_BIALLELIC"
fi

VCF_MD5="$(hash_file "$VCF")"
PVAR_LINES="$(awk 'BEGIN{c=0} !/^#/ {c++} END{print c}' "${OUT_PREFIX}.pvar")"
RUN_TS="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

{
  printf '{\n'
  printf '  "timestamp_utc": "%s",\n' "$RUN_TS"
  printf '  "tools": {\n'
  printf '    "bcftools": "%s",\n' "$BCFTOOLS_VER"
  printf '    "plink2": "%s",\n' "$PLINK2_VER"
  printf '    "tabix": "%s"\n' "$TABIX_VER"
  printf '  },\n'
  printf '  "input": {\n'
  printf '    "vcf": "%s",\n' "$(basename "$VCF")"
  printf '    "vcf_md5_or_sha": "%s",\n' "$VCF_MD5"
  printf '    "samples": %s\n' "$SAMPLE_CT"
  printf '  },\n'
  printf '  "params": {\n'
  printf '    "dosage_field": "%s",\n' "$DOSAGE_FIELD"
  printf '    "info_key": "%s",\n' "$INFO_KEY"
  printf '    "info_min": %s,\n' "$INFO_MIN"
  printf '    "keep_unfiltered": %s,\n' "$( [[ $KEEP_UNFILTERED -eq 1 ]] && echo true || echo false )"
  printf '    "threads": %s,\n' "$THREADS"
  printf '    "extract_snps": "%s"\n' "${EXTRACT_SNPS:-}"
  printf '  },\n'
  printf '  "counts": {\n'
  printf '    "total_records_vcf": %s,\n' "$TOTAL_RECORDS"
  printf '    "snp_biallelic_in_vcf": %s,\n' "$SNP_BIALLELIC"
  printf '    "pass_info_only_in_vcf": %s,\n' "$PASS_INFO_ONLY"
  printf '    "kept_variants_expected": %s,\n' "$KEPT_IN_FILTERED"
  printf '    "kept_variants_in_pvar": %s,\n' "$PVAR_LINES"
  if [[ -n "$SUBSET_PREFIX" ]]; then
    printf '    "subset_snps_kept": %s,\n' "$SUBSET_SNPS_KEPT"
    printf '    "panel_size": %s,\n' "$PANEL_SIZE"
    printf '    "subset_fraction": %s\n' "${SUBSET_FRAC:-0}"
  else
    printf '    "subset_snps_kept": null,\n'
    printf '    "panel_size": null,\n'
    printf '    "subset_fraction": null\n'
  fi
  printf '  }\n'
  printf '}\n'
} > "$qc_json"

# ---- final status evaluation for frontend ----
# Tightened gates live HERE (after QC JSON is written), using both absolute and relative thresholds.
STAGE="finalize"
STATUS="ok"
MSG=""

# Hard failure: nothing kept
if [[ "${PVAR_LINES}" -eq 0 || "${KEPT_IN_FILTERED}" -eq 0 ]]; then
  write_status "failed" "No variants kept after filters or conversion."
  exit 1
fi

# Bishop-2022-style panel (~300k pruned): gates
ABS_FAIL=30000   # <30k => fail
ABS_DEG=80000    # 30k-80k => degraded (unless fraction is OK)
REL_FAIL=0.05    # <5% of panel => fail
REL_DEG=0.20     # 5-20% => degraded

ISSUES=()

if [[ -n "$EXTRACT_SNPS" && "$PANEL_SIZE" -gt 0 ]]; then
  # Evaluate absolute gates
  if (( SUBSET_SNPS_KEPT < ABS_FAIL )); then
    STATUS="failed"
    MSG="Low panel overlap (kept=${SUBSET_SNPS_KEPT} < ${ABS_FAIL})."
  fi
  if [[ "$STATUS" != "failed" ]] && (( SUBSET_SNPS_KEPT < ABS_DEG )); then
    ISSUES+=("low_panel_overlap_abs kept=${SUBSET_SNPS_KEPT} thr=${ABS_DEG}")
  fi

  # Evaluate relative gates with awk for float comparison
  FAIL_REL=$(awk -v f="$SUBSET_FRAC" -v t="$REL_FAIL" 'BEGIN{print (f<t)?"1":"0"}')
  DEG_REL=$(awk -v f="$SUBSET_FRAC" -v t="$REL_DEG"  'BEGIN{print (f<t)?"1":"0"}')

  if [[ "$STATUS" != "failed" && "$FAIL_REL" == "1" ]]; then
    STATUS="failed"
    MSG="Low panel overlap fraction (frac=${SUBSET_FRAC} < ${REL_FAIL})."
  fi
  if [[ "$STATUS" != "failed" && "$DEG_REL" == "1" ]]; then
    ISSUES+=("low_panel_overlap_frac frac=${SUBSET_FRAC} thr=${REL_DEG}")
  fi
else
  # No --extract-snps: coarse sanity check only
  if [[ "$PVAR_LINES" -lt 100000 ]]; then
    ISSUES+=("kept_variants_in_pvar_below_100k")
  fi
fi

if [[ "$STATUS" != "failed" && "${#ISSUES[@]}" -gt 0 ]]; then
  STATUS="degraded"
  MSG="QC warnings: ${ISSUES[*]}"
fi

write_status "$STATUS" "$MSG"

# ------------- cleanup temp -------------
if [[ -n "$TMP_FILT" && -f "$TMP_FILT" ]]; then :; fi

log "Wrote: ${OUT_PREFIX}.pgen .pvar .psam"
[[ -n "$SUBSET_PREFIX" ]] && log "Wrote: ${SUBSET_PREFIX}.pgen .pvar .psam"
log "Wrote: ${qc_json}"
echo "$OUT_PREFIX"
