#!/usr/bin/env bash
set -Eeuo pipefail

# scripts/pipeline/build_refs.sh
# Purpose: Build local reference artifacts from already-downloaded data
#   - ref_bcfs_b38: Convert 1000G VCFs to BCF + CSI (bcftools)
#   - ref_brefs_b38: Convert 1000G VCFs to bref3 (bref3 JAR)

# -----------------------------------------------------------------------------
# Host vs container switch (mirror user.sh/add_pgs.sh behavior)
# -----------------------------------------------------------------------------
if [ -z "${INSIDE_DOCKER:-}" ]; then
  THIS_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
  PROJECT_ROOT=$(cd -- "${THIS_SCRIPT_DIR}/../../" &> /dev/null && pwd)
  DOCKER_IMAGE="${DOCKER_IMAGE:-pgspilot}"

  docker run --rm \
    -e INSIDE_DOCKER=1 \
    -v "${PROJECT_ROOT}:/app" \
    "$DOCKER_IMAGE" \
    /app/scripts/pipeline/build_refs.sh "$@"
  exit 0
fi

# -----------------------------------------------------------------------------
# Container worker
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/config.sh"

GENOME_DIR="${ROOT_DIR}/genome_data"

log(){ printf '%s\n' "$*" >&2; }
die(){ log "✗ $*"; exit 1; }
ensure_dir(){ mkdir -p "$1"; }

SUPPORTED=(
  "ref_bcfs_b38"
  "ref_brefs_b38"
)

list_supported(){ printf '%s\n' "${SUPPORTED[@]}"; }
contains(){ local n="$1"; shift || true; for x in "$@"; do [[ "$x" == "$n" ]] && return 0; done; return 1; }

setup_ref_bcfs_b38(){
  local out_dir="${GENOME_DIR}/ref_bcfs_b38"; ensure_dir "$out_dir"
  local in_dir="${GENOME_DIR}/1000G"
  local threads="${BCF_THREADS:-${THREADS:-4}}"

  log "==> [ref_bcfs_b38] Converting 1000G VCFs → BCF + CSI into $out_dir (parallel: $threads)"

  # Ensure inputs exist
  local missing=()
  local chr
  for chr in $(seq 1 22); do
    local in_vcf="${in_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    if [[ ! -s "$in_vcf" ]]; then
      missing+=("chr${chr}")
    fi
  done
  if (( ${#missing[@]} > 0 )); then
    die "[ref_bcfs_b38] Missing input VCFs for: ${missing[*]} (run download_data.sh --only 1000G first)"
  fi

  convert_one(){
    local chr="$1"
    local in_vcf="${in_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    local base
    base=$(basename "$in_vcf" .vcf.gz)
    local out_bcf="${out_dir}/${base}.bcf"
    local tmp_bcf="${out_bcf}.part"
    local out_csi="${out_bcf}.csi"

    if [[ -s "$out_bcf" && -s "$out_csi" && "$out_bcf" -nt "$in_vcf" && "$out_csi" -nt "$out_bcf" ]]; then
      log "[ref_bcfs_b38 chr${chr}] Up-to-date; skipping."
      return 0
    fi

    log "[ref_bcfs_b38 chr${chr}] Converting → BCF"
    if ! bcftools view -Ob -o "$tmp_bcf" --threads "$threads" "$in_vcf" >/dev/null 2>&1; then
      rm -f "$tmp_bcf" 2>/dev/null || true
      die "[ref_bcfs_b38 chr${chr}] bcftools view failed"
    fi
    if [[ ! -s "$tmp_bcf" ]]; then
      rm -f "$tmp_bcf" 2>/dev/null || true
      die "[ref_bcfs_b38 chr${chr}] Empty BCF after conversion"
    fi
    mv -f "$tmp_bcf" "$out_bcf"

    log "[ref_bcfs_b38 chr${chr}] Indexing (CSI)"
    if ! bcftools index -f -c "$out_bcf" >/dev/null 2>&1; then
      rm -f "$out_csi" 2>/dev/null || true
      die "[ref_bcfs_b38 chr${chr}] bcftools index failed"
    fi
    [[ -s "$out_csi" ]] || die "[ref_bcfs_b38 chr${chr}] Missing CSI index after indexing"

    log "[ref_bcfs_b38 chr${chr}] Done: $(basename "$out_bcf"), $(basename "$out_csi")"
  }

  local -a pids=()
  local -i max_jobs=$(( (threads>0)?threads:1 ))
  for chr in $(seq 1 22); do
    convert_one "$chr" &
    pids+=("$!")
    if (( ${#pids[@]} >= max_jobs )); then
      wait "${pids[0]}" || die "One of the BCF conversions failed."
      unset 'pids[0]'
      pids=(${pids[@]})
    fi
  done
  local pid; for pid in "${pids[@]}"; do wait "$pid" || die "One of the BCF conversions failed."; done
  log "[ref_bcfs_b38] All chromosomes converted and indexed."
}

setup_ref_brefs_b38(){
  local out_dir="${GENOME_DIR}/ref_brefs_b38"; ensure_dir "$out_dir"
  local in_dir="${GENOME_DIR}/1000G"
  local jar_dir="${GENOME_DIR}/jars"
  local bref3_jar="${jar_dir}/bref3.27Feb25.75f.jar"
  local java_mem="${BREF3_MEM:-4g}"
  local threads="${BREF3_THREADS:-${THREADS:-4}}"

  log "==> [ref_brefs_b38] Converting 1000G VCFs → bref3 into $out_dir (parallel: $threads)"
  [[ -s "$bref3_jar" ]] || die "[ref_brefs_b38] Missing $bref3_jar (run download_data.sh --only jars first)"

  local missing=()
  local chr
  for chr in $(seq 1 22); do
    local in_vcf="${in_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    [[ -s "$in_vcf" ]] || missing+=("chr${chr}")
  done
  (( ${#missing[@]} == 0 )) || die "[ref_brefs_b38] Missing input VCFs for: ${missing[*]} (run download_data.sh --only 1000G first)"

  convert_one(){
    local chr="$1"
    local in_vcf="${in_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    local base; base=$(basename "$in_vcf" .vcf.gz)
    local out_prefix="${out_dir}/${base}"; local out_bref3="${out_prefix}.bref3"

    if [[ -s "$out_bref3" && "$out_bref3" -nt "$in_vcf" ]]; then
      log "[ref_brefs_b38 chr${chr}] Up-to-date; skipping."
      return 0
    fi

    local tmpdir; tmpdir=$(mktemp -d)
    local cleanup_tmp(){ rm -rf "$tmpdir" 2>/dev/null || true; }; trap cleanup_tmp RETURN
    local tmp_prefix="${tmpdir}/chr${chr}"; local tmp_bref3="${tmp_prefix}.bref3"

    log "[ref_brefs_b38 chr${chr}] Running bref3 converter"
    if ! java -Xmx"$java_mem" -jar "$bref3_jar" vcf="$in_vcf" out="$tmp_prefix" >/dev/null 2>&1; then
      log "[ref_brefs_b38 chr${chr}] First attempt failed; retrying once…"; sleep 3
      java -Xmx"$java_mem" -jar "$bref3_jar" vcf="$in_vcf" out="$tmp_prefix" >/dev/null 2>&1 || die "[ref_brefs_b38 chr${chr}] bref3 conversion failed"
    fi

    [[ -s "$tmp_bref3" ]] || die "[ref_brefs_b38 chr${chr}] Missing output bref3 from converter"
    mv -f "$tmp_bref3" "$out_bref3"
    log "[ref_brefs_b38 chr${chr}] Done: $(basename "$out_bref3")"
  }

  local -a pids=(); local -i max_jobs=$(( (threads>0)?threads:1 ))
  for chr in $(seq 1 22); do
    convert_one "$chr" &
    pids+=("$!")
    if (( ${#pids[@]} >= max_jobs )); then
      wait "${pids[0]}" || die "One of the bref3 conversions failed."
      unset 'pids[0]'; pids=(${pids[@]})
    fi
  done
  local pid; for pid in "${pids[@]}"; do wait "$pid" || die "One of the bref3 conversions failed."; done
  log "[ref_brefs_b38] All chromosomes converted to bref3."
}

run_selected(){
  local -a selected=("$@")
  local name
  ensure_dir "$GENOME_DIR"
  for name in "${selected[@]}"; do
    case "$name" in
      ref_bcfs_b38)   setup_ref_bcfs_b38 ;;
      ref_brefs_b38)  setup_ref_brefs_b38 ;;
      *) die "Unknown target: $name (use --list)" ;;
    esac
  done
}

ONLY=""; DO_LIST=0; DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY="$2"; shift 2;;
    --list) DO_LIST=1; shift;;
    --dry-run) DRY_RUN=1; shift;;
    --all|*) shift || true;;
  esac
done

if (( DO_LIST )); then list_supported; exit 0; fi
SELECTED=()
if [[ -n "$ONLY" ]]; then IFS=',' read -r -a SELECTED <<< "$ONLY"; else SELECTED=("${SUPPORTED[@]}"); fi
for n in "${SELECTED[@]}"; do contains "$n" "${SUPPORTED[@]}" || die "Unsupported target: $n"; done

log "==> genome_data root: $GENOME_DIR"
if (( DRY_RUN )); then log "[dry-run] Would build: ${SELECTED[*]}"; exit 0; fi
run_selected "${SELECTED[@]}"
log "✓ Reference builds complete."
