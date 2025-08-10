#!/usr/bin/env bash
set -Eeuo pipefail

# scripts/pipeline/download_data.sh
# Purpose: Orchestrate downloading and setting up genome_data subfolders.
# NOTE: Implements robust download for the 1000G subfolder; others are placeholders.

# -----------------------------------------------------------------------------
# Host vs container switch (mirror user.sh/add_pgs.sh behavior)
# -----------------------------------------------------------------------------
if [ -z "${INSIDE_DOCKER:-}" ]; then
  # Resolve repo root and re-invoke inside Docker with the repo mounted at /app
  THIS_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
  PROJECT_ROOT=$(cd -- "${THIS_SCRIPT_DIR}/../../" &> /dev/null && pwd)
  DOCKER_IMAGE="${DOCKER_IMAGE:-pgspilot}"

  docker run --rm \
    -e INSIDE_DOCKER=1 \
    -v "${PROJECT_ROOT}:/app" \
    "$DOCKER_IMAGE" \
    /app/scripts/pipeline/download_data.sh "$@"
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

ensure_dir(){
  local dir="$1"
  mkdir -p "$dir"
}

# Supported subfolders detected in current tree
SUPPORTED=(
  "1000G"
  "beagle_maps_b38"
  "chain"
  "eagle_maps_b38"
  "fasta"
  "jars"
  "ref_bcfs_b38"
  "ref_brefs_b38"
)

list_supported(){
  printf '%s\n' "${SUPPORTED[@]}"
}

contains() {
  local needle="$1"; shift || true
  local item
  for item in "$@"; do [[ "$item" == "$needle" ]] && return 0; done
  return 1
}

# Download logic for 1000G (chromosomes 1..22)
setup_1000G(){
  local dir="${GENOME_DIR}/1000G"; ensure_dir "$dir"
  local base_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
  local threads="${DOWNLOAD_THREADS:-${THREADS:-4}}"

  log "==> [1000G] Downloading chr1..22 VCFs into $dir (parallel: $threads)"

  # Fetch remote Content-Length for a URL (best-effort)
  get_remote_size(){
    local url="$1"
    # wget --spider emits headers on stderr; parse Length or Content-Length
    local size
    size=$(wget --spider --server-response -O - "$url" 2>&1 \
      | awk 'tolower($0) ~ /content-length:/ {gsub("\r","",$2); print $2; exit} /Length: [0-9]+/ {gsub("\r","",$2); print $2; exit}') || true
    printf '%s' "${size:-}"
  }

  validate_header(){
    local file="$1"
    # Ensure compressed and contains VCF header lines near start
    if zgrep -a -m1 '^##fileformat=VCF' "$file" >/dev/null 2>&1; then
      return 0
    fi
    if zgrep -a -m1 '^#CHROM' "$file" >/dev/null 2>&1; then
      return 0
    fi
    return 1
  }

  download_chr(){
    local chr="$1"
    local name="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    local url="${base_url}/${name}"
    local dest="${dir}/${name}"
    local tmp="${dest}.part"
    local name_tbi="${name}.tbi"
    local url_tbi="${base_url}/${name_tbi}"
    local dest_tbi="${dir}/${name_tbi}"
    local tmp_tbi="${dest_tbi}.part"

    # Determine remote size if possible
    local remote_size
    remote_size="$(get_remote_size "$url")"

    # If final file exists and matches expected size (when known), skip
    if [[ -s "$dest" ]]; then
      if [[ -n "$remote_size" ]]; then
        local local_size
        local_size=$(stat -c '%s' "$dest" 2>/dev/null || stat -f '%z' "$dest")
        if [[ "$local_size" == "$remote_size" ]]; then
          log "[1000G chr${chr}] Exists (size OK: $local_size). Skipping."
          return 0
        fi
      fi
      # Fallback: light header validation
      if validate_header "$dest"; then
        log "[1000G chr${chr}] Exists (header OK). Skipping."
        return 0
      else
        log "[1000G chr${chr}] Existing file seems incomplete; will resume."
        mv -f "$dest" "$tmp"
      fi
    fi

    # If tmp exists but final does not, resume into tmp
    if [[ -f "$tmp" && -n "$remote_size" ]]; then
      local part_size
      part_size=$(stat -c '%s' "$tmp" 2>/dev/null || stat -f '%z' "$tmp")
      log "[1000G chr${chr}] Resuming: $part_size / ${remote_size:-?} bytes"
    fi

    # Robust wget with resume and retries
    # Options: -c resume, --tries, --waitretry, --read-timeout, --timeout, --retry-connrefused
    if ! wget -c \
          --tries=10 \
          --waitretry=5 \
          --read-timeout=60 \
          --timeout=60 \
          --retry-connrefused \
          --no-verbose \
          -O "$tmp" "$url"; then
      log "[1000G chr${chr}] Download failed. Retrying once after 10s…"
      sleep 10
      wget -c --tries=10 --waitretry=5 --read-timeout=60 --timeout=60 --retry-connrefused --no-verbose -O "$tmp" "$url"
    fi

    # Verify completion (VCF)
    if [[ -n "$remote_size" ]]; then
      local final_size
      final_size=$(stat -c '%s' "$tmp" 2>/dev/null || stat -f '%z' "$tmp")
      if [[ "$final_size" != "$remote_size" ]]; then
        die "[1000G chr${chr}] Size mismatch after download ($final_size != $remote_size)"
      fi
    fi

    if ! validate_header "$tmp"; then
      die "[1000G chr${chr}] Missing VCF header signature; file may be corrupted."
    fi

    mv -f "$tmp" "$dest"
    log "[1000G chr${chr}] Done → $(basename "$dest")"

    # --- Download matching TBI index ---
    # Determine remote size for TBI if possible
    local remote_size_tbi
    remote_size_tbi="$(get_remote_size "$url_tbi")"

    # Skip if .tbi exists and appears complete
    if [[ -s "$dest_tbi" ]]; then
      if [[ -n "$remote_size_tbi" ]]; then
        local local_size_tbi
        local_size_tbi=$(stat -c '%s' "$dest_tbi" 2>/dev/null || stat -f '%z' "$dest_tbi")
        if [[ "$local_size_tbi" == "$remote_size_tbi" ]]; then
          log "[1000G chr${chr}] TBI exists (size OK: $local_size_tbi). Skipping."
          return 0
        fi
      fi
      if [[ -s "$dest_tbi" ]]; then
        log "[1000G chr${chr}] TBI exists (non-empty). Skipping."
        return 0
      fi
      mv -f "$dest_tbi" "$tmp_tbi"
    fi

    # Resume TBI if partial exists
    if [[ -f "$tmp_tbi" && -n "$remote_size_tbi" ]]; then
      local part_size_tbi
      part_size_tbi=$(stat -c '%s' "$tmp_tbi" 2>/dev/null || stat -f '%z' "$tmp_tbi")
      log "[1000G chr${chr}] TBI resume: $part_size_tbi / ${remote_size_tbi:-?} bytes"
    fi

    if ! wget -c \
          --tries=10 \
          --waitretry=5 \
          --read-timeout=60 \
          --timeout=60 \
          --retry-connrefused \
          --no-verbose \
          -O "$tmp_tbi" "$url_tbi"; then
      log "[1000G chr${chr}] TBI download failed. Retrying once after 10s…"
      sleep 10
      wget -c --tries=10 --waitretry=5 --read-timeout=60 --timeout=60 --retry-connrefused --no-verbose -O "$tmp_tbi" "$url_tbi"
    fi

    if [[ -n "$remote_size_tbi" ]]; then
      local final_size_tbi
      final_size_tbi=$(stat -c '%s' "$tmp_tbi" 2>/dev/null || stat -f '%z' "$tmp_tbi")
      if [[ "$final_size_tbi" != "$remote_size_tbi" ]]; then
        die "[1000G chr${chr}] TBI size mismatch after download ($final_size_tbi != $remote_size_tbi)"
      fi
    fi

    if [[ ! -s "$tmp_tbi" ]]; then
      die "[1000G chr${chr}] TBI file is empty after download."
    fi

    mv -f "$tmp_tbi" "$dest_tbi"
    log "[1000G chr${chr}] TBI Done → $(basename "$dest_tbi")"
  }

  # Concurrency controller
  local -a pids=()
  local -i max_jobs
  max_jobs=$(( threads > 0 ? threads : 1 ))

  local chr
  for chr in $(seq 1 22); do
    download_chr "$chr" &
    pids+=("$!")
    if (( ${#pids[@]} >= max_jobs )); then
      wait "${pids[0]}" || die "One of the downloads failed."
      unset 'pids[0]'
      # compact array
      pids=(${pids[@]})
    fi
  done
  # wait remaining
  local pid
  for pid in "${pids[@]}"; do
    wait "$pid" || die "One of the downloads failed."
  done

  log "[1000G] All chromosomes completed."
}

setup_beagle_maps_b38(){
  local dir="${GENOME_DIR}/beagle_maps_b38"; ensure_dir "$dir"
  local url_zip="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"
  local dest_zip="${dir}/plink.GRCh38.map.zip"
  local tmp_zip="${dest_zip}.part"

  # If sufficient .map files already exist, skip
  local existing_maps
  existing_maps=$(find "$dir" -maxdepth 1 -type f -name '*.map' | wc -l | tr -d ' ')
  if [[ "${existing_maps:-0}" -ge 40 ]]; then
    log "==> [beagle_maps_b38] Found ${existing_maps} .map files; skipping download."
    return 0
  fi

  log "==> [beagle_maps_b38] Downloading Beagle GRCh38 PLINK maps zip…"

  # Get remote size if available
  local remote_size
  remote_size="$(wget --spider --server-response -O - "$url_zip" 2>&1 \
    | awk 'tolower($0) ~ /content-length:/ {gsub("\r","",$2); print $2; exit} /Length: [0-9]+/ {gsub("\r","",$2); print $2; exit}')" || true

  # If completed zip exists with same size, skip downloading
  if [[ -s "$dest_zip" && -n "$remote_size" ]]; then
    local local_size
    local_size=$(stat -c '%s' "$dest_zip" 2>/dev/null || stat -f '%z' "$dest_zip")
    if [[ "$local_size" == "$remote_size" ]]; then
      log "[beagle_maps_b38] Zip exists (size OK: $local_size). Reusing."
    else
      log "[beagle_maps_b38] Existing zip size mismatch; will resume."
      mv -f "$dest_zip" "$tmp_zip"
    fi
  fi

  # Download/resume zip
  if [[ ! -s "$dest_zip" ]]; then
    if ! wget -c \
          --tries=10 \
          --waitretry=5 \
          --read-timeout=60 \
          --timeout=60 \
          --retry-connrefused \
          --no-verbose \
          -O "$tmp_zip" "$url_zip"; then
      log "[beagle_maps_b38] Zip download failed. Retrying once after 10s…"
      sleep 10
      wget -c --tries=10 --waitretry=5 --read-timeout=60 --timeout=60 --retry-connrefused --no-verbose -O "$tmp_zip" "$url_zip"
    fi

    # Verify size if known
    if [[ -n "$remote_size" ]]; then
      local final_size
      final_size=$(stat -c '%s' "$tmp_zip" 2>/dev/null || stat -f '%z' "$tmp_zip")
      if [[ "$final_size" != "$remote_size" ]]; then
        die "[beagle_maps_b38] Zip size mismatch after download ($final_size != $remote_size)"
      fi
    fi

    mv -f "$tmp_zip" "$dest_zip"
  fi

  # Test and extract zip to a temporary directory
  local tmpdir
  tmpdir=$(mktemp -d)
  # Ensure cleanup on function exit
  local cleanup_tmp
  cleanup_tmp(){ rm -rf "$tmpdir" 2>/dev/null || true; }
  trap cleanup_tmp RETURN

  if ! unzip -tq "$dest_zip" >/dev/null; then
    die "[beagle_maps_b38] Zip integrity test failed."
  fi

  unzip -q -o "$dest_zip" -d "$tmpdir"

  # Move .map files into destination directory
  shopt -s nullglob
  local moved=0
  local f
  while IFS= read -r -d '' f; do
    mv -f "$f" "$dir/"
    moved=$((moved+1))
  done < <(find "$tmpdir" -type f -name '*.map' -print0)
  shopt -u nullglob

  if [[ "$moved" -eq 0 ]]; then
    die "[beagle_maps_b38] No .map files found in the archive."
  fi

  log "[beagle_maps_b38] Moved ${moved} .map files to $dir"

  # Optionally keep a readme or manifest if present
  if [[ -f "$tmpdir/README.txt" ]]; then
    mv -f "$tmpdir/README.txt" "$dir/" || true
  fi

  # Optionally remove the zip after successful extraction to save space
  rm -f "$dest_zip" 2>/dev/null || true
}

setup_chain(){
  local dir="${GENOME_DIR}/chain"; ensure_dir "$dir"
  log "==> [chain] Placeholder: liftOver chain files (e.g., hg19ToHg38.over.chain.gz) in $dir"
}

setup_eagle_maps_b38(){
  local dir="${GENOME_DIR}/eagle_maps_b38"; ensure_dir "$dir"
  log "==> [eagle_maps_b38] Placeholder: recombination maps for Eagle (gz + tbi) in $dir"
}

setup_fasta(){
  local dir="${GENOME_DIR}/fasta"; ensure_dir "$dir"
  log "==> [fasta] Placeholder: reference FASTA(s) and .fai index in $dir"
  log "    Expected: Homo_sapiens_assembly19.fasta(.fai), Homo_sapiens_assembly38.fasta(.fai)"
}

setup_jars(){
  local dir="${GENOME_DIR}/jars"; ensure_dir "$dir"
  log "==> [jars] Placeholder: required Java JARs (e.g., picard/htsjdk/others) in $dir"
}

setup_ref_bcfs_b38(){
  local dir="${GENOME_DIR}/ref_bcfs_b38"; ensure_dir "$dir"
  log "==> [ref_bcfs_b38] Placeholder: reference BCFs (.bcf/.csi) by chromosome in $dir"
}

setup_ref_brefs_b38(){
  local dir="${GENOME_DIR}/ref_brefs_b38"; ensure_dir "$dir"
  log "==> [ref_brefs_b38] Placeholder: reference bref3 files by chromosome in $dir"
}

run_selected(){
  local -a selected=("$@")
  local name

  ensure_dir "$GENOME_DIR"

  for name in "${selected[@]}"; do
    case "$name" in
      1000G)             setup_1000G ;;
      beagle_maps_b38)   setup_beagle_maps_b38 ;;
      bin)               setup_bin ;;
      chain)             setup_chain ;;
      eagle_maps_b38)    setup_eagle_maps_b38 ;;
      fasta)             setup_fasta ;;
      jars)              setup_jars ;;
      ref_bcfs_b38)      setup_ref_bcfs_b38 ;;
      ref_brefs_b38)     setup_ref_brefs_b38 ;;
      *) die "Unknown subfolder: $name (use --list)" ;;
    esac
  done
}

# CLI
ONLY=""
DO_LIST=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --only)
      ONLY="$2"; shift 2 ;;
    --list)
      DO_LIST=1; shift ;;
    --dry-run)
      DRY_RUN=1; shift ;;
    --all|*)
      # default is --all; unrecognized positional triggers default path
      shift || true ;;
  esac
done

if (( DO_LIST )); then
  list_supported
  exit 0
fi

# Build selection
SELECTED=()
if [[ -n "$ONLY" ]]; then
  IFS=',' read -r -a SELECTED <<< "$ONLY"
else
  SELECTED=("${SUPPORTED[@]}")
fi

# Validate names
for n in "${SELECTED[@]}"; do
  contains "$n" "${SUPPORTED[@]}" || die "Unsupported subfolder: $n"
done

log "==> genome_data root: $GENOME_DIR"

if (( DRY_RUN )); then
  log "[dry-run] Would set up the following subfolders: ${SELECTED[*]}"
  exit 0
fi

run_selected "${SELECTED[@]}"

log "✓ Placeholders executed. Implement download logic per subfolder next."
