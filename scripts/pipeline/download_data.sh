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
  local url="https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz"
  local dest="${dir}/hg19ToHg38.over.chain.gz"
  local tmp="${dest}.part"

  log "==> [chain] Ensuring liftOver chain at $dest"

  # Best-effort size probe
  local remote_size
  remote_size="$(wget --spider --server-response -O - "$url" 2>&1 \
    | awk 'tolower($0) ~ /content-length:/ {gsub("\r","",$2); print $2; exit} /Length: [0-9]+/ {gsub("\r","",$2); print $2; exit}')" || true

  # If exists and passes gzip -t (and size matches if known), skip
  if [[ -s "$dest" ]]; then
    if gzip -t "$dest" 2>/dev/null; then
      if [[ -n "$remote_size" ]]; then
        local local_size
        local_size=$(stat -c '%s' "$dest" 2>/dev/null || stat -f '%z' "$dest")
        if [[ "$local_size" == "$remote_size" ]]; then
          log "[chain] Existing file valid and size matches ($local_size). Skipping."
          return 0
      
        fi
      else
        log "[chain] Existing file valid (gzip test passed). Skipping."
        return 0
      fi
    else
      log "[chain] Existing file failed gzip test; will resume/overwrite."
      mv -f "$dest" "$tmp"
    fi
  fi

  # Resume or fresh download
  if ! wget -c \
        --tries=10 \
        --waitretry=5 \
        --read-timeout=60 \
        --timeout=60 \
        --retry-connrefused \
        --no-verbose \
        -O "$tmp" "$url"; then
    log "[chain] Download failed. Retrying once after 10s…"
    sleep 10
    wget -c --tries=10 --waitretry=5 --read-timeout=60 --timeout=60 --retry-connrefused --no-verbose -O "$tmp" "$url"
  fi

  # Verify size if known and gzip integrity
  if [[ -n "$remote_size" ]]; then
    local final_size
    final_size=$(stat -c '%s' "$tmp" 2>/dev/null || stat -f '%z' "$tmp")
    if [[ "$final_size" != "$remote_size" ]]; then
      die "[chain] Size mismatch after download ($final_size != $remote_size)"
    fi
  fi

  if ! gzip -t "$tmp" 2>/dev/null; then
    die "[chain] Gzip integrity check failed after download."
  fi

  mv -f "$tmp" "$dest"
  log "[chain] Ready: $(basename "$dest")"
}

setup_eagle_maps_b38(){
  local dir="${GENOME_DIR}/eagle_maps_b38"; ensure_dir "$dir"
  local url="https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
  local dest="${dir}/genetic_map_hg38_withX.txt.gz"
  local tmp="${dest}.part"

  log "==> [eagle_maps_b38] Ensuring Eagle genetic map at $dest"

  # Best-effort size probe
  local remote_size
  remote_size="$(wget --spider --server-response -O - "$url" 2>&1 \
    | awk 'tolower($0) ~ /content-length:/ {gsub("\r","",$2); print $2; exit} /Length: [0-9]+/ {gsub("\r","",$2); print $2; exit}')" || true

  # If exists and passes gzip -t (and size matches if known), skip
  if [[ -s "$dest" ]]; then
    if gzip -t "$dest" 2>/dev/null; then
      if [[ -n "$remote_size" ]]; then
        local local_size
        local_size=$(stat -c '%s' "$dest" 2>/dev/null || stat -f '%z' "$dest")
        if [[ "$local_size" == "$remote_size" ]]; then
          log "[eagle_maps_b38] Existing file valid and size matches ($local_size). Skipping."
          return 0
        fi
      else
        log "[eagle_maps_b38] Existing file valid (gzip test). Skipping."
        return 0
      fi
    else
      log "[eagle_maps_b38] Existing file failed gzip test; will resume/overwrite."
      mv -f "$dest" "$tmp"
    fi
  fi

  # Resume or fresh download
  if ! wget -c \
        --tries=10 \
        --waitretry=5 \
        --read-timeout=60 \
        --timeout=60 \
        --retry-connrefused \
        --no-verbose \
        -O "$tmp" "$url"; then
    log "[eagle_maps_b38] Download failed. Retrying once after 10s…"
    sleep 10
    wget -c --tries=10 --waitretry=5 --read-timeout=60 --timeout=60 --retry-connrefused --no-verbose -O "$tmp" "$url"
  fi

  # Verify size if known and gzip integrity
  if [[ -n "$remote_size" ]]; then
    local final_size
    final_size=$(stat -c '%s' "$tmp" 2>/dev/null || stat -f '%z' "$tmp")
    if [[ "$final_size" != "$remote_size" ]]; then
      die "[eagle_maps_b38] Size mismatch after download ($final_size != $remote_size)"
    fi
  fi

  if ! gzip -t "$tmp" 2>/dev/null; then
    die "[eagle_maps_b38] Gzip integrity check failed after download."
  fi

  mv -f "$tmp" "$dest"
  log "[eagle_maps_b38] Ready: $(basename "$dest")"
}

setup_fasta(){
  local dir="${GENOME_DIR}/fasta"; ensure_dir "$dir"
  log "==> [fasta] Ensuring GRCh38 and GRCh19 FASTA + .fai in $dir"

  # Helper: download .gz with resume, verify gzip, decompress atomically to .fasta
  download_and_unpack(){
    local url="$1"      # source .gz URL
    local dest_fasta="$2"  # final .fasta path

    local dest_gz="${dest_fasta}.gz"
    local tmp_gz="${dest_gz}.part"
    local tmp_fa="${dest_fasta}.part"

    # If FASTA exists and looks valid, just ensure index
    if [[ -s "$dest_fasta" ]]; then
      if head -n1 "$dest_fasta" | grep -q '^>' ; then
        :
      else
        log "[fasta] Existing $(basename "$dest_fasta") failed header check; will re-download."
        rm -f "$dest_fasta"
      fi
    fi

    if [[ ! -s "$dest_fasta" ]]; then
      # Download with resume and retries (follow redirects)
      if ! wget -c -L \
            --tries=10 \
            --waitretry=5 \
            --read-timeout=60 \
            --timeout=60 \
            --retry-connrefused \
            --no-verbose \
            -O "$tmp_gz" "$url"; then
        log "[fasta] Download failed for $(basename "$dest_fasta"). Retrying once after 10s…"
        sleep 10
        wget -c -L --tries=10 --waitretry=5 --read-timeout=60 --timeout=60 --retry-connrefused --no-verbose -O "$tmp_gz" "$url"
      fi

      # Validate and decompress
      if ! gzip -t "$tmp_gz" 2>/dev/null; then
        die "[fasta] Gzip integrity failed for $(basename "$dest_fasta")."
      fi
      gzip -cd "$tmp_gz" > "$tmp_fa"

      if [[ ! -s "$tmp_fa" ]] || ! head -n1 "$tmp_fa" | grep -q '^>' ; then
        die "[fasta] Decompressed FASTA invalid for $(basename "$dest_fasta")."
      fi

      mv -f "$tmp_fa" "$dest_fasta"
      rm -f "$tmp_gz" "$dest_gz" 2>/dev/null || true
      log "[fasta] Ready: $(basename "$dest_fasta")"
    fi

    # Build index if missing or stale
    if [[ ! -s "${dest_fasta}.fai" || "${dest_fasta}.fai" -ot "$dest_fasta" ]]; then
      log "[fasta] Indexing $(basename "$dest_fasta")"
      samtools faidx "$dest_fasta"
    fi
  }

  # Run both downloads (GRCh38 and GRCh19) possibly in parallel
  local -a pids=()

  download_and_unpack \
    "https://github.com/broadinstitute/gatk/raw/refs/heads/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz?download=" \
    "${dir}/Homo_sapiens_assembly38.fasta" &
  pids+=("$!")

  download_and_unpack \
    "https://github.com/broadinstitute/gatk/raw/refs/heads/master/src/test/resources/large/Homo_sapiens_assembly19.fasta.gz?download=" \
    "${dir}/Homo_sapiens_assembly19.fasta" &
  pids+=("$!")

  local pid
  for pid in "${pids[@]}"; do
    wait "$pid"
  done

  log "[fasta] GRCh38 and GRCh19 FASTA + .fai ready."
}

setup_jars(){
  local dir="${GENOME_DIR}/jars"; ensure_dir "$dir"
  log "==> [jars] Ensuring Beagle and bref3 JARs in $dir"

  # Define artifacts
  local url_beagle="https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar"
  local url_bref3="https://faculty.washington.edu/browning/beagle/bref3.27Feb25.75f.jar"
  local dest_beagle="${dir}/beagle.27Feb25.75f.jar"
  local dest_bref3="${dir}/bref3.27Feb25.75f.jar"

  # Helper to download with resume and verify ZIP/JAR integrity
  fetch_jar(){
    local url="$1"; local dest="$2"
    local tmp="${dest}.part"

    # Probe remote size
    local remote_size
    remote_size="$(wget --spider --server-response -O - "$url" 2>&1 \
      | awk 'tolower($0) ~ /content-length:/ {gsub("\r","",$2); print $2; exit} /Length: [0-9]+/ {gsub("\r","",$2); print $2; exit}')" || true

    # Skip if existing passes a zip integrity check and matches size (if known)
    if [[ -s "$dest" ]]; then
      if unzip -tq "$dest" >/dev/null 2>&1; then
        if [[ -n "$remote_size" ]]; then
          local local_size
          local_size=$(stat -c '%s' "$dest" 2>/dev/null || stat -f '%z' "$dest")
          if [[ "$local_size" == "$remote_size" ]]; then
            log "[jars] $(basename "$dest") present and valid; skipping."
            return 0
          fi
        else
          log "[jars] $(basename "$dest") valid (zip test); skipping."
          return 0
        fi
      else
        log "[jars] $(basename "$dest") failed integrity test; will resume/overwrite."
        mv -f "$dest" "$tmp"
      fi
    fi

    # Download with resume and retries (follow redirects just in case)
    if ! wget -c -L \
          --tries=10 \
          --waitretry=5 \
          --read-timeout=60 \
          --timeout=60 \
          --retry-connrefused \
          --no-verbose \
          -O "$tmp" "$url"; then
      log "[jars] Download failed for $(basename "$dest"). Retrying once after 10s…"
      sleep 10
      wget -c -L --tries=10 --waitretry=5 --read-timeout=60 --timeout=60 --retry-connrefused --no-verbose -O "$tmp" "$url"
    fi

    # Verify size if known
    if [[ -n "$remote_size" ]]; then
      local final_size
      final_size=$(stat -c '%s' "$tmp" 2>/dev/null || stat -f '%z' "$tmp")
      if [[ "$final_size" != "$remote_size" ]]; then
        die "[jars] Size mismatch after download for $(basename "$dest") ($final_size != $remote_size)"
      fi
    fi

    # Integrity check as ZIP
    if ! unzip -tq "$tmp" >/dev/null 2>&1; then
      die "[jars] Integrity check failed for $(basename "$dest")."
    fi

    mv -f "$tmp" "$dest"
    log "[jars] Ready: $(basename "$dest")"
  }

  # Download in parallel
  local -a pids=()
  fetch_jar "$url_beagle" "$dest_beagle" & pids+=("$!")
  fetch_jar "$url_bref3" "$dest_bref3" & pids+=("$!")
  local pid
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
}

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
    die "[ref_bcfs_b38] Missing input VCFs for: ${missing[*]} (run --only 1000G first)"
  fi

  convert_one(){
    local chr="$1"
    local in_vcf="${in_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    local base
    base=$(basename "$in_vcf" .vcf.gz)
    local out_bcf="${out_dir}/${base}.bcf"
    local tmp_bcf="${out_bcf}.part"
    local out_csi="${out_bcf}.csi"

    # Skip if up-to-date
    if [[ -s "$out_bcf" && -s "$out_csi" && "$out_bcf" -nt "$in_vcf" && "$out_csi" -nt "$out_bcf" ]]; then
      log "[ref_bcfs_b38 chr${chr}] Up-to-date; skipping."
      return 0
    fi

    # Convert to BCF (atomic)
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

    # Index with CSI
    log "[ref_bcfs_b38 chr${chr}] Indexing (CSI)"
    if ! bcftools index -f -c "$out_bcf" >/dev/null 2>&1; then
      rm -f "$out_csi" 2>/dev/null || true
      die "[ref_bcfs_b38 chr${chr}] bcftools index failed"
    fi

    # Sanity check: index exists
    if [[ ! -s "$out_csi" ]]; then
      die "[ref_bcfs_b38 chr${chr}] Missing CSI index after indexing"
    fi

    log "[ref_bcfs_b38 chr${chr}] Done: $(basename "$out_bcf"), $(basename "$out_csi")"
  }

  # Concurrency controller
  local -a pids=()
  local -i max_jobs
  max_jobs=$(( threads > 0 ? threads : 1 ))

  for chr in $(seq 1 22); do
    convert_one "$chr" &
    pids+=("$!")
    if (( ${#pids[@]} >= max_jobs )); then
      wait "${pids[0]}" || die "One of the BCF conversions failed."
      unset 'pids[0]'
      pids=(${pids[@]})
    fi
  done
  local pid
  for pid in "${pids[@]}"; do
    wait "$pid" || die "One of the BCF conversions failed."
  done

  log "[ref_bcfs_b38] All chromosomes converted and indexed."
}

setup_ref_brefs_b38(){
  local out_dir="${GENOME_DIR}/ref_brefs_b38"; ensure_dir "$out_dir"
  local in_dir="${GENOME_DIR}/1000G"
  local jar_dir="${GENOME_DIR}/jars"
  local bref3_jar="${jar_dir}/bref3.27Feb25.75f.jar"
  local java_mem="${BREF3_MEM:-4g}"
  local threads="${BREF3_THREADS:-${THREADS:-4}}"  # not all bref3 ops are multithreaded, but kept for future compat

  log "==> [ref_brefs_b38] Converting 1000G VCFs → bref3 into $out_dir"

  [[ -s "$bref3_jar" ]] || die "[ref_brefs_b38] Missing $bref3_jar (run --only jars first)"

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
    die "[ref_brefs_b38] Missing input VCFs for: ${missing[*]} (run --only 1000G first)"
  fi

  convert_one(){
    local chr="$1"
    local in_vcf="${in_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    local base
    base=$(basename "$in_vcf" .vcf.gz)
    local out_prefix="${out_dir}/${base}"
    local out_bref3="${out_prefix}.bref3"

    # Skip if up-to-date
    if [[ -s "$out_bref3" && "$out_bref3" -nt "$in_vcf" ]]; then
      log "[ref_brefs_b38 chr${chr}] Up-to-date; skipping."
      return 0
    fi

    # Work in a temp dir for atomicity
    local tmpdir
    tmpdir=$(mktemp -d)
    local cleanup_tmp
    cleanup_tmp(){ rm -rf "$tmpdir" 2>/dev/null || true; }
    trap cleanup_tmp RETURN

    local tmp_prefix="${tmpdir}/chr${chr}"
    local tmp_bref3="${tmp_prefix}.bref3"

    log "[ref_brefs_b38 chr${chr}] Running bref3 converter"
    # Primary invocation pattern used by beagle's bref3 utility
    if ! java -Xmx"$java_mem" -jar "$bref3_jar" vcf="$in_vcf" out="$tmp_prefix" >/dev/null 2>&1; then
      # Retry once (transient issues)
      log "[ref_brefs_b38 chr${chr}] First attempt failed; retrying once…"
      sleep 3
      java -Xmx"$java_mem" -jar "$bref3_jar" vcf="$in_vcf" out="$tmp_prefix" >/dev/null 2>&1 || die "[ref_brefs_b38 chr${chr}] bref3 conversion failed"
    fi

    [[ -s "$tmp_bref3" ]] || die "[ref_brefs_b38 chr${chr}] Missing output bref3 from converter"

    mv -f "$tmp_bref3" "$out_bref3"
    log "[ref_brefs_b38 chr${chr}] Done: $(basename "$out_bref3")"
  }

  # Concurrency controller (limit parallel Java processes)
  local -a pids=()
  local -i max_jobs
  max_jobs=$(( threads > 0 ? threads : 1 ))

  for chr in $(seq 1 22); do
    convert_one "$chr" &
    pids+=("$!")
    if (( ${#pids[@]} >= max_jobs )); then
      wait "${pids[0]}" || die "One of the bref3 conversions failed."
      unset 'pids[0]'
      pids=(${pids[@]})
    fi
  done
  local pid
  for pid in "${pids[@]}"; do
    wait "$pid" || die "One of the bref3 conversions failed."
  done

  log "[ref_brefs_b38] All chromosomes converted to bref3."
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
