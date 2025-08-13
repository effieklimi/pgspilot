#!/usr/bin/env bash
set -Eeuo pipefail

# scripts/pipeline/download_data.sh
# Purpose: Download and stage network-fetched genome_data artifacts only.
# NOTE: Local reference builds (BCF/bref3) have moved to build_refs.sh.

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
    local log_live="${dest}.log"
    local log_fail="${dest}.fail.log"
    local name_tbi="${name}.tbi"
    local url_tbi="${base_url}/${name_tbi}"
    local dest_tbi="${dir}/${name_tbi}"
    local tmp_tbi="${dest_tbi}.part"

    # Determine remote size if possible
    local remote_size
    remote_size="$(get_remote_size "$url")"

    # Ensure VCF exists and is complete; do not return early to allow TBI check
    local need_vcf_download=1
    if [[ -s "$dest" ]]; then
      local vcf_ok=0
      if [[ -n "$remote_size" ]]; then
        local local_size
        local_size=$(stat -c '%s' "$dest" 2>/dev/null || stat -f '%z' "$dest")
        if [[ "$local_size" == "$remote_size" ]]; then
          vcf_ok=1
        fi
      fi
      if (( vcf_ok == 0 )) && validate_header "$dest"; then
        vcf_ok=1
      fi
      if (( vcf_ok == 1 )); then
        log "[1000G chr${chr}] VCF exists and looks complete."
        need_vcf_download=0
        # Clean up any leftover partial to avoid confusing resume logs next time
        [[ -f "$tmp" ]] && rm -f "$tmp" 2>/dev/null || true
      else
        log "[1000G chr${chr}] Existing VCF seems incomplete; will resume."
        mv -f "$dest" "$tmp"
      fi
    fi

    # If tmp exists but final does not, resume into tmp
    if [[ -f "$tmp" && -n "$remote_size" ]]; then
      local part_size
      part_size=$(stat -c '%s' "$tmp" 2>/dev/null || stat -f '%z' "$tmp")
      log "[1000G chr${chr}] Resuming: $part_size / ${remote_size:-?} bytes"
    fi

    if (( need_vcf_download == 1 )); then
      # Initialize per-VCF log and track progress by monitoring the tmp file size
      : > "$log_live"
      printf '[start] %s url=%s\n' "$(date -u +%FT%TZ)" "$url" >> "$log_live"
      if [[ -n "$remote_size" ]]; then printf '[remote] size=%s\n' "$remote_size" >> "$log_live"; fi
      if [[ -f "$tmp" ]]; then
        local part_size0
        part_size0=$(stat -c '%s' "$tmp" 2>/dev/null || stat -f '%z' "$tmp")
        printf '[resume] %s part_size=%s\n' "$(date -u +%FT%TZ)" "$part_size0" >> "$log_live"
      fi

      # Helper: run wget once with a progress monitor
      run_wget_with_progress(){
        : "${1:?missing url}" "${2:?missing tmp}" "${3:?missing log}"
        wget -c \
          --tries=10 \
          --waitretry=5 \
          --read-timeout=60 \
          --timeout=60 \
          --retry-connrefused \
          --no-verbose \
          -O "$2" "$1" >>"$3" 2>&1 &
        local wpid=$!
        (
          set +e; set +u
          local last_size=0
          while kill -0 "$wpid" 2>/dev/null; do
            local size
            size=$( (wc -c < "$2") 2>/dev/null || echo 0 )
            local delta=$(( size - last_size ))
            printf '[progress] %s size=%s delta=%s\n' "$(date -u +%FT%TZ)" "$size" "$delta" >> "$3"
            last_size=$size
            sleep 10
          done
        ) & local mon_pid=$!
        wait "$wpid"; local rc=$?
        kill "$mon_pid" 2>/dev/null || true
        return $rc
      }

      if ! run_wget_with_progress "$url" "$tmp" "$log_live"; then
        log "[1000G chr${chr}] Download failed. Retrying once after 10s…"
        printf '[retry] %s waiting-before-retry\n' "$(date -u +%FT%TZ)" >> "$log_live"
        sleep 10
        run_wget_with_progress "$url" "$tmp" "$log_live" || { mv -f "$log_live" "$log_fail" 2>/dev/null || true; die "[1000G chr${chr}] Download failed on retry."; }
      fi

      # Verify completion (VCF)
      if [[ -n "$remote_size" ]]; then
        local final_size
        final_size=$(stat -c '%s' "$tmp" 2>/dev/null || stat -f '%z' "$tmp")
        if [[ "$final_size" != "$remote_size" ]]; then
          mv -f "$log_live" "$log_fail" 2>/dev/null || true
          die "[1000G chr${chr}] Size mismatch after download ($final_size != $remote_size)"
        fi
      fi

      if ! validate_header "$tmp"; then
        mv -f "$log_live" "$log_fail" 2>/dev/null || true
        die "[1000G chr${chr}] Missing VCF header signature; file may be corrupted."
      fi

      if [[ -f "$tmp" ]]; then
        mv -f "$tmp" "$dest"
        log "[1000G chr${chr}] Done → $(basename "$dest")"
        rm -f "$log_live" 2>/dev/null || true
      fi
    fi

    # --- Download matching TBI index ---
    # Determine remote size for TBI if possible
    local remote_size_tbi
    remote_size_tbi="$(get_remote_size "$url_tbi")"

    # Skip/resume logic for TBI
    if [[ -s "$dest_tbi" ]]; then
      if [[ -n "$remote_size_tbi" ]]; then
        local local_size_tbi
        local_size_tbi=$(stat -c '%s' "$dest_tbi" 2>/dev/null || stat -f '%z' "$dest_tbi")
        if [[ "$local_size_tbi" == "$remote_size_tbi" ]]; then
          log "[1000G chr${chr}] TBI exists (size OK: $local_size_tbi). Skipping."
          :
        else
          log "[1000G chr${chr}] TBI exists but size mismatch; will resume."
          mv -f "$dest_tbi" "$tmp_tbi"
        fi
      else
        log "[1000G chr${chr}] TBI exists (non-empty). Skipping."
        :
      fi
    fi

    # Resume TBI if partial exists
    if [[ -f "$tmp_tbi" && -n "$remote_size_tbi" ]]; then
      local part_size_tbi
      part_size_tbi=$(stat -c '%s' "$tmp_tbi" 2>/dev/null || stat -f '%z' "$tmp_tbi")
      log "[1000G chr${chr}] TBI resume: $part_size_tbi / ${remote_size_tbi:-?} bytes"
    fi

    if [[ ! -s "$dest_tbi" ]]; then
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

      if [[ -n "$remote_size_tbi" && -f "$tmp_tbi" ]]; then
        local final_size_tbi
        final_size_tbi=$(stat -c '%s' "$tmp_tbi" 2>/dev/null || stat -f '%z' "$tmp_tbi")
        if [[ "$final_size_tbi" != "$remote_size_tbi" ]]; then
          die "[1000G chr${chr}] TBI size mismatch after download ($final_size_tbi != $remote_size_tbi)"
        fi
      fi

      if [[ -f "$tmp_tbi" && ! -s "$tmp_tbi" ]]; then
        die "[1000G chr${chr}] TBI file is empty after download."
      fi

      if [[ -f "$tmp_tbi" ]]; then
        mv -f "$tmp_tbi" "$dest_tbi"
        log "[1000G chr${chr}] TBI Done → $(basename "$dest_tbi")"
      fi
    fi
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

  # If required chromosome maps already exist and are non-empty, skip
  have_chr_map(){
    local chr="$1"
    # Prefer exact known naming; fall back to any map containing chr token
    if find "$dir" -maxdepth 1 -type f -name "plink.chr${chr}.GRCh38.map" -size +0c | grep -q .; then
      return 0
    fi
    if find "$dir" -maxdepth 1 -type f -name "*chr${chr}*.map" -size +0c | grep -q .; then
      return 0
    fi
    return 1
  }

  local all_ok=1
  local c
  for c in $(seq 1 22); do
    if ! have_chr_map "$c"; then all_ok=0; break; fi
  done
  if (( all_ok == 1 )); then
    local existing_maps
    existing_maps=$(find "$dir" -maxdepth 1 -type f -name '*.map' | wc -l | tr -d ' ')
    log "==> [beagle_maps_b38] Required maps present (${existing_maps} files). Skipping download."
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
  # Ensure cleanup on function return. With `set -u`, expand the value now so
  # the trap doesn't reference an out-of-scope local later.
  # Expand the directory path now and embed it quoted in the trap command,
  # so that when the RETURN trap runs under `set -u`, it does not reference
  # the local variable (which would be unset at that time).
  trap "rm -rf \"$tmpdir\" 2>/dev/null || true" RETURN

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

  # Final verification: ensure required maps exist after extraction
  all_ok=1
  for c in $(seq 1 22); do
    if ! have_chr_map "$c"; then all_ok=0; break; fi
  done
  if (( all_ok == 0 )); then
    die "[beagle_maps_b38] Missing one or more required chromosome maps after extraction."
  fi

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
  # Note: GRCh19 here refers to hg19/GRCh37 content, matching historical filename
  log "==> [fasta] Ensuring GRCh38 and GRCh37 (hg19) FASTA + .fai in $dir"

  # Helper: download .gz with resume, verify size + gzip, then decompress atomically to .fasta
  download_and_unpack(){
    local url="$1"          # source .gz URL
    local dest_fasta="$2"   # final .fasta path

    local dest_gz="${dest_fasta}.gz"
    local tmp_gz="${dest_gz}.part"
    local tmp_fa="${dest_fasta}.part"
    local log_live="${dest_fasta}.log"
    local log_fail="${dest_fasta}.fail.log"

    # Best-effort remote size probe
    get_remote_size(){
      local u="$1"
      wget --spider --server-response -O - "$u" 2>&1 \
        | awk 'tolower($0) ~ /content-length:/ {gsub("\r","",$2); print $2; exit} /Length: [0-9]+/ {gsub("\r","",$2); print $2; exit}' || true
    }

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
      # Initialize per-FASTA log and track progress by monitoring the tmp file size
      : > "$log_live"
      printf '[start] %s url=%s\n' "$(date -u +%FT%TZ)" "$url" >> "$log_live"
      local remote_size
      remote_size="$(get_remote_size "$url")"
      if [[ -n "$remote_size" ]]; then printf '[remote] size=%s\n' "$remote_size" >> "$log_live"; fi
      if [[ -f "$tmp_gz" ]]; then
        local part_size0
        part_size0=$(stat -c '%s' "$tmp_gz" 2>/dev/null || stat -f '%z' "$tmp_gz")
        printf '[resume] %s part_size=%s\n' "$(date -u +%FT%TZ)" "$part_size0" >> "$log_live"
        # If the existing partial is already larger than the reported remote size, start clean
        if [[ -n "$remote_size" && "$part_size0" -gt "$remote_size" ]]; then
          log "[fasta] Partial exceeds remote size; restarting clean for $(basename "$dest_fasta")."
          rm -f "$tmp_gz"
        fi
      fi

      run_wget_with_progress(){
        : "${1:?missing url}" "${2:?missing tmp_gz}" "${3:?missing log}" "${4:-1}"
        local do_resume="$4"
        local resume_flag=()
        if [[ "$do_resume" == 1 ]]; then resume_flag=(-c); fi
        wget "${resume_flag[@]}" -L \
          --tries=10 \
          --waitretry=5 \
          --read-timeout=60 \
          --timeout=60 \
          --retry-connrefused \
          --no-verbose \
          -O "$2" "$1" >>"$3" 2>&1 &
        local wpid=$!
        (
          set +e; set +u
          local last_size=0
          while kill -0 "$wpid" 2>/dev/null; do
            local size
            size=$( (wc -c < "$2") 2>/dev/null || echo 0 )
            local delta=$(( size - last_size ))
            printf '[progress] %s size=%s delta=%s\n' "$(date -u +%FT%TZ)" "$size" "$delta" >> "$3"
            last_size=$size
            sleep 10
          done
        ) & local mon_pid=$!
        wait "$wpid"; local rc=$?
        kill "$mon_pid" 2>/dev/null || true
        return $rc
      }

      # First attempt with resume enabled
      if ! run_wget_with_progress "$url" "$tmp_gz" "$log_live" 1; then
        log "[fasta] Download failed for $(basename "$dest_fasta"). Retrying once clean after 10s…"
        printf '[retry] %s waiting-before-retry\n' "$(date -u +%FT%TZ)" >> "$log_live"
        sleep 10
        rm -f "$tmp_gz" 2>/dev/null || true
        run_wget_with_progress "$url" "$tmp_gz" "$log_live" 0 || { mv -f "$log_live" "$log_fail" 2>/dev/null || true; die "[fasta] Download failed for $(basename "$dest_fasta") on retry."; }
      fi

      # Verify size if known
      if [[ -n "$remote_size" ]]; then
        local final_size
        final_size=$(stat -c '%s' "$tmp_gz" 2>/dev/null || stat -f '%z' "$tmp_gz")
        if [[ "$final_size" != "$remote_size" ]]; then
          log "[fasta] Size mismatch after download for $(basename "$dest_fasta") ($final_size != $remote_size). Retrying once clean…"
          printf '[retry] %s size-mismatch-clean-retry\n' "$(date -u +%FT%TZ)" >> "$log_live"
          rm -f "$tmp_gz" 2>/dev/null || true
          run_wget_with_progress "$url" "$tmp_gz" "$log_live" 0 || { mv -f "$log_live" "$log_fail" 2>/dev/null || true; die "[fasta] Download failed for $(basename "$dest_fasta") after clean retry."; }
          # Re-check size if known
          final_size=$(stat -c '%s' "$tmp_gz" 2>/dev/null || stat -f '%z' "$tmp_gz")
          if [[ "$final_size" != "$remote_size" ]]; then
            mv -f "$log_live" "$log_fail" 2>/dev/null || true
            die "[fasta] Size mismatch persists for $(basename "$dest_fasta") ($final_size != $remote_size)."
          fi
        fi
      fi

      # Validate and decompress
      if ! gzip -t "$tmp_gz" 2>/dev/null; then
        log "[fasta] Gzip integrity failed for $(basename "$dest_fasta"). Retrying once clean…"
        printf '[retry] %s gzip-test-clean-retry\n' "$(date -u +%FT%TZ)" >> "$log_live"
        rm -f "$tmp_gz" 2>/dev/null || true
        run_wget_with_progress "$url" "$tmp_gz" "$log_live" 0 || { mv -f "$log_live" "$log_fail" 2>/dev/null || true; die "[fasta] Download failed for $(basename "$dest_fasta") on gzip clean retry."; }
        if ! gzip -t "$tmp_gz" 2>/dev/null; then
          mv -f "$log_live" "$log_fail" 2>/dev/null || true
          die "[fasta] Gzip integrity failed for $(basename "$dest_fasta") after clean retry."
        fi
      fi
      printf '[decompress] %s start\n' "$(date -u +%FT%TZ)" >> "$log_live"
      gzip -cd "$tmp_gz" > "$tmp_fa"

      if [[ ! -s "$tmp_fa" ]] || ! head -n1 "$tmp_fa" | grep -q '^>' ; then
        mv -f "$log_live" "$log_fail" 2>/dev/null || true
        die "[fasta] Decompressed FASTA invalid for $(basename "$dest_fasta")."
      fi

      mv -f "$tmp_fa" "$dest_fasta"
      rm -f "$tmp_gz" "$dest_gz" 2>/dev/null || true
      log "[fasta] Ready: $(basename "$dest_fasta")"
      rm -f "$log_live" 2>/dev/null || true
    fi

    # Build index if missing or stale
    if [[ ! -s "${dest_fasta}.fai" || "${dest_fasta}.fai" -ot "$dest_fasta" ]]; then
      log "[fasta] Indexing $(basename "$dest_fasta")"
      samtools faidx "$dest_fasta"
    fi
  }

  # Run both downloads (GRCh38 and GRCh37) in parallel from UCSC (stable source)
  local -a pids=()

  download_and_unpack \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" \
    "${dir}/Homo_sapiens_assembly38.fasta" &
  pids+=("$!")

  download_and_unpack \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz" \
    "${dir}/Homo_sapiens_assembly19.fasta" &
  pids+=("$!")

  local pid rc=0
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then rc=1; fi
  done
  (( rc == 0 )) || die "[fasta] One or more FASTA downloads failed."

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

run_selected(){
  local -a selected=("$@")
  local name

  ensure_dir "$GENOME_DIR"

  for name in "${selected[@]}"; do
    case "$name" in
      1000G)             setup_1000G ;;
      beagle_maps_b38)   setup_beagle_maps_b38 ;;
      bin)               ensure_dir "${GENOME_DIR}/bin" ;;
      chain)             setup_chain ;;
      eagle_maps_b38)    setup_eagle_maps_b38 ;;
      fasta)             setup_fasta ;;
      jars)              setup_jars ;;
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

log "✓ Downloads complete."
