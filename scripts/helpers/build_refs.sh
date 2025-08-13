#!/usr/bin/env bash
set -Eeuo pipefail

# Host vs container switch (mirror user.sh/add_pgs.sh behavior)
if [ -z "${INSIDE_DOCKER:-}" ]; then
  THIS_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
  PROJECT_ROOT=$(cd -- "${THIS_SCRIPT_DIR}/../../" &> /dev/null && pwd)
  DOCKER_IMAGE="${DOCKER_IMAGE:-pgspilot}"

  docker run --rm \
    -e INSIDE_DOCKER=1 \
    -v "${PROJECT_ROOT}:/app" \
    "$DOCKER_IMAGE" \
    /app/scripts/helpers/build_refs.sh "$@"
  exit 0
fi

# Container worker
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# shellcheck source=/dev/null
source "/app/scripts/pipeline/config.sh"

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
    local out_prefix="${out_dir}/${base}"
    local tmp_bcf="${out_bcf}.part.$$" # Use PID to avoid conflicts
    local out_csi="${out_bcf}.csi"

    if [[ -s "$out_bcf" && -s "$out_csi" && "$out_bcf" -nt "$in_vcf" && "$out_csi" -nt "$out_bcf" ]]; then
      log "[ref_bcfs_b38 chr${chr}] Up-to-date; skipping."
      return 0
    fi

    # Basic input sanity: ensure input has non-zero variant records; ensure index exists
    local in_size
    in_size=$(stat -c '%s' "$in_vcf" 2>/dev/null || stat -f '%z' "$in_vcf")
    if [[ -z "$in_size" || "$in_size" -le 0 ]]; then
      die "[ref_bcfs_b38 chr${chr}] Input VCF appears empty: $in_vcf"
    fi
    [[ -f "${in_vcf}.tbi" ]] || tabix -f -p vcf "$in_vcf" >/dev/null 2>&1 || true

    log "[ref_bcfs_b38 chr${chr}] Converting → BCF"
    # Progress logging similar to bref3: write live logs and monitor output size
    local log_live="${out_prefix}.bcf.log"
    local log_fail="${out_prefix}.bcf.fail.log"
    : > "$log_live"
    printf '[start] %s bcftools view → %s\n' "$(date -u +%FT%TZ)" "$(basename "$tmp_bcf")" >> "$log_live"

    # First attempt
    bcftools view -Ob -o "$tmp_bcf" --threads "$threads" "$in_vcf" >>"$log_live" 2>&1 &
    local bcf_pid=$!
    (
      set +e; set +u
      local last_size=0
      local progress_file="$tmp_bcf"
      while kill -0 "$bcf_pid" 2>/dev/null; do
        local size
        size=$( (wc -c < "$progress_file") 2>/dev/null || echo 0 )
        local delta=$(( size - last_size ))
        printf '[progress] %s size=%s delta=%s\n' "$(date -u +%FT%TZ)" "$size" "$delta" >> "$log_live"
        last_size=$size
        sleep 15
      done
    ) & local mon_pid=$!
    wait "$bcf_pid"; local rc=$?
    kill "$mon_pid" 2>/dev/null || true
    if (( rc != 0 )); then
      rm -f "$tmp_bcf" 2>/dev/null || true
      mv -f "$log_live" "$log_fail" 2>/dev/null || true
      if [[ -s "$log_fail" ]]; then
        log "[ref_bcfs_b38 chr${chr}] ---- tail of $(basename "$log_fail") ----"; tail -n 50 "$log_fail" >&2 || true; log "[ref_bcfs_b38 chr${chr}] ----------------------------------------"
      fi
      die "[ref_bcfs_b38 chr${chr}] bcftools view failed"
    fi
    if [[ ! -s "$tmp_bcf" ]]; then
      log "[ref_bcfs_b38 chr${chr}] Empty BCF after first attempt; retrying once..."
      printf '[retry] %s empty-after-first-attempt; retrying\n' "$(date -u +%FT%TZ)" >> "$log_live"
      rm -f "$tmp_bcf" 2>/dev/null || true
      bcftools view -Ob -o "$tmp_bcf" --threads "$threads" "$in_vcf" >>"$log_live" 2>&1 &
      bcf_pid=$!
      (
        set +e; set +u
        local last_size=0
        local progress_file="$tmp_bcf"
        while kill -0 "$bcf_pid" 2>/dev/null; do
          local size
          size=$( (wc -c < "$progress_file") 2>/dev/null || echo 0 )
          local delta=$(( size - last_size ))
          printf '[progress] %s size=%s delta=%s\n' "$(date -u +%FT%TZ)" "$size" "$delta" >> "$log_live"
          last_size=$size
          sleep 15
        done
      ) & mon_pid=$!
      wait "$bcf_pid"; rc=$?
      kill "$mon_pid" 2>/dev/null || true
      if (( rc != 0 )); then
        rm -f "$tmp_bcf" 2>/dev/null || true
        mv -f "$log_live" "$log_fail" 2>/dev/null || true
        if [[ -s "$log_fail" ]]; then
          log "[ref_bcfs_b38 chr${chr}] ---- tail of $(basename "$log_fail") ----"; tail -n 50 "$log_fail" >&2 || true; log "[ref_bcfs_b38 chr${chr}] ----------------------------------------"
        fi
        die "[ref_bcfs_b38 chr${chr}] bcftools view failed on retry"
      fi
      if [[ ! -s "$tmp_bcf" ]]; then
        local nvar
        nvar=$(zgrep -vc '^#' "$in_vcf" 2>/dev/null || echo 0)
        mv -f "$log_live" "$log_fail" 2>/dev/null || true
        die "[ref_bcfs_b38 chr${chr}] Empty BCF after conversion (input variants: $nvar)."
      fi
    fi

    local out_rows
    out_rows=$(bcftools view -H "$tmp_bcf" | wc -l | awk '{print $1}')
    if [[ -z "$out_rows" || "$out_rows" -eq 0 ]]; then
      rm -f "$tmp_bcf" 2>/dev/null || true
      local nvar
      nvar=$(zgrep -vc '^#' "$in_vcf" 2>/dev/null || echo 0)
      die "[ref_bcfs_b38 chr${chr}] BCF appears to have no variant rows (input variants: $nvar)."
    fi
    mv -f "$tmp_bcf" "$out_bcf"
    printf '[done] %s convert-complete size=%s\n' "$(date -u +%FT%TZ)" "$(stat -c '%s' "$out_bcf" 2>/dev/null || stat -f '%z' "$out_bcf")" >> "$log_live"

    log "[ref_bcfs_b38 chr${chr}] Indexing (CSI)"
    printf '[index] %s start\n' "$(date -u +%FT%TZ)" >> "$log_live"
    rm -f "$out_csi" 2>/dev/null || true
    
    if ! bcftools index -f -c "$out_bcf" >>"$log_live" 2>&1; then
      mv -f "$log_live" "$log_fail" 2>/dev/null || true
      die "[ref_bcfs_b38 chr${chr}] bcftools index failed"
    fi
    
    sleep 0.1
    if [[ ! -s "$out_csi" ]]; then
      # Retry indexing once more
      log "[ref_bcfs_b38 chr${chr}] Index file not found, retrying..."
      printf '[index] %s retry\n' "$(date -u +%FT%TZ)" >> "$log_live"
      sleep 1
      if ! bcftools index -f -c "$out_bcf" >>"$log_live" 2>&1; then
        mv -f "$log_live" "$log_fail" 2>/dev/null || true
        die "[ref_bcfs_b38 chr${chr}] bcftools index failed on retry"
      fi
      sleep 0.1
      [[ -s "$out_csi" ]] || die "[ref_bcfs_b38 chr${chr}] Missing CSI index after indexing retry"
    fi
    printf '[index] %s done\n' "$(date -u +%FT%TZ)" >> "$log_live"

    log "[ref_bcfs_b38 chr${chr}] Done: $(basename "$out_bcf"), $(basename "$out_csi")"
    rm -f "$log_live" 2>/dev/null || true
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
  local java_mem="${BREF3_MEM:-8g}"
  local threads="${BREF3_THREADS:-${THREADS:-2}}"
  local do_prefilter="${BREF3_FILTER:-0}"
  local max_indel_len="${BREF3_MAX_INDEL_LEN:-50}"

  to_mb(){
    local v="$1"; v=${v,,}
    if [[ "$v" =~ ^([0-9]+)g$ ]]; then echo $(( ${BASH_REMATCH[1]} * 1024 ));
    elif [[ "$v" =~ ^([0-9]+)m$ ]]; then echo ${BASH_REMATCH[1]};
    elif [[ "$v" =~ ^[0-9]+$ ]]; then echo "$v"; else echo 0; fi
  }
  local xmx_mb; xmx_mb=$(to_mb "$java_mem")
  local cgroup_limit_bytes=""; [[ -r /sys/fs/cgroup/memory.max ]] && cgroup_limit_bytes=$(cat /sys/fs/cgroup/memory.max 2>/dev/null || true)
  if [[ -z "$cgroup_limit_bytes" || "$cgroup_limit_bytes" == "max" ]]; then
    [[ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]] && cgroup_limit_bytes=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || true)
  fi
  local mem_total_mb=0
  if [[ -n "$cgroup_limit_bytes" && "$cgroup_limit_bytes" != "max" && "$cgroup_limit_bytes" =~ ^[0-9]+$ ]]; then
    mem_total_mb=$(( cgroup_limit_bytes / 1024 / 1024 ))
  else
    local mt_kb; mt_kb=$(awk '/MemTotal:/ {print $2}' /proc/meminfo 2>/dev/null || echo 0)
    [[ -n "$mt_kb" ]] && mem_total_mb=$(( mt_kb / 1024 ))
  fi
  local safe_jobs="$threads"
  if (( xmx_mb > 0 && mem_total_mb > 0 )); then
    local budget_mb=$(( mem_total_mb * 70 / 100 ))
    local cap=$(( budget_mb / xmx_mb ))
    if (( cap < 1 )); then cap=1; fi
    if (( threads > cap )); then
      log "[ref_brefs_b38] Capping parallel jobs from $threads → $cap based on memory limit ${mem_total_mb}MB and Xmx=${xmx_mb}MB"
      safe_jobs=$cap
    fi
  fi
  local -i max_jobs_global=$(( (safe_jobs>0)?safe_jobs:1 ))

  if (( mem_total_mb > 0 && xmx_mb > 0 )); then
    local budget_mb=$(( mem_total_mb * 70 / 100 ))
    local per_job_budget_mb=$(( budget_mb / max_jobs_global ))
    local headroom_mb=256
    if (( per_job_budget_mb > headroom_mb )); then
      per_job_budget_mb=$(( per_job_budget_mb - headroom_mb ))
    fi
    if (( per_job_budget_mb < 512 )); then
      per_job_budget_mb=512
    fi
    if (( xmx_mb > per_job_budget_mb )); then
      log "[ref_brefs_b38] Reducing Java Xmx from ${xmx_mb}MB → ${per_job_budget_mb}MB per job based on memory limit"
      java_mem="${per_job_budget_mb}m"
      xmx_mb=$per_job_budget_mb
    fi
  fi

  log "==> [ref_brefs_b38] Converting 1000G VCFs → bref3 into $out_dir (parallel: $max_jobs_global)"
  [[ -s "$bref3_jar" ]] || die "[ref_brefs_b38] Missing $bref3_jar (run download_data.sh --only jars first)"

  if command -v java >/dev/null 2>&1; then
    local jver_raw jver_str jmajor
    jver_raw=$(java -version 2>&1 | head -n1 || true)
    jver_str=$(printf '%s' "$jver_raw" | awk -F\" '/version/ {print $2}' || true)
    jmajor=$(printf '%s' "$jver_str" | awk -F. '{print $1}')
    # If version is like 1.8, jmajor will be 1; modern JDKs are 11,17,21…
    if [[ -n "$jmajor" && "$jmajor" =~ ^[0-9]+$ ]]; then
      if (( jmajor < 17 )); then
        die "[ref_brefs_b38] Java $jver_str detected; bref3 requires Java 17+. Update the Docker image (uses OpenJDK 17) and retry."
      fi
    fi
  else
    die "[ref_brefs_b38] 'java' not found in PATH; install OpenJDK 17 and retry."
  fi

  local bref3_mode="positional"
  local bref3_help
  bref3_help=$(java -jar "$bref3_jar" help 2>&1 || true)
  if printf '%s' "$bref3_help" | grep -qi 'vcf='; then
    bref3_mode="kv"
  else
    bref3_mode="positional"
  fi
  log "[ref_brefs_b38] Detected bref3 mode: $bref3_mode"

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
    local base
    base=$(basename "$in_vcf" .vcf.gz)
    local out_prefix="${out_dir}/${base}"
    local out_bref3="${out_prefix}.bref3"

    if [[ -s "$out_bref3" && "$out_bref3" -nt "$in_vcf" ]]; then
      log "[ref_brefs_b38 chr${chr}] Up-to-date; skipping."
      return 0
    fi

    local tmpdir
    tmpdir=$(mktemp -d)
    local cleanup_tmp
    cleanup_tmp(){ rm -rf "$tmpdir" 2>/dev/null || true; }
    trap cleanup_tmp EXIT

    local tmp_prefix="${tmpdir}/chr${chr}"
    local tmp_bref3="${tmp_prefix}.bref3"

    [[ -f "${in_vcf}.tbi" ]] || tabix -f -p vcf "$in_vcf" >/dev/null 2>&1 || true

    local vcf_for_bref="$in_vcf"
    local tmp_vcfgz="${tmpdir}/chr${chr}.biallelic.snvindel.vcf.gz"
    if [[ "$do_prefilter" == "1" ]]; then
      if ! command -v bcftools >/dev/null 2>&1; then
        die "[ref_brefs_b38 chr${chr}] bcftools not found; required for BREF3_FILTER=1"
      fi
      log "[ref_brefs_b38 chr${chr}] Prefiltering to biallelic SNVs/short INDELs (≤${max_indel_len}bp)"
      if ! bcftools view -m2 -M2 -v snps,indels "$in_vcf" \
          | bcftools filter -e "strlen(REF) > ${max_indel_len} || strlen(ALT) > ${max_indel_len} || ALT ~ '<'" - \
          | bgzip -c > "$tmp_vcfgz"; then
        die "[ref_brefs_b38 chr${chr}] Prefiltering with bcftools failed"
      fi
      tabix -f -p vcf "$tmp_vcfgz" >/dev/null 2>&1 || true
      vcf_for_bref="$tmp_vcfgz"
    fi

    log "[ref_brefs_b38 chr${chr}] Running bref3 converter"
    local log_live="${out_prefix}.bref3.log"
    local log_fail="${out_prefix}.bref3.fail.log"
    local run_bref3
    if [[ "$bref3_mode" == "kv" ]]; then
      run_bref3=(java -XX:+ExitOnOutOfMemoryError -Xmx"$java_mem" -jar "$bref3_jar" vcf="$vcf_for_bref" out="$tmp_prefix")
    else
      if [[ "$vcf_for_bref" == *.gz ]]; then
        run_bref3=(bash -c 'zcat "$1" | java -XX:+ExitOnOutOfMemoryError -Xmx"$2" -jar "$3" > "$4"' _ "$vcf_for_bref" "$java_mem" "$bref3_jar" "$tmp_bref3")
      else
        run_bref3=(bash -c 'cat "$1" | java -XX:+ExitOnOutOfMemoryError -Xmx"$2" -jar "$3" > "$4"' _ "$vcf_for_bref" "$java_mem" "$bref3_jar" "$tmp_bref3")
      fi
    fi
    : > "$log_live"
    "${run_bref3[@]}" >>"$log_live" 2>&1 &
    local bref3_pid=$!
    (
      set +e
      set +u
      local last_size=0
      local progress_file
      if [[ "$bref3_mode" == "kv" ]]; then
        progress_file="${tmp_prefix}.bref3"
      else
        progress_file="$tmp_bref3"
      fi
      while kill -0 "$bref3_pid" 2>/dev/null; do
        local size
        size=$( (wc -c < "$progress_file") 2>/dev/null || echo 0 )
        local delta=$(( size - last_size ))
        printf '[progress] %s size=%s delta=%s\n' "$(date -u +%FT%TZ)" "$size" "$delta" >> "$log_live"
        last_size=$size
        sleep 15
      done
    ) & local mon_pid=$!
    wait "$bref3_pid"
    local rc=$?
    kill "$mon_pid" 2>/dev/null || true
    if (( rc != 0 )); then
      log "[ref_brefs_b38 chr${chr}] First attempt failed; preparing minimal-filter retry…"
      sleep 1
      local retry_vcf="$vcf_for_bref"
      local minimal_filtered="${tmpdir}/chr${chr}.minimal.vcf.gz"
      if [[ "$do_prefilter" == "0" ]] && command -v bcftools >/dev/null 2>&1; then
        log "[ref_brefs_b38 chr${chr}] Applying minimal filter to drop symbolic/breakend alleles for retry"
        if ! bcftools filter -e "ALT ~ '<' || ALT == '*' || ALT ~ '\\[' || ALT ~ '\\]'" "$vcf_for_bref" | bgzip -c > "$minimal_filtered" 2>>"$log_live"; then
          log "[ref_brefs_b38 chr${chr}] Minimal filter step failed; retrying without filter"
        else
          tabix -f -p vcf "$minimal_filtered" >/dev/null 2>&1 || true
          retry_vcf="$minimal_filtered"
        fi
      fi
      if [[ "$bref3_mode" == "kv" ]]; then
        run_bref3=(java -XX:+ExitOnOutOfMemoryError -Xmx"$java_mem" -jar "$bref3_jar" vcf="$retry_vcf" out="$tmp_prefix")
      else
        if [[ "$retry_vcf" == *.gz ]]; then
          run_bref3=(bash -c 'zcat "$1" | java -XX:+ExitOnOutOfMemoryError -Xmx"$2" -jar "$3" > "$4"' _ "$retry_vcf" "$java_mem" "$bref3_jar" "$tmp_bref3")
        else
          run_bref3=(bash -c 'cat "$1" | java -XX:+ExitOnOutOfMemoryError -Xmx"$2" -jar "$3" > "$4"' _ "$retry_vcf" "$java_mem" "$bref3_jar" "$tmp_bref3")
        fi
      fi
      "${run_bref3[@]}" >>"$log_live" 2>&1 &
      bref3_pid=$!
      (
        set +e
        set +u
        local last_size=0
        local progress_file
        if [[ "$bref3_mode" == "kv" ]]; then
          progress_file="${tmp_prefix}.bref3"
        else
          progress_file="$tmp_bref3"
        fi
        while kill -0 "$bref3_pid" 2>/dev/null; do
          local size
          size=$( (wc -c < "$progress_file") 2>/dev/null || echo 0 )
          local delta=$(( size - last_size ))
          printf '[progress] %s size=%s delta=%s\n' "$(date -u +%FT%TZ)" "$size" "$delta" >> "$log_live"
          last_size=$size
          sleep 15
        done
      ) & mon_pid=$!
      wait "$bref3_pid"; rc=$?
      kill "$mon_pid" 2>/dev/null || true
      if (( rc != 0 )); then
        # Preserve logs for debugging
        mv -f "$log_live" "$log_fail" 2>/dev/null || true
        # Surface last lines for quick diagnosis
        if [[ -s "$log_fail" ]]; then
          log "[ref_brefs_b38 chr${chr}] ---- tail of $(basename "$log_fail") ----"
          tail -n 50 "$log_fail" >&2 || true
          log "[ref_brefs_b38 chr${chr}] ----------------------------------------"
        fi
        die "[ref_brefs_b38 chr${chr}] bref3 conversion failed (see $(basename "$log_fail"))"
      fi
    fi

    [[ -s "$tmp_bref3" ]] || { die "[ref_brefs_b38 chr${chr}] Missing output bref3 from converter"; }

    mv -f "$tmp_bref3" "$out_bref3"
    rm -f "$log_live" 2>/dev/null || true
    log "[ref_brefs_b38 chr${chr}] Done: $(basename "$out_bref3")"
  }

  local -a pids=(); local -i max_jobs=$max_jobs_global
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
