#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pgs_fetch_hmpos_and_filter.py — Fetch hmPOS for a PGS ID (prefer GRCh38), liftover to GRCh38 if
needed, then apply subpopulation-specific MAF filtering of palindromic SNPs for ALL subpops
(EUR, AFR, EAS, AMR, SAS). Emits one TSV per subpop (+ metadata).

What it does:
  1) Resolve and download the PGS Catalog *harmonized* scoring file (hmPOS), preferring GRCh38.
     If GRCh38 is unavailable, it fetches GRCh37 and liftover→GRCh38 once.
  2) Use a GRCh38 FASTA (fixed path) to fetch the reference base and orient alleles.
  3) Convert OR/HR -> ln(OR/HR); leave beta/logOR unchanged.
  4) Build a candidate set on GRCh38 (ref/alt/effect_allele/beta) ONCE.
  5) For each subpopulation (EUR/AFR/EAS/AMR/SAS), apply subpop MAF filtering for palindromic SNPs.
  6) Write: <out_dir>/<PGSID>.<SUBPOP>.b38.tsv[.gz] (+ .meta.json per subpop)

Dependencies:
    pip install "pandas>=2.0" "numpy>=1.23" "requests>=2.28" "pyfaidx>=0.7" "pyliftover>=0.4"

Fixed paths (override via ENV if ever needed):
    FASTA38_PATH        = /app/genome_data/fasta/Homo_sapiens_assembly38.fasta
    CHAIN37TO38_PATH    = /app/genome_data/chain/hg19ToHg38.over.chain.gz
    MAF_DIR_DEFAULT     = /app/pca_model
    OUT_DIR_DEFAULT     = /app/pgs/weights/harmonized
"""

from __future__ import annotations

import argparse, io, os, sys, json, gzip, logging, math, hashlib, tempfile, datetime as dt
from typing import Dict, Optional, Tuple, List

import numpy as np
import pandas as pd
import requests
from pyfaidx import Fasta

try:
    from pyliftover import LiftOver
except Exception:
    LiftOver = None  # raise later only if we actually need liftover

# ── Fixed locations (no CLI flags; overridable via ENV) ─────────────────

ENV = os.environ
FASTA38_PATH        = ENV.get("PGS_FASTA38",      "/app/genome_data/fasta/Homo_sapiens_assembly38.fasta")
CHAIN37TO38_PATH    = ENV.get("PGS_CHAIN_37_38",  "/app/genome_data/chain/hg19ToHg38.over.chain.gz")
MAF_DIR_DEFAULT     = ENV.get("PGS_MAF_DIR",      "/app/pca_model")
OUT_DIR_DEFAULT     = ENV.get("PGS_OUT_DIR",      "/app/pgs/weights/harmonized")

AMBIG_MIN_DEFAULT = float(ENV.get("PGS_AMBIG_MAF_MIN", "0.45"))
AMBIG_MAX_DEFAULT = float(ENV.get("PGS_AMBIG_MAF_MAX", "0.55"))

SUBPOPS = ("EUR","AFR","EAS","AMR","SAS")

LOG = logging.getLogger("pgs_fetch_hmpos")
VALID_A = {"A","C","G","T"}
COMPLEMENT = str.maketrans("ACGTacgt","TGCAtgca")

# ── Utils ───────────────────────────────────────────────────────────────

def setup_logging(verbosity: int) -> None:
    level = logging.WARNING if verbosity <= 0 else (logging.INFO if verbosity == 1 else logging.DEBUG)
    logging.basicConfig(level=level, format="%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%dT%H:%M:%S")

def to_ucsc(ch: str) -> str:
    s = str(ch).strip()
    return s if s.startswith("chr") else ("chr"+s if s else s)

def palindromic(a1: str, a2: str) -> bool:
    s = {a1,a2}
    return s == {"A","T"} or s == {"C","G"}

def parse_header_meta(lines: List[str]) -> Dict[str,str]:
    meta = {}
    for line in lines:
        if not line.startswith("#"): break
        if "=" in line:
            k, v = line.lstrip("#").split("=", 1)
            meta[k.strip().lower()] = v.strip()
    return meta

def md5_bytes(b: bytes) -> str:
    h = hashlib.md5(); h.update(b); return h.hexdigest()

def md5_file(path: str) -> Optional[str]:
    try:
        with open(path, "rb") as f:
            h = hashlib.md5()
            for chunk in iter(lambda: f.read(1024*1024), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None

def atomic_write(path: str, data: str) -> None:
    d = os.path.dirname(path) or "."
    os.makedirs(d, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=d, encoding="utf-8", newline="") as tmp:
        tmp.write(data)
        tmp_path = tmp.name
    os.replace(tmp_path, path)

def atomic_write_df_tsv(df: pd.DataFrame, path: str, gzip_out: bool) -> str:
    final = path + ".gz" if gzip_out and not path.endswith(".gz") else path
    d = os.path.dirname(final) or "."
    os.makedirs(d, exist_ok=True)
    mode = "wb" if final.endswith(".gz") else "w"
    with tempfile.NamedTemporaryFile(mode, delete=False, dir=d) as tmp:
        tmp_path = tmp.name
        if final.endswith(".gz"):
            with gzip.GzipFile(fileobj=tmp, mode="wb") as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as txt:
                    df.to_csv(txt, sep="\t", index=False)
        else:
            with open(tmp_path, "w", encoding="utf-8", newline="") as txt:
                df.to_csv(txt, sep="\t", index=False)
    os.replace(tmp_path, final)
    return final

# ── PGS Catalog hmPOS resolution ───────────────────────────────────────

REST_BASE = "https://www.pgscatalog.org/rest"
FTP_BASE  = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"

def _find_hmpos_url_in_json(obj, pgs_id: str, build: str) -> Optional[str]:
    want_suffix = f"{pgs_id}_hmPOS_{build}.txt.gz"
    hit = None
    def walk(x):
        nonlocal hit
        if hit is not None: return
        if isinstance(x, dict):
            for _, v in x.items(): walk(v)
        elif isinstance(x, list):
            for v in x: walk(v)
        elif isinstance(x, str):
            if x.endswith(want_suffix) and "/Harmonized/" in x:
                hit = x
    walk(obj)
    return hit

def resolve_hmpos_urls(pgs_id: str) -> Tuple[Optional[str], Optional[str], str]:
    """Return (url_b38, url_b37, resolution_mode)."""
    try:
        r = requests.get(f"{REST_BASE}/score/{pgs_id}", timeout=30)
        if r.ok:
            j = r.json()
            u38 = _find_hmpos_url_in_json(j, pgs_id, "GRCh38")
            u37 = _find_hmpos_url_in_json(j, pgs_id, "GRCh37")
            if u38 or u37:
                return u38, u37, "REST"
    except Exception as e:
        LOG.debug("REST lookup failed (%s): %s", pgs_id, e)
    # FTP fallbacks
    u38 = f"{FTP_BASE}/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_GRCh38.txt.gz"
    u37 = f"{FTP_BASE}/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_GRCh37.txt.gz"
    return u38, u37, "FTP"

# ── MAF handling ───────────────────────────────────────────────────────

def maf_path_for(subpop: str, maf_dir: str) -> str:
    return os.path.join(maf_dir, f"maf_{subpop.upper()}.grch38.tsv.gz")

def load_maf_table(path: str) -> Dict[Tuple[str, int, str, str], float]:
    if not os.path.exists(path):
        raise SystemExit(f"[ERROR] MAF file not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str).rename(columns=str.lower)
    def pick(*names):
        for n in names:
            if n in df.columns: return n
        return None
    c_chr = pick("chr","chrom","chromosome")
    c_pos = pick("pos","position","chr_position")
    c_ref = pick("ref","ref_allele","allele1","a1")
    c_alt = pick("alt","alt_allele","allele2","a2")
    c_maf = pick("maf")
    c_af  = pick("af","a1_freq","alt_af","af_alt","allele2_af")
    if not all([c_chr,c_pos,c_ref,c_alt]) or (c_maf is None and c_af is None):
        raise SystemExit("[ERROR] MAF file must have chr,pos,ref,alt and maf or af.")
    out: Dict[Tuple[str,int,str,str], float] = {}
    for _, r in df.iterrows():
        try:
            chrom = to_ucsc(r[c_chr])
            pos   = int(float(r[c_pos]))
            ref   = str(r[c_ref]).upper()
            alt   = str(r[c_alt]).upper()
            maf   = float(r[c_maf]) if (c_maf and str(r[c_maf]) != "") else float(r[c_af])
            maf   = maf if maf <= 0.5 else 1.0 - maf
        except Exception:
            continue
        if ref in VALID_A and alt in VALID_A and ref != alt:
            out[(chrom,pos,ref,alt)] = maf
            out[(chrom,pos,alt,ref)] = maf
    return out

# ── Core helpers ───────────────────────────────────────────────────────

class InputError(RuntimeError): ...
class ConfigError(RuntimeError): ...

def fetch_hmpos_dataframe(url: str) -> Tuple[pd.DataFrame, Dict[str,str], dict]:
    LOG.info("Downloading: %s", url)
    r = requests.get(url, timeout=60, stream=True)
    if not r.ok:
        raise InputError(f"Could not download hmPOS file (HTTP {r.status_code}): {url}")
    content = r.content
    md5 = md5_bytes(content)
    head_lines = []
    with gzip.GzipFile(fileobj=io.BytesIO(content), mode="rb") as gz:
        while True:
            pos = gz.tell()
            line = gz.readline()
            if not line: break
            s = line.decode("utf-8", errors="replace")
            if s.startswith("#"):
                head_lines.append(s)
            else:
                gz.seek(pos)
                tail = gz.read()
                joined = ("".join(head_lines).encode("utf-8")) + line + tail
                # Ignore header metadata lines that start with '#'
                df = pd.read_csv(io.BytesIO(joined), sep="\t", dtype=str, comment="#", low_memory=False)
                return df.fillna(""), parse_header_meta(head_lines), {"url": url, "md5": md5}
    raise InputError("Empty or malformed hmPOS file")

def detect_cols(df: pd.DataFrame) -> Tuple[str,str,str,str,str]:
    cols = {c.lower(): c for c in df.columns}
    def have(*names):
        for n in names:
            if n.lower() in cols: return cols[n.lower()]
        return None
    chr_c = have("hm_chr","chr_name","chromosome","chr")
    pos_c = have("hm_pos","chr_position","position","pos")
    ea_c  = have("effect_allele","ea")
    oa_c  = have("other_allele","oa","reference_allele","ref_allele")
    wt_c  = have("effect_weight","beta","logor","log_or","or","odds_ratio","hr","hazard_ratio")
    if None in (chr_c, pos_c, ea_c, oa_c, wt_c):
        raise InputError("hmPOS/PGS file missing required columns (hm_chr/pos, effect_allele/other_allele, weight).")
    return chr_c, pos_c, ea_c, oa_c, wt_c

def liftover_pos(lifter: Optional["LiftOver"], chrom: str, pos1: int) -> Tuple[Optional[str], Optional[int], str]:
    if lifter is None:
        return chrom, pos1, "+"
    res = lifter.convert_coordinate(chrom, pos1 - 1)
    if not res:
        return None, None, "+"
    new_ch, new_pos0, strand, _ = res[0]
    if not str(new_ch).startswith("chr"):
        new_ch = "chr" + str(new_ch)
    return str(new_ch), int(new_pos0) + 1, strand

def safe_float(x: str) -> Optional[float]:
    try:
        v = float(x)
        if math.isfinite(v): return v
    except Exception:
        pass
    return None

# ── Main ───────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser("Fetch hmPOS for a PGS ID, liftover to GRCh38 if needed, then subpop-MAF filter palindromic SNPs for ALL subpops.")
    p.add_argument("--pgs-id", required=True, help="PGS identifier, e.g. PGS000001")
    p.add_argument("--maf-dir", default=MAF_DIR_DEFAULT, help="Directory containing maf_<POP>.grch38.tsv.gz")
    p.add_argument("--ambig-min", type=float, default=AMBIG_MIN_DEFAULT, help="Lower bound of ambiguous MAF window")
    p.add_argument("--ambig-max", type=float, default=AMBIG_MAX_DEFAULT, help="Upper bound of ambiguous MAF window")
    p.add_argument("--out-dir", default=OUT_DIR_DEFAULT, help="Where to write outputs")
    p.add_argument("--gzip", action="store_true", help="Gzip the output TSVs")
    g = p.add_mutually_exclusive_group()
    g.add_argument("-q","--quiet", action="store_true", help="Only warnings and errors")
    g.add_argument("-v","--verbose", action="count", default=1, help="Increase verbosity (repeat for DEBUG)")
    args = p.parse_args()

    setup_logging(0 if args.quiet else (args.verbose if isinstance(args.verbose, int) else 1))

    # Resolve and fetch hmPOS (prefer 38, else 37)
    pgs_id = args.pgs_id.strip()
    url38, url37, resolved_via = resolve_hmpos_urls(pgs_id)

    src_build = None
    try:
        if not url38:
            raise InputError("No GRCh38 hmPOS URL found via resolver.")
        df, meta_hdr, transfer = fetch_hmpos_dataframe(url38)
        src_build = "GRCh38"
    except Exception as e38:
        LOG.info("Could not fetch GRCh38 hmPOS (%s). Trying GRCh37.", e38)
        if not url37:
            raise SystemExit(f"[ERROR] Neither GRCh38 nor GRCh37 hmPOS URL available for {pgs_id}.")
        df, meta_hdr, transfer = fetch_hmpos_dataframe(url37)
        src_build = "GRCh37"

    # Fixed FASTA & chain
    if not os.path.exists(FASTA38_PATH) or not os.path.exists(FASTA38_PATH + ".fai"):
        raise SystemExit(f"[ERROR] FASTA or index missing: {FASTA38_PATH} (.fai required)")
    fa = Fasta(FASTA38_PATH, rebuild=False, sequence_always_upper=True)
    fasta38_md5 = md5_file(FASTA38_PATH)

    chain_md5 = None
    lifter = None
    liftover_applied = (src_build == "GRCh37")
    if liftover_applied:
        if LiftOver is None:
            raise SystemExit("[ERROR] pyliftover not installed but required for GRCh37→GRCh38 liftover.")
        if not os.path.exists(CHAIN37TO38_PATH):
            raise SystemExit(f"[ERROR] Liftover chain not found: {CHAIN37TO38_PATH}")
        lifter = LiftOver(CHAIN37TO38_PATH)
        chain_md5 = md5_file(CHAIN37TO38_PATH)

    # Detect columns and weight type
    chr_c, pos_c, ea_c, oa_c, wt_c = detect_cols(df)
    weight_type_src = (meta_hdr.get("weight_type") or "").lower()

    # ── Build candidate rows ONCE on GRCh38 ─────────────────────────────
    base_qc = dict(
        input_rows=int(df.shape[0]),
        dropped_bad_pos=0,
        dropped_liftover=0,
        dropped_allele_invalid=0,
        dropped_ref_fetch=0,
        dropped_ref_mismatch=0,
        dropped_weight_invalid=0
    )
    candidates: List[Tuple[str,int,str,str,str,float,bool]] = []  # (chr,pos,ref,alt,eff,beta,is_pal)

    for _, r in df.iterrows():
        chrom = to_ucsc(str(r[chr_c]))
        try:
            pos = int(float(r[pos_c]))
        except Exception:
            base_qc["dropped_bad_pos"] += 1
            continue

        strand = "+"
        if liftover_applied:
            ch38, p38, strand = liftover_pos(lifter, chrom, pos)
            if ch38 is None or p38 is None:
                base_qc["dropped_liftover"] += 1
                continue
            chrom, pos = ch38, p38

        ea = str(r[ea_c]).upper()
        oa = str(r[oa_c]).upper()
        if ea not in VALID_A or oa not in VALID_A or ea == oa:
            base_qc["dropped_allele_invalid"] += 1
            continue
        if strand == "-":
            ea = ea.translate(COMPLEMENT)
            oa = oa.translate(COMPLEMENT)

        # reference base
        try:
            ref = fa[chrom][pos-1:pos].seq.upper()
        except Exception:
            ref = None
        if ref not in VALID_A:
            base_qc["dropped_ref_fetch"] += 1
            continue

        # orient vs ref
        if ref == ea:
            alt, eff = oa, ea
        elif ref == oa:
            alt, eff = ea, ea
        else:
            base_qc["dropped_ref_mismatch"] += 1
            continue

        # effect size
        w = safe_float(str(r[wt_c]).strip())
        if w is None:
            base_qc["dropped_weight_invalid"] += 1
            continue
        wt_type = weight_type_src or wt_c.lower()
        if any(t in wt_type for t in (" or", "odds_ratio", "odds ratio", "hr", "hazard")) or wt_type in {"or","odds_ratio","hr"}:
            if w <= 0:
                base_qc["dropped_weight_invalid"] += 1
                continue
            w = float(np.log(w))

        candidates.append((chrom, pos, ref, alt, eff, float(w), palindromic(ref, alt)))

    if not candidates:
        raise SystemExit("[ERROR] No variants survived base harmonization.")

    # ── Per-subpopulation filtering + write ─────────────────────────────
    os.makedirs(args.out_dir, exist_ok=True)
    outputs = []

    for sub in SUBPOPS:
        maf_path = maf_path_for(sub, args.maf_dir)
        if not os.path.exists(maf_path):
            LOG.warning("MAF file missing for %s: %s (skipping this subpop)", sub, maf_path)
            continue
        maf = load_maf_table(maf_path)

        sub_qc = dict(
            subpopulation=sub,
            n_candidates=len(candidates),
            dropped_pal_maf_missing=0,
            dropped_pal_maf_ambig=0
        )

        out_rows = []
        for chrom, pos, ref, alt, eff, beta, is_pal in candidates:
            if is_pal:
                m = maf.get((chrom, pos, ref, alt))
                if m is None:
                    sub_qc["dropped_pal_maf_missing"] += 1
                    continue
                if args.ambig_min <= m <= args.ambig_max:
                    sub_qc["dropped_pal_maf_ambig"] += 1
                    continue
            out_rows.append((chrom, pos, ref, alt, eff, beta))

        if not out_rows:
            LOG.warning("All variants dropped after subpop filtering for %s.", sub)
            continue

        out_df = pd.DataFrame(out_rows, columns=["chr","pos","ref","alt","effect_allele","beta"])
        out_tsv = os.path.join(args.out_dir, f"{pgs_id}.{sub}.b38.tsv")
        final_tsv = atomic_write_df_tsv(out_df, out_tsv, gzip_out=bool(args.gzip))

        meta_out = {
            "pgs_id": pgs_id,
            "subpopulation": sub,
            "hmpos_source_build": src_build,               # GRCh38 or GRCh37
            "liftover_applied": bool(liftover_applied),
            "liftover_chain": CHAIN37TO38_PATH if liftover_applied else "",
            "liftover_chain_md5": md5_file(CHAIN37TO38_PATH) if liftover_applied else "",
            "hmpos_source_url": transfer["url"],
            "hmpos_download_md5": transfer["md5"],
            "resolved_via": resolved_via,                  # REST or FTP
            "fasta38": FASTA38_PATH,
            "fasta38_md5": md5_file(FASTA38_PATH),
            "maf_file": os.path.abspath(maf_path),
            "ambiguous_maf_window": [args.ambig_min, args.ambig_max],
            "weight_type_src": (meta_hdr.get("weight_type") or wt_c),
            "n_input_rows": base_qc["input_rows"],
            "n_base_kept": len(candidates),
            "n_kept_rows": int(out_df.shape[0]),
            "base_drop_counters": base_qc,
            "subpop_drop_counters": sub_qc,
            "out_path": os.path.abspath(final_tsv),
            "timestamp_utc": dt.datetime.utcnow().isoformat(timespec="seconds")+"Z",
            "notes": "Used hmPOS (GRCh38 preferred; lifted from GRCh37 if needed). Applied subpopulation-specific MAF filtering to palindromic SNPs."
        }
        atomic_write(final_tsv + ".meta.json", json.dumps(meta_out, indent=2))
        LOG.info("[OK] %s %s → %s (n=%d)", pgs_id, sub, final_tsv, out_df.shape[0])
        outputs.append(final_tsv)

    if not outputs:
        raise SystemExit("[ERROR] No subpopulation outputs were produced (check MAF files).")
    LOG.info("Done. Wrote %d subpopulation files.", len(outputs))

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
