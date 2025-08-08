#!/usr/bin/env python3
"""
pgs_harmonize_mvp.py — GRCh38-only, MAF-aware, per-subpopulation runs.

Run this ONCE per subpopulation (EUR/AFR/EAS/AMR/SAS). The script:
  • Reads a PGS Catalog scoring file (.txt/.tsv[.gz]) on GRCh37 or GRCh38.
  • Liftover to GRCh38 (fixed chain path).
  • Checks and orients alleles to GRCh38 reference FASTA (fixed path).
  • Converts OR/HR to ln(OR/HR); leaves beta/logOR as-is.
  • Handles palindromic SNPs (A/T, C/G) using a per-subpop MAF file on GRCh38:
      - Keep only if MAF < ambig_min or MAF > ambig_max (default 0.45–0.55).
      - Drop if MAF missing.
  • Writes: weights_hm/<PGSID>.<SUBPOP>.b38.tsv (+ .meta.json)
  • Appends one row to a fixed TSV registry.

Minimal dependencies:
    pip install pandas numpy pyfaidx pyliftover
"""

import argparse, gzip, io, os, sys, json, datetime
from typing import Tuple, Optional, Dict, Tuple as TupleType

import numpy as np
import pandas as pd
from pyfaidx import Fasta

try:
    from pyliftover import LiftOver
except ImportError:
    LiftOver = None

# ── FIXED PATHS (edit once for your deployment) ────────────────────────
FASTA38_PATH      = "/ref/GRCh38.no_alt.fa"                       # must have .fai
CHAIN37TO38_PATH  = "/ref/hg19ToHg38.over.chain.gz"               # required for GRCh37/HG19 input
MAF_DIR           = "/ref/maf"                                    # contains maf_EUR.grch38.tsv.gz, etc.
REGISTRY_PATH     = "/app/weights_hm/pgs_harmonization_registry.tsv"
DEFAULT_OUT_DIR   = "/app/weights_hm"

# Ambiguity window defaults (can override via CLI)
AMBIG_MAF_MIN_DEFAULT = 0.45
AMBIG_MAF_MAX_DEFAULT = 0.55

VALID_A = {"A","C","G","T"}
COMPLEMENT = str.maketrans("ACGTacgt","TGCAtgca")

# ── helpers ────────────────────────────────────────────────────────────

def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def read_text(path: str) -> str:
    with open_text(path) as f:
        return f.read()

def parse_meta(raw: str) -> dict:
    meta = {}
    for line in raw.splitlines():
        if not line.startswith("#"): break
        if "=" in line:
            k, v = line.lstrip("#").split("=", 1)
            meta[k.strip().lower()] = v.strip()
    return meta

def to_ucsc(ch: str) -> str:
    s = str(ch)
    return s if s.startswith("chr") else f"chr{s}"

def liftover_pos(lifter: Optional["LiftOver"], chrom: str, pos1: int) -> Tuple[Optional[str], Optional[int], str]:
    if lifter is None:
        return chrom, pos1, "+"
    res = lifter.convert_coordinate(chrom, pos1-1)
    if not res:
        return None, None, "+"
    new_ch, new_pos0, strand, _ = res[0]
    if not str(new_ch).startswith("chr"):
        new_ch = "chr" + str(new_ch)
    return new_ch, int(new_pos0)+1, strand

def fetch_ref(fa: Fasta, chrom: str, pos1: int) -> Optional[str]:
    try:
        return fa[chrom][pos1-1:pos1].seq.upper()
    except Exception:
        return None

def maf_path_for(subpop: str) -> str:
    subpop = subpop.upper()
    return os.path.join(MAF_DIR, f"maf_{subpop}.grch38.tsv.gz")

def load_maf_table(path: str) -> Dict[TupleType[str, int, str, str], float]:
    with open_text(path) as fh:
        df = pd.read_csv(fh, sep="\t", dtype=str).rename(columns=str.lower)
    def pick(*names):
        for n in names:
            if n in df.columns: return n
        return None
    c_chr = pick("chr","chrom","chromosome")
    c_pos = pick("pos","position")
    c_ref = pick("ref","ref_allele","allele1")
    c_alt = pick("alt","alt_allele","allele2")
    c_maf = pick("maf")
    c_af  = pick("af","a1_freq","alt_af","af_alt","allele2_af")
    if not all([c_chr,c_pos,c_ref,c_alt]) or (c_maf is None and c_af is None):
        raise SystemExit("[ERROR] MAF file must have chr,pos,ref,alt and maf or af.")
    out: Dict[TupleType[str,int,str,str], float] = {}
    for _, r in df.iterrows():
        try:
            chrom = r[c_chr]; chrom = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
            pos   = int(float(r[c_pos]))
            ref   = str(r[c_ref]).upper()
            alt   = str(r[c_alt]).upper()
            maf   = float(r[c_maf]) if (c_maf and r[c_maf] != "") else float(r[c_af])
            maf   = maf if maf <= 0.5 else 1.0 - maf
        except Exception:
            continue
        if ref in VALID_A and alt in VALID_A and ref != alt:
            out[(chrom,pos,ref,alt)] = maf
            out[(chrom,pos,alt,ref)] = maf  # order-agnostic
    return out

def append_registry(row: dict):
    cols = [
        "timestamp_iso", "pgs_id", "pgs_url", "trait",
        "weight_type_src", "genome_build_src",
        "ancestry", "maf_file", "ambiguous_maf_window",
        "source_file", "out_path", "n_input_rows", "n_harmonised_rows"
    ]
    os.makedirs(os.path.dirname(REGISTRY_PATH), exist_ok=True) if os.path.dirname(REGISTRY_PATH) else None
    new = not os.path.exists(REGISTRY_PATH)
    with open(REGISTRY_PATH, "a") as fh:
        if new:
            fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")

# ── main ───────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser("PGS harmoniser (GRCh38, per-subpop, fixed refs)")
    p.add_argument("--in", dest="inp", required=True, help="PGS scoring file (.txt/.tsv[.gz])")
    p.add_argument("--subpop", required=True, choices=["EUR","AFR","EAS","AMR","SAS"], help="Subpopulation for MAF handling and output naming")
    p.add_argument("--out-dir", default=DEFAULT_OUT_DIR, help=f"Output directory (default: {DEFAULT_OUT_DIR})")
    p.add_argument("--ambig-maf-min", type=float, default=AMBIG_MAF_MIN_DEFAULT, help="Lower bound for ambiguous palindrome drop window")
    p.add_argument("--ambig-maf-max", type=float, default=AMBIG_MAF_MAX_DEFAULT, help="Upper bound for ambiguous palindrome drop window")
    args = p.parse_args()

    # Load scoring file
    raw  = read_text(args.inp)
    meta = parse_meta(raw)
    df   = pd.read_csv(io.StringIO(raw), sep="\t", comment="#", dtype=str).fillna("")

    # Columns (v2 + common v1 variants)
    COLS = {
        "chr": ["chr_name","chromosome","chr"],
        "pos": ["chr_position","position","pos"],
        "ea" : ["effect_allele","EA"],
        "oa" : ["other_allele","OA"],
        "wt" : ["effect_weight","beta","OR","HR","logOR","log_or","odds_ratio"],
    }
    pick = lambda cands: next((c for c in cands if c in df.columns), None)
    chr_c, pos_c, ea_c, oa_c, wt_c = [pick(COLS[k]) for k in ("chr","pos","ea","oa","wt")]
    if None in (chr_c, pos_c, ea_c, oa_c, wt_c):
        sys.exit("[ERROR] Scoring file lacks required columns (chr/pos/EA/OA/weight).")

    # Build / liftover
    build  = (meta.get("genome_build","") or "").upper()
    lifter = None
    if build in ("GRCH37","HG19","B37"):
        if LiftOver is None:
            sys.exit("pyliftover not installed – cannot liftover GRCh37/HG19 input.")
        if not os.path.exists(CHAIN37TO38_PATH):
            sys.exit(f"[ERROR] liftover chain not found: {CHAIN37TO38_PATH}")
        lifter = LiftOver(CHAIN37TO38_PATH)

    # References & MAF
    if not os.path.exists(FASTA38_PATH) or not os.path.exists(FASTA38_PATH + ".fai"):
        sys.exit(f"[ERROR] FASTA or .fai missing: {FASTA38_PATH}")
    fa = Fasta(FASTA38_PATH, rebuild=False, sequence_always_upper=True)

    maf_path = maf_path_for(args.subpop)
    if not os.path.exists(maf_path):
        sys.exit(f"[ERROR] MAF file for {args.subpop} not found: {maf_path}\n"
                 f"Hint: expected pattern maf_<POP>.grch38.tsv.gz under {MAF_DIR}")
    maf_lookup = load_maf_table(maf_path)
    maf_min, maf_max = args.ambig_maf_min, args.ambig_maf_max

    # Output naming: <PGSID>.<SUBPOP>.b38.tsv (or stem if PGSID missing)
    pgs_id  = meta.get("pgs_id", meta.get("pgsid","")).strip()
    trait   = meta.get("trait_reported", meta.get("trait_mapped","")).strip()
    weight_type = (meta.get("weight_type","beta") or "beta").lower()
    pgs_url = f"https://www.pgscatalog.org/score/{pgs_id}/" if pgs_id else ""
    stem_in = os.path.basename(args.inp).rsplit(".", 1)[0]
    stem    = pgs_id if pgs_id else stem_in
    os.makedirs(args.out_dir, exist_ok=True)
    out_tsv = os.path.join(args.out_dir, f"{stem}.{args.subpop}.b38.tsv")

    # Process rows
    out_rows = []
    for _, r in df.iterrows():
        chrom = to_ucsc(r[chr_c])
        try:
            pos = int(float(r[pos_c]))
        except ValueError:
            continue

        chrom, pos, strand = liftover_pos(lifter, chrom, pos)
        if chrom is None or pos is None:
            continue

        ea, oa = str(r[ea_c]).upper(), str(r[oa_c]).upper()
        if ea not in VALID_A or oa not in VALID_A or ea == oa:
            continue
        if strand == "-":
            ea, oa = ea.translate(COMPLEMENT), oa.translate(COMPLEMENT)

        ref = fetch_ref(fa, chrom, pos)
        if ref not in VALID_A:
            continue

        if ref == ea:
            alt, eff = oa, ea
        elif ref == oa:
            alt, eff = ea, ea
        else:
            continue

        # Palindromic handling
        pal = (set((ref, alt)) == {"A","T"}) or (set((ref, alt)) == {"C","G"})
        if pal:
            maf = maf_lookup.get((chrom, pos, ref, alt))
            if maf is None or (maf_min <= maf <= maf_max):
                continue

        # Effect size → beta
        try:
            w = float(r[wt_c])
        except ValueError:
            continue
        if weight_type in ("or","odds ratio","odds_ratio","hr","hazard ratio","hazard_ratio"):
            w = np.log(w)

        out_rows.append((chrom, pos, ref, alt, eff, w))

    if not out_rows:
        sys.exit("[ERROR] No variants survived harmonisation.")

    out_df = pd.DataFrame(out_rows, columns=["chr","pos","ref","alt","effect_allele","beta"])
    out_df.to_csv(out_tsv, sep="\t", index=False)

    meta_out = {
        "source_file": os.path.basename(args.inp),
        "n_input_rows": int(df.shape[0]),
        "n_harmonised_rows": int(out_df.shape[0]),
        "genome_build_src": build,
        "ancestry_tag": args.subpop,
        "maf_file": maf_path,
        "ambiguous_maf_window": [maf_min, maf_max],
        "weight_type_src": weight_type,
        "pgs_id": pgs_id,
        "trait": trait,
        "out_path": out_tsv,
    }
    with open(out_tsv + ".meta.json", "w") as fh:
        json.dump(meta_out, fh, indent=2)

    # Append to fixed registry
    append_registry({
        "timestamp_iso": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "pgs_id": pgs_id,
        "pgs_url": pgs_url,
        "trait": trait,
        "weight_type_src": weight_type,
        "genome_build_src": build,
        "ancestry": args.subpop,
        "maf_file": maf_path,
        "ambiguous_maf_window": f"{maf_min},{maf_max}",
        "source_file": os.path.basename(args.inp),
        "out_path": os.path.abspath(out_tsv),
        "n_input_rows": int(df.shape[0]),
        "n_harmonised_rows": int(out_df.shape[0]),
    })

    print(f"[OK] {stem} {args.subpop}: wrote {out_tsv}  (n={len(out_rows)})")

if __name__ == "__main__":
    main()
