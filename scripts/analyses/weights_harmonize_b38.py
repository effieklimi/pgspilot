#!/usr/bin/env python3
"""
pgs_harmonize_mvp.py  ── minimum-viable PGS weight harmoniser (GRCh38, MAF-aware)

• Accepts ONE PGS Catalog scoring file (.txt/.tsv[.gz]) that may be on
  GRCh37/hg19 or GRCh38.
• Emits a GRCh38-based weight file with six columns:
    chr  pos  ref  alt  effect_allele  beta
• Liftover 37→38 if needed (requires UCSC hg19ToHg38 chain file).
• Ensures REF/ALT match the GRCh38 reference FASTA; drops any rows that
  fail to reconcile or that are not biallelic SNVs.
• Converts OR / HR to ln(OR/HR); leaves beta or logOR untouched.
• MAF-aware handling of palindromic SNPs (A/T or C/G):
    - If --maf-file provided: keep palindromes only when MAF is outside
      [--ambig-maf-min, --ambig-maf-max] (default 0.45–0.55).
      If MAF is missing → drop.
    - If --maf-file not provided: drop all palindromic SNPs (safe default).

Dependencies (pure-Python):
    pip install pandas numpy pyfaidx pyliftover

Example:
    python pgs_harmonize_mvp.py \
        --in PGS000300.txt.gz \
        --out PGS000300.b38.tsv \
        --fasta38 /ref/GRCh38.no_alt.fa \
        --chain37to38 /ref/hg19ToHg38.over.chain.gz \
        --maf-file /ref/maf_EUR.grch38.tsv.gz \
        --ancestry EUR

Writes:
  - PGS000300.b38.tsv
  - PGS000300.b38.tsv.meta.json
"""

import argparse, gzip, io, os, sys, json
from typing import Tuple, Optional, Dict, Tuple as TupleType

import numpy as np
import pandas as pd
from pyfaidx import Fasta

# pyliftover is required only if any input is on GRCh37/HG19
try:
    from pyliftover import LiftOver
except ImportError:
    LiftOver = None

VALID_A = {"A", "C", "G", "T"}
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

# ───────────────────────────── helpers ──────────────────────────────── #

def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def read_text(path: str) -> str:
    """Return full text, transparently handling .gz input."""
    with open_text(path) as f:
        return f.read()

def parse_meta(raw: str) -> dict:
    """Parse #key=value header lines into a dict (lower-cased keys)."""
    meta = {}
    for line in raw.splitlines():
        if not line.startswith("#"):
            break
        if "=" in line:
            k, v = line.lstrip("#").split("=", 1)
            meta[k.strip().lower()] = v.strip()
    return meta

def to_ucsc(ch: str) -> str:
    chs = str(ch)
    return chs if chs.startswith("chr") else f"chr{chs}"

def liftover_pos(
    lifter: Optional["LiftOver"], chrom: str, pos1: int
) -> Tuple[Optional[str], Optional[int], str]:
    """Return (chrom38, pos1_38, strand). If no lifter, return input."""
    if lifter is None:
        return chrom, pos1, "+"
    res = lifter.convert_coordinate(chrom, pos1 - 1)
    if not res:
        return None, None, "+"
    new_ch, new_pos0, strand, _ = res[0]
    if not str(new_ch).startswith("chr"):
        new_ch = "chr" + str(new_ch)
    return new_ch, int(new_pos0) + 1, strand

def fetch_ref(fa: Fasta, chrom: str, pos1: int) -> Optional[str]:
    try:
        return fa[chrom][pos1 - 1 : pos1].seq.upper()
    except Exception:
        return None

def load_maf_table(path: Optional[str]) -> Dict[TupleType[str, int, str, str], float]:
    """
    Load TSV[.gz] with columns (case-insensitive):
        chr, pos, ref, alt, maf   OR   chr, pos, ref, alt, af
    We store maf = min(af, 1-af). Keys must be GRCh38.
    Dict is keyed by (chr,pos,ref,alt) and also (chr,pos,alt,ref) to be order-agnostic.
    """
    if not path:
        return {}
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
        raise SystemExit("[ERROR] --maf-file missing required columns (need chr,pos,ref,alt and maf/af).")
    out: Dict[TupleType[str,int,str,str], float] = {}
    for _, r in df.iterrows():
        try:
            chrom = r[c_chr]
            chrom = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
            pos = int(float(r[c_pos]))
            ref = str(r[c_ref]).upper()
            alt = str(r[c_alt]).upper()
            maf = float(r[c_maf]) if c_maf and r[c_maf] != "" else float(r[c_af])
            maf = maf if maf <= 0.5 else 1.0 - maf
        except Exception:
            continue
        if ref in VALID_A and alt in VALID_A and ref != alt:
            out[(chrom,pos,ref,alt)] = maf
            out[(chrom,pos,alt,ref)] = maf
    return out

# ───────────────────────────── main ─────────────────────────────────── #

def main():
    ap = argparse.ArgumentParser("PGS MVP harmoniser → GRCh38 (MAF-aware palindromic handling)")
    ap.add_argument("--in", dest="inp", required=True, help="PGS scoring file (.txt/.tsv[.gz])")
    ap.add_argument("--out", dest="out", required=True, help="Output .tsv path")
    ap.add_argument("--fasta38", required=True, help="GRCh38 reference FASTA (indexed .fai)")
    ap.add_argument("--chain37to38", help="hg19ToHg38.over.chain(.gz) for liftover (required if input is GRCh37/HG19)")
    ap.add_argument("--maf-file", help="Per-ancestry MAF TSV[.gz] with chr,pos,ref,alt and maf/af on GRCh38 (for palindrome handling)")
    ap.add_argument("--ambig-maf-min", type=float, default=0.45, help="Lower bound for ambiguous palindrome drop window")
    ap.add_argument("--ambig-maf-max", type=float, default=0.55, help="Upper bound for ambiguous palindrome drop window")
    ap.add_argument("--ancestry", default="UNK", help="Optional ancestry tag for bookkeeping (e.g., EUR/AFR/EAS)")
    args = ap.parse_args()

    raw = read_text(args.inp)
    meta = parse_meta(raw)
    df = pd.read_csv(io.StringIO(raw), sep="\t", comment="#", dtype=str).fillna("")

    # Column mapping (v2 and most v1 files)
    COLS = {
        "chr": ["chr_name", "chromosome", "chr"],
        "pos": ["chr_position", "position", "pos"],
        "ea": ["effect_allele", "EA"],
        "oa": ["other_allele", "OA"],
        "wt": ["effect_weight", "beta", "OR", "HR", "logOR", "log_or", "odds_ratio"],
    }
    def pick(cands): return next((c for c in cands if c in df.columns), None)

    chr_c, pos_c, ea_c, oa_c, wt_c = [pick(COLS[k]) for k in ("chr", "pos", "ea", "oa", "wt")]
    if None in (chr_c, pos_c, ea_c, oa_c, wt_c):
        sys.exit("[ERROR] Scoring file lacks required columns – aborting.")

    build = (meta.get("genome_build","") or "").upper()
    lifter = None
    if build in ("GRCH37","HG19","B37") or args.chain37to38:
        if LiftOver is None:
            sys.exit("pyliftover not installed – cannot liftover GRCh37 input.")
        if not args.chain37to38 and build in ("GRCH37","HG19","B37"):
            sys.exit("[ERROR] --chain37to38 is required for GRCh37/HG19 input.")
        if args.chain37to38:
            lifter = LiftOver(args.chain37to38)

    fa = Fasta(args.fasta38, rebuild=False, sequence_always_upper=True)
    maf_lookup = load_maf_table(args.maf_file)
    maf_min, maf_max = args.ambig_maf_min, args.ambig_maf_max

    weight_type = (meta.get("weight_type","beta") or "beta").lower()
    out_rows = []

    for _, row in df.iterrows():
        chrom = to_ucsc(row[chr_c])
        try:
            pos = int(float(row[pos_c]))
        except ValueError:
            continue

        chrom, pos, strand = liftover_pos(lifter, chrom, pos)
        if chrom is None or pos is None:
            continue

        ea, oa = row[ea_c].upper(), row[oa_c].upper()
        if ea not in VALID_A or oa not in VALID_A or ea == oa:
            continue
        if strand == "-":
            ea, oa = ea.translate(COMPLEMENT), oa.translate(COMPLEMENT)

        ref = fetch_ref(fa, chrom, pos)
        if ref not in VALID_A:
            continue

        # Orient alleles to reference
        if ref == ea:
            alt, eff = oa, ea
        elif ref == oa:
            alt, eff = ea, ea
        else:
            # neither EA nor OA matches reference base
            continue

        # Drop non-SNV (already ensured A/C/G/T) and handle palindromes via MAF
        pal = (set((ref, alt)) == {"A","T"}) or (set((ref, alt)) == {"C","G"})
        if pal:
            if maf_lookup:
                maf = maf_lookup.get((chrom, pos, ref, alt))
                if maf is None or (maf_min <= maf <= maf_max):
                    continue
            else:
                # No MAF file supplied → drop all palindromes (safe default)
                continue

        # Effect size → beta
        try:
            w = float(row[wt_c])
        except ValueError:
            continue
        wt = weight_type
        if wt in ("or","odds ratio","odds_ratio","hr","hazard ratio","hazard_ratio"):
            w = np.log(w)  # ln(OR/HR)
        # keep beta/logOR as-is

        out_rows.append((chrom, pos, ref, alt, eff, w))

    if not out_rows:
        sys.exit("[ERROR] No variants survived harmonisation.")

    out_df = pd.DataFrame(out_rows, columns=["chr","pos","ref","alt","effect_allele","beta"])
    out_df.to_csv(args.out, sep="\t", index=False)

    with open(args.out + ".meta.json", "w") as fh:
        json.dump(
            {
                "source_file": os.path.basename(args.inp),
                "n_input_rows": int(df.shape[0]),
                "n_harmonised_rows": int(out_df.shape[0]),
                "genome_build_src": build,
                "ancestry_tag": args.ancestry,
                "maf_file": args.maf_file or "",
                "ambiguous_maf_window": [maf_min, maf_max],
                "weight_type_src": weight_type,
            },
            fh,
            indent=2,
        )

if __name__ == "__main__":
    main()
