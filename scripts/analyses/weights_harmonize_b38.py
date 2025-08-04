#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
weights_harmonize_b38.py

Convert PGS Catalog scoring files to a clean, GRCh38-harmonized format:
chr, pos, ref, alt, effect_allele, beta (+ metadata columns).

- Supports most PGS v2 files and many older variants.
- Liftover GRCh37→GRCh38 (requires chain file) and orients alleles to + strand.
- Ensures ref matches GRCh38 FASTA; drops sites that don't reconcile.
- Drops non-SNV and ambiguous palindromic SNPs (configurable).
- Converts OR/HR to log(OR/HR); keeps Beta/logOR as-is.

Usage:
  python weights_harmonize_b38.py \
      --in weights_raw/PGS000300.txt \
      --out weights_hm/ \
      --fasta38 /path/to/GRCh38.no_alt.fa \
      --chain37to38 /path/to/hg19ToHg38.over.chain.gz \
      --ancestry EUR \
      --drop-mhc

This will produce:
- weights_hm/PGS000300.harmonized.b38.tsv
- weights_hm/PGS000300.drops.tsv
- weights_hm/PGS000300.harmonized.meta.json

You can also pass a directory to --in; all *.txt/tsv files will be processed.

Output:
  weights_hm/<basename>.harmonized.b38.tsv
  weights_hm/<basename>.drops.tsv
  
  
What it does
- Ingests PGS Catalog scoring files (v2-style and many v1 variants).
- Parses metadata (e.g., pgs_id, weight_type, genome_build).
- Liftover GRCh37 → GRCh38 (if needed) using a chain file.
- Checks the GRCh38 reference base and orients alleles to chr:pos:ref:alt (forward strand).
- Converts OR/HR → log(OR/HR) so you always output beta.
- Drops anything non-biallelic SNP (indels, long alleles), ambiguous palindromic SNPs (configurable), and rows that don’t match the reference after liftover.
- Writes a harmonized TSV plus a per-file drop log with reasons.
- Dependencies: pandas, numpy, pyfaidx (for FASTA), and pyliftover
- Install: pip install pandas numpy pyfaidx pyliftover
"""

import argparse
import gzip
import io
import json
import os
import sys
from typing import Dict, Optional, Tuple, List

import numpy as np
import pandas as pd
from pyfaidx import Fasta
try:
    from pyliftover import LiftOver
except Exception as e:
    LiftOver = None


# ----------------------------- Helpers -------------------------------- #

CHR_REMAP = {
    "23": "X", "24": "Y", "M": "MT",
    "x": "X", "y": "Y", "mt": "MT", "m": "MT"
}

AUTOSOMES = {str(i) for i in range(1, 23)}
MHC_CHR = "chr6"
MHC_START_38 = 25_000_000
MHC_END_38   = 34_000_000

VALID_A = {"A", "C", "G", "T"}

COMPLEMENT = str.maketrans({"A":"T","T":"A","C":"G","G":"C",
                            "a":"t","t":"a","c":"g","g":"c"})

def read_text(path: str) -> str:
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")) as f:
        return f.read()

def iter_files(path: str) -> List[str]:
    if os.path.isdir(path):
        out = []
        for fn in os.listdir(path):
            if fn.lower().endswith((".txt", ".tsv", ".txt.gz", ".tsv.gz")):
                out.append(os.path.join(path, fn))
        return sorted(out)
    return [path]

def parse_metadata(raw_text: str) -> Dict[str, str]:
    meta = {}
    for line in raw_text.splitlines():
        if not line.startswith("#"):
            break
        line = line.strip()
        if line.startswith("##") or line.startswith("#"):
            if "=" in line:
                key = line.lstrip("#").split("=", 1)[0].strip()
                val = line.lstrip("#").split("=", 1)[1].strip()
                meta[key] = val
    # Normalize common keys
    meta["pgs_id"] = meta.get("pgs_id", meta.get("PGSID", meta.get("pgs", "")))
    meta["weight_type"] = meta.get("weight_type", "").strip()
    meta["genome_build"] = meta.get("genome_build", "").strip()
    meta["trait_reported"] = meta.get("trait_reported", meta.get("trait_mapped", ""))
    return meta

def load_table(raw_text: str) -> pd.DataFrame:
    # pandas can skip '#' comments, header is the first non-comment line
    buf = io.StringIO(raw_text)
    df = pd.read_csv(buf, sep="\t", comment="#", dtype=str)
    # strip whitespace in all string cells
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
    return df

def pick_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None

def is_snv(a1: str, a2: str) -> bool:
    return (a1 in VALID_A and a2 in VALID_A and len(a1) == 1 and len(a2) == 1 and a1 != a2)

def is_palindromic(a1: str, a2: str) -> bool:
    # A/T or T/A or C/G or G/C
    return (a1 == "A" and a2 == "T") or (a1 == "T" and a2 == "A") or \
           (a1 == "C" and a2 == "G") or (a1 == "G" and a2 == "C")

def with_chr_prefix(ch: str) -> str:
    ch0 = str(ch)
    ch0 = CHR_REMAP.get(ch0, ch0)
    if not ch0.startswith("chr"):
        ch0 = "chr" + ch0
    return ch0

def to_ucsc_chrom(ch: str) -> str:
    # For liftover library that expects 'chrN'
    return with_chr_prefix(ch)

def parse_int(x: str) -> Optional[int]:
    try:
        return int(float(x))
    except Exception:
        return None

def convert_weight_to_beta(w: float, weight_type: str) -> Optional[float]:
    wt = (weight_type or "").strip().lower()
    if wt in ("beta", "effect", "weight", "weights", "beta coefficient", "beta_coefficient"):
        return w
    # Sometimes 'LogOR', 'logOR', 'log odds ratio'
    if wt in ("logor", "log or", "lnor", "log_odds", "log_odds_ratio", "log-odds"):
        return w
    # Odds ratio / Hazard ratio -> log()
    if wt in ("or", "odds ratio", "odds_ratio", "hazard ratio", "hazard_ratio", "hr"):
        try:
            return float(np.log(w))
        except Exception:
            return None
    # Unknown -> assume beta
    return w

def fetch_ref_base(fa: Fasta, chrom: str, pos1: int) -> Optional[str]:
    try:
        base = fa[chrom][pos1-1:pos1].seq.upper()
        if base in VALID_A:
            return base
        return None
    except Exception:
        return None

def liftover_pos(lifter: LiftOver, chrom: str, pos1: int) -> Optional[Tuple[str, int, str]]:
    """
    Returns (chrom38, pos1_38, strand) or None if no mapping.
    pyliftover returns 0-based coords; we return 1-based.
    """
    try:
        results = lifter.convert_coordinate(chrom, pos1-1)
    except Exception:
        results = None
    if not results:
        return None
    # Pick the first mapping (best score)
    new_chrom, new_pos0, strand, _ = results[0]
    return new_chrom, int(new_pos0) + 1, strand

def orient_to_ref(ref_base: str, ea: str, oa: str) -> Optional[Tuple[str, str, str]]:
    """
    Given ref_base (A/C/G/T) and (effect_allele, other_allele) on + strand,
    return (ref, alt, effect_allele) if they can be reconciled. If not, None.
    """
    ea = ea.upper(); oa = oa.upper()
    if not is_snv(ea, oa):
        return None
    # Two possibilities: ref==ea or ref==oa
    if ref_base == ea:
        alt = oa
        return ref_base, alt, ea  # effect allele equals REF
    if ref_base == oa:
        alt = ea
        return ref_base, alt, ea  # effect allele equals ALT
    # Not matching reference
    return None

def drop_reason(row: dict, reason: str, extra: dict = None) -> dict:
    d = dict(row)
    d["drop_reason"] = reason
    if extra:
        d.update({f"ctx_{k}": v for k, v in extra.items()})
    return d


# ----------------------------- Main pipeline ------------------------------ #

def harmonize_file(
    in_path: str,
    out_dir: str,
    fasta38_path: str,
    chain37to38: Optional[str],
    ancestry: str,
    drop_palindromic: bool,
    drop_mhc: bool,
    include_chroms: str = "1-22,X,Y",
) -> Tuple[str, str]:
    """
    Returns (harmonized_tsv_path, drops_tsv_path)
    """
    os.makedirs(out_dir, exist_ok=True)
    raw_text = read_text(in_path)
    meta = parse_metadata(raw_text)

    df = load_table(raw_text)
    orig_cols = set(df.columns)

    # Column mapping (handle v1/v2 variants)
    col_chr = pick_column(df, ["chr_name", "chromosome", "chr", "hm_chr"])
    col_pos = pick_column(df, ["chr_position", "position", "pos", "hm_pos"])
    col_ea  = pick_column(df, ["effect_allele", "effectAllele", "EA", "hm_effect_allele"])
    col_oa  = pick_column(df, ["other_allele", "non_effect_allele", "otherAllele", "OA", "hm_other_allele"])
    col_w   = pick_column(df, ["effect_weight", "weight", "beta", "logOR", "log_or", "odds_ratio", "OR", "HR"])
    col_rsid= pick_column(df, ["rsID", "rsid", "variant_id", "rs_number"])

    required = [col_chr, col_pos, col_ea, col_oa, col_w]
    if any(c is None for c in required):
        missing = ["chr","pos","effect_allele","other_allele","effect_weight"]
        raise RuntimeError(
            f"{os.path.basename(in_path)}: missing required columns. "
            f"Found: {sorted(orig_cols)}. Need candidates for {missing}"
        )

    # Genome build
    build = (meta.get("genome_build") or "").upper()
    if build not in ("GRCH38", "HG38", "GRCH37", "HG19", ""):
        # sometimes it's blank; allow it
        pass

    # Prepare liftover if needed
    lifter = None
    needs_liftover = False
    if build in ("GRCH37", "HG19", "B37") or build == "":
        # If build unknown, we'll still allow liftover if provided (safer default).
        if chain37to38 is None:
            if build == "GRCH37" or build == "HG19":
                raise RuntimeError("This file is GRCh37 but --chain37to38 not provided.")
        else:
            if LiftOver is None:
                raise RuntimeError("pyliftover not available; install or omit GRCh37 files.")
            lifter = LiftOver(chain37to38)
            needs_liftover = True

    # Load FASTA
    fa = Fasta(fasta38_path, rebuild=False, sequence_always_upper=True)

    # Which chromosomes to keep
    allowed = set()
    for tok in include_chroms.split(","):
        tok = tok.strip()
        if "-" in tok:
            a, b = tok.split("-")
            try:
                a_i = 23 if a.upper() == "X" else 24 if a.upper() == "Y" else int(a)
                b_i = 23 if b.upper() == "X" else 24 if b.upper() == "Y" else int(b)
                for i in range(a_i, b_i + 1):
                    allowed.add("chr" + ("X" if i == 23 else "Y" if i == 24 else str(i)))
            except Exception:
                pass
        else:
            allowed.add(with_chr_prefix(tok))

    rows_out = []
    rows_drop = []

    weight_type = meta.get("weight_type", "Beta")

    for idx, r in df.iterrows():
        try:
            ch_raw = str(r[col_chr]).strip()
            pos_raw = r[col_pos]
            ea = str(r[col_ea]).upper().strip()
            oa = str(r[col_oa]).upper().strip()
            w  = r[col_w]
        except Exception:
            rows_drop.append(drop_reason(r, "missing_required_fields"))
            continue

        # Basic sanity
        pos = parse_int(pos_raw)
        if pos is None:
            rows_drop.append(drop_reason(r, "bad_position", {"pos": pos_raw}))
            continue

        # Normalize chrom tokens and allow 1..22, X, Y
        ch_norm = CHR_REMAP.get(ch_raw, ch_raw)
        ch_ucsc = with_chr_prefix(ch_norm)
        # Liftover if needed
        strand = "+"
        if needs_liftover:
            lo = liftover_pos(lifter, to_ucsc_chrom(ch_norm), pos)
            if lo is None:
                rows_drop.append(drop_reason(r, "liftover_failed"))
                continue
            ch_ucsc, pos, strand = lo
            if not ch_ucsc.startswith("chr"):
                ch_ucsc = "chr" + ch_ucsc.lstrip("chr")

        # Keep only requested chroms
        if allowed and ch_ucsc not in allowed:
            rows_drop.append(drop_reason(r, "chromosome_excluded", {"chrom": ch_ucsc}))
            continue

        # Only handle SNVs for pilot
        if not is_snv(ea, oa):
            rows_drop.append(drop_reason(r, "non_snv"))
            continue

        # If liftover landed on negative strand, complement alleles to + strand
        if strand == "-":
            ea = ea.translate(COMPLEMENT).upper()
            oa = oa.translate(COMPLEMENT).upper()

        # Optionally drop MHC region (GRCh38)
        if drop_mhc and ch_ucsc == MHC_CHR and MHC_START_38 <= pos <= MHC_END_38:
            rows_drop.append(drop_reason(r, "mhc_region"))
            continue

        # Drop ambiguous palindromic SNPs (A/T or C/G)
        if drop_palindromic and is_palindromic(ea, oa):
            rows_drop.append(drop_reason(r, "palindromic_ambiguous"))
            continue

        # Fetch reference base from GRCh38
        ref_base = fetch_ref_base(fa, ch_ucsc, pos)
        if ref_base is None or ref_base not in VALID_A:
            rows_drop.append(drop_reason(r, "ref_fetch_failed", {"chrom": ch_ucsc, "pos": pos}))
            continue

        # Orient to reference
        oriented = orient_to_ref(ref_base, ea, oa)
        if oriented is None:
            rows_drop.append(drop_reason(r, "ref_mismatch", {"ref38": ref_base, "ea": ea, "oa": oa}))
            continue
        ref38, alt38, ea38 = oriented  # ea38 equals ref or alt

        # Parse/convert weight to beta
        try:
            w_val = float(str(w))
        except Exception:
            rows_drop.append(drop_reason(r, "bad_weight", {"weight": w}))
            continue

        beta = convert_weight_to_beta(w_val, weight_type)
        if beta is None or not np.isfinite(beta):
            rows_drop.append(drop_reason(r, "weight_conversion_failed", {"weight": w, "type": weight_type}))
            continue

        # Assemble output row
        rows_out.append({
            "pgs_id": meta.get("pgs_id", ""),
            "trait": meta.get("trait_reported", ""),
            "weight_type_src": weight_type,
            "ancestry": ancestry,
            "chr": ch_ucsc,
            "pos": int(pos),
            "ref": ref38,
            "alt": alt38,
            "effect_allele": ea38,   # for scoring, use DS if effect==ALT else 2-DS
            "beta": float(beta),
            "rsid": r.get(col_rsid, "") if col_rsid else "",
        })

    # Write outputs
    base = os.path.basename(in_path)
    stem = base.replace(".txt.gz","").replace(".tsv.gz","").replace(".txt","").replace(".tsv","")
    out_path = os.path.join(out_dir, f"{stem}.harmonized.b38.tsv")
    drop_path = os.path.join(out_dir, f"{stem}.drops.tsv")

    if rows_out:
        out_df = pd.DataFrame(rows_out, columns=[
            "pgs_id","trait","ancestry","weight_type_src","chr","pos","ref","alt","effect_allele","beta","rsid"
        ])
        out_df.sort_values(["chr","pos"], inplace=True, key=lambda s: s.map(chrom_sort_key))
        out_df.to_csv(out_path, sep="\t", index=False)
    else:
        # still write an empty file with header
        pd.DataFrame(columns=[
            "pgs_id","trait","ancestry","weight_type_src","chr","pos","ref","alt","effect_allele","beta","rsid"
        ]).to_csv(out_path, sep="\t", index=False)

    if rows_drop:
        dd = pd.DataFrame(rows_drop)
        # Keep original columns if they existed
        dd.to_csv(drop_path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["drop_reason"]).to_csv(drop_path, sep="\t", index=False)

    # Write a tiny meta sidecar for auditability
    meta_out = {
        "source_file": base,
        "pgs_id": meta.get("pgs_id",""),
        "weight_type_src": weight_type,
        "genome_build_src": meta.get("genome_build",""),
        "n_input_rows": int(df.shape[0]),
        "n_harmonized_rows": int(len(rows_out)),
        "n_dropped_rows": int(len(rows_drop)),
    }
    with open(os.path.join(out_dir, f"{stem}.harmonized.meta.json"), "w") as f:
        json.dump(meta_out, f, indent=2)

    return out_path, drop_path


def chrom_sort_key(series: pd.Series) -> pd.Series:
    """
    Map chr tokens to sortable integers.
    """
    def one(x: str) -> int:
        x = str(x)
        if not x.startswith("chr"):
            x = "chr" + x
        tok = x[3:]
        if tok == "X": return 23
        if tok == "Y": return 24
        if tok in ("M","MT"): return 25
        try:
            return int(tok)
        except Exception:
            return 100
    return series.map(one)


def main():
    p = argparse.ArgumentParser(description="Harmonize PGS weights to GRCh38 chr:pos:ref:alt with effect-allele betas.")
    p.add_argument("--in", dest="inp", required=True, help="Input PGS scoring file (.txt/.tsv[.gz]) or a directory")
    p.add_argument("--out", dest="out_dir", required=True, help="Output directory")
    p.add_argument("--fasta38", required=True, help="GRCh38 FASTA (indexed; e.g., GRCh38.no_alt.fa)")
    p.add_argument("--chain37to38", default=None, help="hg19ToHg38.over.chain(.gz) for liftover. Required if source is GRCh37.")
    p.add_argument("--ancestry", default="UNK", help="Ancestry tag to embed in output (e.g., EUR, AFR, EAS).")
    p.add_argument("--include-chroms", default="1-22,X", help="Chromosomes to include, default autosomes+X.")
    p.add_argument("--drop-palindromic", action="store_true", help="Drop all palindromic A/T, C/G SNPs (default).")
    p.add_argument("--keep-palindromic", dest="drop_palindromic", action="store_false", help="Keep palindromic SNPs (not recommended without MAF gating).")
    p.add_argument("--drop-mhc", action="store_true", help="Drop MHC region (chr6:25-34Mb, GRCh38).")
    args = p.parse_args()

    files = iter_files(args.inp)
    if not files:
        sys.exit("No input files found.")

    for f in files:
        print(f"[+] Harmonizing {os.path.basename(f)}")
        try:
            harmonize_file(
                f, args.out_dir, args.fasta38, args.chain37to38,
                ancestry=args.ancestry,
                drop_palindromic=args.drop_palindromic if hasattr(args, "drop_palindromic") else True,
                drop_mhc=args.drop_mhc,
                include_chroms=args.include_chroms
            )
        except Exception as e:
            print(f"[!] Failed on {f}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()

"""
Optional next tweaks (if you want them)
MAF-aware palindromic keep: add --maf /path/to/1kg.maf.tsv to only drop palindromics with MAF ~0.5.

Ambiguous position resolution: prefer mappings where the lifted reference base matches either allele; the current “first liftover result” works well in practice with UCSC chains, but you can add a tie-break if needed.

Include-only autosomes: change --include-chroms (default 1-22,X) to just 1-22.
"""