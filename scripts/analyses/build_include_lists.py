#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_include_lists.py

Purpose:
  Freeze the fixed variant list(s) you'll score for each trait × ancestry,
  by intersecting harmonized weights with scorable sites (+ optional site mask).

Inputs:
  --weights <file_or_dir>            # harmonized weight TSVs from weights_harmonize_b38.py
  --scorable <scorable_sites.b38.tsv>
  --site-mask <site_mask.tsv>        # optional, presence means 'keep'; can also include avgDR2
  --mask-min-dr2 <float>             # optional threshold if site-mask has a column like avgDR2
  --include-chroms "1-22,X"          # default: autosomes + X
  --drop-mhc                         # optional: drop MHC (chr6:25-34Mb, GRCh38)
  --union-by-trait                   # optional: emit one include list per trait (union across ancestries)
  --out <traits_dir>                 # default: ./traits

Outputs (per trait × ancestry OR per trait if union):
  traits/<trait>__<ANC>.include.tsv                  # columns: chr pos ref alt effect_allele
  traits/<trait>__<ANC>.include.meta.json           # summary counts & provenance
  traits/<trait>__<ANC>.include.drops.tsv           # rows excluded and reasons

Notes:
  - The include list is what prevents per-user variant drift:
    all users get scored against the same frozen keys.
  - If you use --union-by-trait, effect_allele conflicts across ancestries are
    resolved by simple majority rule (ties → 'ALT'); logged in the drops file.
"""

import argparse
import os
import sys
import json
import gzip
import io
import re
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

MHC_CHR = "chr6"
MHC_START_38 = 25_000_000
MHC_END_38   = 34_000_000

def with_chr_prefix(ch: str) -> str:
    ch = str(ch).strip()
    remap = {"23": "X", "24": "Y", "M": "MT", "x":"X", "y":"Y", "mt":"MT", "m":"MT"}
    ch = remap.get(ch, ch)
    return ch if ch.startswith("chr") else f"chr{ch}"

def chrom_sort_key(series: pd.Series) -> pd.Series:
    def one(x: str) -> int:
        x = str(x)
        if not x.startswith("chr"):
            x = "chr" + x
        tok = x[3:]
        if tok == "X": return 23
        if tok == "Y": return 24
        if tok in ("M", "MT"): return 25
        try:
            return int(tok)
        except Exception:
            return 100
    return series.map(one)

def parse_tsv(path: str) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        df = pd.read_csv(f, sep="\t", dtype=str)
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
    return df

def load_scorable_sites(path: str, include_chroms: List[str]) -> pd.DataFrame:
    df = parse_tsv(path)
    # Flexible column names
    cols = {c.lower(): c for c in df.columns}
    def need(*names):
        for n in names:
            if n in cols:
                return cols[n]
        raise RuntimeError(f"{path}: required column missing; looked for one of {names}, found {df.columns.tolist()}")
    c_chrom = need("chr", "chrom", "chromosome")
    c_pos   = need("pos", "position")
    c_ref   = need("ref", "reference", "ref_allele", "reference_allele")
    c_alt   = need("alt", "alternate", "alt_allele", "alternate_allele")

    df = df[[c_chrom, c_pos, c_ref, c_alt]].rename(columns={
        c_chrom: "chr", c_pos: "pos", c_ref: "ref", c_alt: "alt"
    })
    df["chr"] = df["chr"].apply(with_chr_prefix)
    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()

    if include_chroms:
        allowed = set(include_chroms)
        df = df[df["chr"].isin(allowed)].copy()

    # Keep only biallelic SNVs (single base A/C/G/T)
    mask_snv = df["ref"].str.len().eq(1) & df["alt"].str.len().eq(1) & (df["ref"] != df["alt"])
    df = df[mask_snv].copy()

    return df

def load_site_mask(path: str, include_chroms: List[str],
                   mask_min_dr2: Optional[float]) -> pd.DataFrame:
    df = parse_tsv(path)
    cols = {c.lower(): c for c in df.columns}
    def pick(*names):
        for n in names:
            if n in cols:
                return cols[n]
        return None
    c_chrom = pick("chr", "chrom", "chromosome")
    c_pos   = pick("pos", "position")
    c_ref   = pick("ref", "reference", "ref_allele", "reference_allele")
    c_alt   = pick("alt", "alternate", "alt_allele", "alternate_allele")
    if not all([c_chrom, c_pos, c_ref, c_alt]):
        raise RuntimeError(f"{path}: site mask must include chr,pos,ref,alt columns")

    df = df[[c_chrom, c_pos, c_ref, c_alt] + [c for c in df.columns if c not in (c_chrom,c_pos,c_ref,c_alt)]]
    df = df.rename(columns={c_chrom:"chr", c_pos:"pos", c_ref:"ref", c_alt:"alt"})
    df["chr"] = df["chr"].apply(with_chr_prefix)
    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()

    if include_chroms:
        allowed = set(include_chroms)
        df = df[df["chr"].isin(allowed)].copy()

    # If a DR2-ish column exists and threshold is provided, filter
    dr2_col = None
    for cand in ["avgdr2", "mean_dr2", "dr2", "avg_info", "info", "avg_rsqr", "rsq"]:
        if cand in (c.lower() for c in df.columns):
            for c in df.columns:
                if c.lower() == cand:
                    dr2_col = c
                    break
            if dr2_col:
                break

    if dr2_col and (mask_min_dr2 is not None):
        # non-numeric become NaN → dropped
        df[dr2_col] = pd.to_numeric(df[dr2_col], errors="coerce")
        before = df.shape[0]
        df = df[df[dr2_col] >= mask_min_dr2].copy()
        after = df.shape[0]
        print(f"[site-mask] Applied {dr2_col}>={mask_min_dr2} → kept {after}/{before}")

    # Keep only SNVs
    mask_snv = df["ref"].str.len().eq(1) & df["alt"].str.len().eq(1) & (df["ref"] != df["alt"])
    df = df[mask_snv].copy()

    return df[["chr","pos","ref","alt"]].drop_duplicates()

def list_weight_files(path: str) -> List[str]:
    if os.path.isdir(path):
        outs = []
        for fn in os.listdir(path):
            if fn.endswith(".tsv") or fn.endswith(".tsv.gz"):
                outs.append(os.path.join(path, fn))
        return sorted(outs)
    return [path]

def sanitize(s: str) -> str:
    s = s.strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.+-]", "", s)
    return s or "UNK"

def in_mhc(chrom: str, pos: int) -> bool:
    return chrom == MHC_CHR and (MHC_START_38 <= pos <= MHC_END_38)

def parse_include_chroms(arg: str) -> List[str]:
    out = set()
    for tok in arg.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok:
            a, b = tok.split("-")
            def to_num(x):
                x = x.strip().upper()
                return 23 if x == "X" else 24 if x == "Y" else int(x)
            ai, bi = to_num(a), to_num(b)
            for i in range(ai, bi + 1):
                label = "X" if i == 23 else "Y" if i == 24 else str(i)
                out.add(with_chr_prefix(label))
        else:
            out.add(with_chr_prefix(tok))
    return sorted(out, key=lambda c: int(c[3:]) if c[3:].isdigit() else (23 if c.endswith("X") else 24 if c.endswith("Y") else 100))

def load_weights(path: str) -> pd.DataFrame:
    df = parse_tsv(path)
    # Expect standard columns from harmonizer
    required = ["chr","pos","ref","alt","effect_allele"]
    for req in required:
        if req not in df.columns:
            raise RuntimeError(f"{path}: missing required column '{req}'. Columns present: {df.columns.tolist()}")
    # Optional metadata
    for c in ["trait","ancestry","pgs_id"]:
        if c not in df.columns:
            df[c] = "UNK"
    # Normalize types/formats
    df["chr"] = df["chr"].apply(with_chr_prefix)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    df["effect_allele"] = df["effect_allele"].str.upper()
    # Only biallelic SNVs
    mask_snv = df["ref"].str.len().eq(1) & df["alt"].str.len().eq(1) & (df["ref"] != df["alt"])
    df = df[mask_snv].copy()
    # Drop rows with missing pos
    df = df[df["pos"].notna()].copy()
    df["pos"] = df["pos"].astype(int)
    return df

def write_include(out_dir: str, trait: str, anc: str,
                  df_keep: pd.DataFrame,
                  drops: pd.DataFrame,
                  sources: List[str]):
    os.makedirs(out_dir, exist_ok=True)
    trait_tag = sanitize(trait)
    anc_tag = sanitize(anc)
    base = f"{trait_tag}__{anc_tag}" if anc_tag != "ALL" else f"{trait_tag}"
    out_tsv = os.path.join(out_dir, f"{base}.include.tsv")
    out_meta = os.path.join(out_dir, f"{base}.include.meta.json")
    out_drop = os.path.join(out_dir, f"{base}.include.drops.tsv")

    # Sort & output include list
    if not df_keep.empty:
        df_keep = df_keep[["chr","pos","ref","alt","effect_allele"]].drop_duplicates()
        df_keep.sort_values(["chr","pos"], inplace=True, key=chrom_sort_key)
        df_keep.to_csv(out_tsv, sep="\t", index=False)
    else:
        # still write header
        pd.DataFrame(columns=["chr","pos","ref","alt","effect_allele"]).to_csv(out_tsv, sep="\t", index=False)

    # Drops
    if drops is not None and not drops.empty:
        drops.to_csv(out_drop, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["chr","pos","ref","alt","reason","source_file"]).to_csv(out_drop, sep="\t", index=False)

    # Meta
    meta = {
        "trait": trait,
        "ancestry": anc,
        "n_included": int(df_keep.shape[0]),
        "n_dropped": int(0 if drops is None else drops.shape[0]),
        "sources": [os.path.basename(s) for s in sources]
    }
    with open(out_meta, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"[✓] Wrote {out_tsv}  (n={meta['n_included']})")

def main():
    ap = argparse.ArgumentParser(description="Build frozen include lists per trait×ancestry from harmonized weights, scorable sites, and optional site mask.")
    ap.add_argument("--weights", required=True, help="Harmonized weights file or directory")
    ap.add_argument("--scorable", required=True, help="scorable_sites.b38.tsv")
    ap.add_argument("--site-mask", default=None, help="Optional site_mask.tsv (presence means keep). May include avgDR2.")
    ap.add_argument("--mask-min-dr2", type=float, default=None, help="If site-mask has a DR2/INFO column, require value ≥ this threshold.")
    ap.add_argument("--include-chroms", default="1-22,X", help='Chromosomes to include, e.g., "1-22,X" (default).')
    ap.add_argument("--drop-mhc", action="store_true", help="Drop MHC region (chr6:25-34Mb, GRCh38).")
    ap.add_argument("--union-by-trait", action="store_true", help="Emit a single include list per trait (union across ancestries).")
    ap.add_argument("--out", default="traits", help="Output directory (default: ./traits)")
    args = ap.parse_args()

    include_chroms = parse_include_chroms(args.include_chroms)

    # Load scorable sites
    print("[*] Loading scorable sites…")
    scorable = load_scorable_sites(args.scorable, include_chroms)
    scorable["key"] = scorable["chr"] + ":" + scorable["pos"].astype(str) + ":" + scorable["ref"] + ":" + scorable["alt"]
    scorable_keys = set(scorable["key"].tolist())
    print(f"    scorable sites: {len(scorable_keys):,}")

    # Load site mask (optional)
    mask_keys = None
    if args.site_mask:
        print("[*] Loading site mask…")
        mask_df = load_site_mask(args.site_mask, include_chroms, args.mask_min_dr2)
        mask_df["key"] = mask_df["chr"] + ":" + mask_df["pos"].astype(str) + ":" + mask_df["ref"] + ":" + mask_df["alt"]
        mask_keys = set(mask_df["key"].tolist())
        print(f"    mask sites: {len(mask_keys):,}")

    # Read harmonized weights
    files = list_weight_files(args.weights)
    if not files:
        sys.exit("No weight files found.")

    print(f"[*] Reading {len(files)} harmonized weight file(s)…")
    all_rows = []
    file_of_row = []  # to track provenance for drops
    for f in files:
        df = load_weights(f)
        # Keep only requested chroms
        if include_chroms:
            df = df[df["chr"].isin(include_chroms)].copy()

        # Drop MHC region if requested
        if args.drop_mhc:
            mask_not_mhc = ~( (df["chr"] == MHC_CHR) & (df["pos"].between(MHC_START_38, MHC_END_38)) )
            df = df[mask_not_mhc].copy()

        # Build key for intersection
        df["key"] = df["chr"] + ":" + df["pos"].astype(str) + ":" + df["ref"] + ":" + df["alt"]
        all_rows.append(df)
        file_of_row.extend([f] * df.shape[0])

    if not all_rows:
        sys.exit("No rows found in weight files.")

    W = pd.concat(all_rows, ignore_index=True)
    W["source_file"] = file_of_row

    # Intersect with scorable sites
    in_scorable = W["key"].isin(scorable_keys)
    # Intersect with site mask if provided
    if mask_keys is not None:
        in_mask = W["key"].isin(mask_keys)
    else:
        in_mask = pd.Series([True] * W.shape[0])

    keep_mask = in_scorable & in_mask
    drops = W.loc[~keep_mask, ["chr","pos","ref","alt","trait","ancestry","source_file"]].copy()
    drops["reason"] = drops.apply(
        lambda r: "not_in_scorable" if r["chr"]+":"+str(r["pos"])+":"+r["ref"]+":"+r["alt"] not in scorable_keys
                  else "not_in_mask",
        axis=1
    )
    keep = W.loc[keep_mask, ["trait","ancestry","chr","pos","ref","alt","effect_allele","source_file","key"]].copy()

    # Deduplicate within trait×ancestry on key; if multiple effect_alleles conflict, resolve by majority
    def resolve_effect_allele(df_sub: pd.DataFrame) -> pd.DataFrame:
        # df_sub: rows for a single (trait, ancestry, key) possibly from multiple files
        ea_counts = df_sub["effect_allele"].value_counts()
        if ea_counts.shape[0] == 1:
            ea = ea_counts.index[0]
            src = ";".join(sorted(set(df_sub["source_file"].tolist())))
            return pd.DataFrame({
                "chr": [df_sub["chr"].iloc[0]],
                "pos": [df_sub["pos"].iloc[0]],
                "ref": [df_sub["ref"].iloc[0]],
                "alt": [df_sub["alt"].iloc[0]],
                "effect_allele": [ea],
                "source_files": [src],
                "conflict": [False]
            })
        # Conflict → choose majority; tie → prefer ALT (deterministic), else REF
        alt = df_sub["alt"].iloc[0]
        ref = df_sub["ref"].iloc[0]
        ea_major = ea_counts.idxmax()
        # If tie (freqs equal), prefer ALT if present
        if ea_counts.nunique() == 1 and len(ea_counts) == 2 and ea_counts.iloc[0] == ea_counts.iloc[1]:
            ea_major = alt if alt in ea_counts.index else (ref if ref in ea_counts.index else ea_major)
        src = ";".join(sorted(set(df_sub["source_file"].tolist())))
        out = pd.DataFrame({
            "chr": [df_sub["chr"].iloc[0]],
            "pos": [df_sub["pos"].iloc[0]],
            "ref": [ref],
            "alt": [alt],
            "effect_allele": [ea_major],
            "source_files": [src],
            "conflict": [True]
        })
        return out

    # Build outputs
    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    if args.union_by_trait:
        # One include per trait; union across ancestries
        traits = sorted(keep["trait"].unique())
        for t in traits:
            K = keep[keep["trait"] == t].copy()
            if K.empty:
                continue
            # Resolve by key (ignoring ancestry)
            grouped = K.groupby(["key"], sort=False)
            rows = []
            confs = []
            for key, df_sub in grouped:
                resolved = resolve_effect_allele(df_sub)
                resolved["trait"] = t
                resolved["ancestry"] = "ALL"
                rows.append(resolved)
                confs.append(resolved[["chr","pos","ref","alt","effect_allele","conflict","source_files"]])
            inc = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["chr","pos","ref","alt","effect_allele"])
            # Prepare drops file: include both intersection drops and allele conflicts (flagged)
            # Mark conflicts separately for transparency
            conf_df = pd.concat(confs, ignore_index=True) if confs else pd.DataFrame(columns=["chr","pos","ref","alt","effect_allele","conflict","source_files"])
            conflict_rows = conf_df[conf_df["conflict"]].copy()
            conflict_rows = conflict_rows.rename(columns={"source_files":"source_file"})
            conflict_rows["reason"] = "effect_allele_conflict_across_inputs"

            drops_t = drops[drops["trait"] == t].copy()
            drops_t = drops_t[["chr","pos","ref","alt","reason","source_file"]]

            all_drops = pd.concat([drops_t, conflict_rows[["chr","pos","ref","alt","reason","source_file"]]], ignore_index=True) \
                           .drop_duplicates()

            write_include(out_dir, trait=t, anc="ALL",
                          df_keep=inc[["chr","pos","ref","alt","effect_allele"]],
                          drops=all_drops,
                          sources=sorted(set(keep[keep["trait"] == t]["source_file"].tolist())))
    else:
        # One include per trait × ancestry
        pairs = keep[["trait","ancestry"]].drop_duplicates().sort_values(["trait","ancestry"])
        for _, row in pairs.iterrows():
            t, a = row["trait"], row["ancestry"]
            K = keep[(keep["trait"] == t) & (keep["ancestry"] == a)].copy()
            if K.empty:
                continue
            grouped = K.groupby(["key"], sort=False)
            rows = []
            confs = []
            for key, df_sub in grouped:
                resolved = resolve_effect_allele(df_sub)
                resolved["trait"] = t
                resolved["ancestry"] = a
                rows.append(resolved)
                confs.append(resolved[["chr","pos","ref","alt","effect_allele","conflict","source_files"]])
            inc = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["chr","pos","ref","alt","effect_allele"])
            conf_df = pd.concat(confs, ignore_index=True) if confs else pd.DataFrame(columns=["chr","pos","ref","alt","effect_allele","conflict","source_files"])
            conflict_rows = conf_df[conf_df["conflict"]].copy()
            conflict_rows = conflict_rows.rename(columns={"source_files":"source_file"})
            conflict_rows["reason"] = "effect_allele_conflict_across_inputs"

            drops_ta = drops[(drops["trait"] == t) & (drops["ancestry"] == a)].copy()
            drops_ta = drops_ta[["chr","pos","ref","alt","reason","source_file"]]

            all_drops = pd.concat([drops_ta, conflict_rows[["chr","pos","ref","alt","reason","source_file"]]], ignore_index=True) \
                           .drop_duplicates()

            write_include(out_dir, trait=t, anc=a,
                          df_keep=inc[["chr","pos","ref","alt","effect_allele"]],
                          drops=all_drops,
                          sources=sorted(set(K["source_file"].tolist())))

if __name__ == "__main__":
    main()
