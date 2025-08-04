#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_calibration.py  — MVP

Purpose:
  Create reference distributions F_k for every trait × ancestry you plan to report.
  For each trait and each ancestry with weights, compute PRS S_k on 1000G,
  then save ECDF (quantiles) and Normal params (mu, sigma).

Inputs:
  --weights_dir  weights_hm/                       # harmonized weights (from weights_harmonize_b38.py)
  --traits_dir   traits/                           # include lists (*.include.tsv)
  --vcf_pattern  "/path/1KG/ALL.chr{chr}.vcf.gz"   # {chr} replaced by 1..22 (we'll try with/without 'chr')
  --labels       /path/1KG_labels.tsv              # columns: sample, super_pop (AFR/AMR/EAS/EUR/SAS)

Optional:
  --out          calibration/                      # output base dir
  --anc          EUR,AFR,...                       # restrict to these ancestries (comma-separated)
  --quantiles    1001                              # number of ECDF points between 0..1
  --min-variants 20                                # skip if fewer variants end up used
  --include-chroms "1-22"                          # default autosomes only

Outputs (per trait in {out}/<trait>/):
  calib__<ANC>.npz     # keys: mu, sigma, q_probs, q_values
  meta.json            # summary of counts & exclusions per trait and ancestry

Notes:
  - Uses GT (hard calls) from 1KG VCFs; DS not required.
  - Includes only exact allele matches (chr:pos:ref:alt).
  - For duplicate weight rows on the same key, the first is kept; the rest are logged/dropped.
"""

import argparse
import os
import sys
import json
import re
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
from cyvcf2 import VCF

VALID_SUPERPOPS = {"AFR", "AMR", "EAS", "EUR", "SAS"}
AUTOSOMES = [str(i) for i in range(1, 23)]
VALID_A = {"A","C","G","T"}

# -------------------------- small utils -------------------------- #

def with_chr_prefix(ch: str) -> str:
    ch = str(ch).strip()
    return ch if ch.startswith("chr") else f"chr{ch}"

def chrom_key(ch: str) -> int:
    ch = with_chr_prefix(ch)
    try:
        return int(ch[3:])
    except Exception:
        return 99

def parse_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)

def list_files(dir_path: str, suffix: str) -> List[str]:
    return sorted([os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(suffix)])

def load_labels(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    sample_col = "sample" if "sample" in df.columns else ("iid" if "iid" in df.columns else None)
    if sample_col is None or "super_pop" not in df.columns:
        raise RuntimeError("labels file must have columns: sample (or IID), superpop")
    df = df[[sample_col, "super_pop"]].rename(columns={sample_col: "sample"})
    df["super_pop"] = df["super_pop"].str.upper()
    df = df[df["super_pop"].isin(VALID_SUPERPOPS)].copy()
    return df

def open_vcf_for_chr(pattern: str, ch: str) -> Optional[VCF]:
    # try both 'N' and 'chrN'
    path1 = pattern.replace("{chr}", ch)
    path2 = pattern.replace("{chr}", with_chr_prefix(ch))
    if os.path.exists(path1):
        return VCF(path1)
    if os.path.exists(path2):
        return VCF(path2)
    # also try the case where pattern already had chr and we pass plain N
    path3 = pattern.replace("{chr}", ch.replace("chr",""))
    if os.path.exists(path3):
        return VCF(path3)
    return None

def ensure_same_sample_order(vcfs: List[VCF]) -> List[str]:
    if not vcfs:
        raise RuntimeError("No VCFs opened.")
    base = vcfs[0].samples
    for v in vcfs[1:]:
        if v.samples != base:
            raise RuntimeError("Sample order differs across chromosome VCFs. Merge first or harmonize order.")
    return base

def parse_trait_anc_from_include(filename: str) -> Tuple[str, str]:
    """
    traits/<trait>__<ANC>.include.tsv  -> (trait, ANC)
    traits/<trait>.include.tsv         -> (trait, 'ALL')
    """
    stem = os.path.basename(filename)
    m = re.match(r"(.+)\.include\.tsv$", stem)
    if not m:
        raise RuntimeError(f"Include filename not recognized: {filename}")
    core = m.group(1)
    if "__" in core:
        trait, anc = core.split("__", 1)
        return trait, anc
    return core, "ALL"

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
                return 23 if x == "X" else int(x)
            for i in range(to_num(a), to_num(b) + 1):
                out.add(with_chr_prefix("X" if i == 23 else str(i)))
        else:
            out.add(with_chr_prefix(tok))
    return sorted(out, key=lambda c: chrom_key(c))

# -------------------------- loading weights & includes -------------------------- #

def load_weights_dir(weights_dir: str) -> pd.DataFrame:
    files = [os.path.join(weights_dir, f) for f in os.listdir(weights_dir) if f.endswith(".tsv") or f.endswith(".tsv.gz")]
    if not files:
        raise RuntimeError(f"No TSV files found in {weights_dir}")
    parts = []
    for f in files:
        df = pd.read_csv(f, sep="\t", dtype=str)
        for col in ["trait","ancestry","chr","pos","ref","alt","effect_allele","beta"]:
            if col not in df.columns:
                raise RuntimeError(f"{f}: missing required column '{col}'")
        df["source_file"] = os.path.basename(f)
        parts.append(df[["trait","ancestry","chr","pos","ref","alt","effect_allele","beta","pgs_id","source_file"] if "pgs_id" in df.columns else ["trait","ancestry","chr","pos","ref","alt","effect_allele","beta","source_file"]])
    W = pd.concat(parts, ignore_index=True)
    W["chr"] = W["chr"].astype(str).apply(with_chr_prefix)
    W["pos"] = pd.to_numeric(W["pos"], errors="coerce").astype("Int64")
    W["ref"] = W["ref"].str.upper()
    W["alt"] = W["alt"].str.upper()
    W["effect_allele"] = W["effect_allele"].str.upper()
    W["beta"] = pd.to_numeric(W["beta"], errors="coerce")
    # keep only valid SNVs with numeric beta
    snv_mask = W["ref"].str.len().eq(1) & W["alt"].str.len().eq(1) & (W["ref"] != W["alt"])
    W = W[snv_mask & W["pos"].notna() & W["beta"].notna()].copy()
    W["key"] = W["chr"] + ":" + W["pos"].astype(int).astype(str) + ":" + W["ref"] + ":" + W["alt"]
    # sanitize ancestry labels
    W["ancestry"] = W["ancestry"].str.upper()
    return W

def load_include(path: str, include_chroms: List[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in ["chr","pos","ref","alt"]:
        if col not in df.columns:
            raise RuntimeError(f"{path}: include file must have columns chr,pos,ref,alt[,effect_allele]")
    df = df[["chr","pos","ref","alt"] + ([c for c in df.columns if c == "effect_allele"])]
    df["chr"] = df["chr"].astype(str).apply(with_chr_prefix)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    df = df[df["chr"].isin(include_chroms)].copy()
    snv_mask = df["ref"].str.len().eq(1) & df["alt"].str.len().eq(1) & (df["ref"] != df["alt"])
    df = df[snv_mask & df["pos"].notna()].copy()
    df["key"] = df["chr"] + ":" + df["pos"].astype(int).astype(str) + ":" + df["ref"] + ":" + df["alt"]
    return df

# -------------------------- PRS scoring on reference -------------------------- #

def compute_reference_scores(
    vcf_pattern: str,
    include_keys: pd.DataFrame,
    weights_tbl: pd.DataFrame,
    label_samples: List[str],
    vcfs_cache: Optional[Dict[str, VCF]] = None
) -> Tuple[np.ndarray, Dict[str, int], Dict[str, int]]:
    """
    Stream across chromosomes, add beta * dosage_effect to each labeled sample.
    dosage_effect = DS_ALT if effect_allele==ALT else (2 - DS_ALT), but we use GT hard calls.
    Returns:
      scores: np.ndarray (n_samples,)
      used_counts: dict with 'expected', 'found', 'allele_mismatch', 'not_found'
      per_chrom_found: dict chr -> found variants
    """
    # Build per-chr map of target variants with their weight and effect allele
    per_chr: Dict[str, Dict[int, Tuple[str,str,float,str]]] = {}
    for _, r in weights_tbl.iterrows():
        ch = r["chr"]; pos = int(r["pos"]); ref = r["ref"]; alt = r["alt"]
        beta = float(r["beta"]); ea = r["effect_allele"]
        per_chr.setdefault(ch, {})[pos] = (ref, alt, beta, ea)

    expected = weights_tbl.shape[0]
    found = 0
    allele_mismatch = 0
    not_found = 0

    # Open one VCF to pull sample order
    v0 = open_vcf_for_chr(vcf_pattern, "1")
    if v0 is None:
        v0 = open_vcf_for_chr(vcf_pattern, "chr1")
    if v0 is None:
        raise RuntimeError("Could not open chr1 VCF with the provided --vcf_pattern.")
    all_samples = v0.samples
    v0.close()

    # Map label samples to indices in VCF order
    sample_index = np.array([all_samples.index(s) for s in label_samples], dtype=np.int64)
    n = len(sample_index)
    scores = np.zeros(n, dtype=np.float64)

    per_chrom_found = {}

    # Iterate chromosomes in order
    for ch in sorted(per_chr.keys(), key=chrom_key):
        keep_positions = per_chr[ch]  # pos -> (ref, alt, beta, ea)
        pos_set = set(keep_positions.keys())
        v = open_vcf_for_chr(vcf_pattern, ch)
        if v is None:
            v = open_vcf_for_chr(vcf_pattern, ch.replace("chr",""))
        if v is None:
            # No file: count all as not found
            not_found += len(pos_set)
            per_chrom_found[ch] = 0
            continue

        found_here = 0
        for var in v:
            pos = var.POS
            if pos not in pos_set:
                continue
            want_ref, want_alt, beta, ea = keep_positions[pos]
            # only biallelic SNVs
            if len(var.ALT) != 1 or len(var.REF) != 1:
                allele_mismatch += 1
                continue
            ref, alt = var.REF, var.ALT[0]
            if not (ref == want_ref and alt == want_alt):
                # strict: skip if allele order differs
                allele_mismatch += 1
                continue

            # Hard-call dosage of ALT
            G = np.asarray(var.genotypes, dtype=np.int16)  # (n_all, 3)
            gt = (G[:,0] + G[:,1]).astype(np.float32)
            # Missing set to NaN then impute to mean
            miss = (G[:,0] < 0) | (G[:,1] < 0)
            gt[miss] = np.nan
            gt = gt[sample_index]  # labeled subset

            # Impute missing to mean (2p) within labeled subset
            if np.any(np.isnan(gt)):
                m = np.nanmean(gt)
                if np.isnan(m):
                    # all missing (unlikely) -> skip
                    continue
                gt = np.where(np.isnan(gt), m, gt)

            # Transform to effect-allele dosage
            if ea == alt:
                eff = gt  # ALT dosage
            elif ea == ref:
                eff = 2.0 - gt  # REF dosage
            else:
                # should not happen for valid SNVs
                allele_mismatch += 1
                continue

            scores += beta * eff
            found += 1
            found_here += 1

        per_chrom_found[ch] = found_here
        # Count remaining positions in this chromosome as not found
        not_found += len(pos_set) - found_here
        v.close()

    used_counts = {
        "expected": int(expected),
        "found": int(found),
        "allele_mismatch": int(allele_mismatch),
        "not_found": int(not_found)
    }
    return scores, used_counts, per_chrom_found

# -------------------------- main -------------------------- #

def main():
    ap = argparse.ArgumentParser(description="Build ECDF calibration per trait × ancestry from 1000G.")
    ap.add_argument("--weights_dir", required=True, help="Directory with harmonized weights (*.tsv)")
    ap.add_argument("--traits_dir", required=True, help="Directory with include lists (*.include.tsv)")
    ap.add_argument("--vcf_pattern", required=True, help="Path pattern with '{chr}' placeholder to 1000G VCFs")
    ap.add_argument("--labels", required=True, help="Labels file with columns: sample, super_pop")
    ap.add_argument("--out", default="calibration", help="Output base directory (default: calibration/)")
    ap.add_argument("--anc", default=None, help="Comma-separated ancestries to calibrate (e.g., EUR,AFR). Default: all present in weights.")
    ap.add_argument("--quantiles", type=int, default=1001, help="Number of ECDF points (default 1001)")
    ap.add_argument("--min-variants", type=int, default=20, help="Skip calibrations with fewer than this many variants used")
    ap.add_argument("--include-chroms", default="1-22", help='Chromosomes to include (default "1-22")')
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    include_chroms = parse_include_chroms(args.include_chroms)

    # Load labels (reference samples and superpops)
    labels_df = load_labels(args.labels)

    # Open a probe VCF just to validate and to lock sample list
    v_probe = open_vcf_for_chr(args.vcf_pattern, "1")
    if v_probe is None:
        v_probe = open_vcf_for_chr(args.vcf_pattern, "chr1")
    if v_probe is None:
        raise RuntimeError("Could not open chr1 VCF with the provided --vcf_pattern")
    v_samples = v_probe.samples
    v_probe.close()

    # Intersect labels with VCF samples, keep sorted by VCF order
    label_mask = labels_df["sample"].isin(set(v_samples))
    labels_df = labels_df[label_mask].copy()
    if labels_df.empty:
        raise RuntimeError("No overlap between label samples and VCF samples.")

    labels_df["vcf_idx"] = labels_df["sample"].map({s:i for i,s in enumerate(v_samples)})
    labels_df.sort_values("vcf_idx", inplace=True)
    label_samples = labels_df["sample"].tolist()

    # Load all harmonized weights into one table
    W = load_weights_dir(args.weights_dir)

    # Ancestry filter (optional)
    if args.anc:
        keep_anc = set(a.strip().upper() for a in args.anc.split(","))
        W = W[W["ancestry"].isin(keep_anc)].copy()
        if W.empty:
            sys.exit("No weights left after ancestry filter.")

    # Find include lists
    include_files = [os.path.join(args.traits_dir, f) for f in os.listdir(args.traits_dir) if f.endswith(".include.tsv")]
    if not include_files:
        raise RuntimeError(f"No include lists (*.include.tsv) found in {args.traits_dir}")

    # Group: trait -> ancestries available in weights
    traits_avail = sorted(W["trait"].unique())
    # For logging & meta
    global_meta = {}

    for inc_path in include_files:
        trait_name, inc_anc = parse_trait_anc_from_include(inc_path)
        if trait_name not in traits_avail:
            print(f"[skip] No weights for trait '{trait_name}' in {args.weights_dir}.")
            continue

        # Determine ancestries to calibrate for this trait
        trait_weights = W[W["trait"] == trait_name].copy()
        ancs_this_trait = sorted(trait_weights["ancestry"].unique())
        if args.anc:
            # already filtered
            pass

        # Load include set (keys only)
        include_df = load_include(inc_path, include_chroms)
        include_keys = set(include_df["key"].tolist())
        n_expected_keys = len(include_keys)
        if n_expected_keys == 0:
            print(f"[skip] Include list {os.path.basename(inc_path)} is empty.")
            continue

        out_trait_dir = os.path.join(args.out, trait_name)
        os.makedirs(out_trait_dir, exist_ok=True)

        trait_meta = {"include_file": os.path.basename(inc_path), "ancestries": {}}

        for anc in ancs_this_trait:
            # If the include file encodes a specific ancestry (not 'ALL'), only build that one
            if inc_anc != "ALL" and anc != inc_anc:
                continue

            Wta = trait_weights[trait_weights["ancestry"] == anc].copy()
            if Wta.empty:
                continue

            # Restrict weights to include keys and requested chromosomes
            Wta = Wta[Wta["chr"].isin(include_chroms)].copy()
            Wta = Wta[Wta["key"].isin(include_keys)].copy()
            if Wta.empty:
                print(f"[warn] No overlapping variants for trait={trait_name}, anc={anc} after intersect with include.")
                continue

            # Deduplicate by key (keep first; log drops)
            dup_mask = Wta.duplicated(subset=["key"], keep="first")
            n_dups = int(dup_mask.sum())
            if n_dups > 0:
                print(f"[info] {n_dups} duplicate weight rows for {trait_name} {anc} → keeping first.")
            Wta_dedup = Wta[~dup_mask].copy()

            # Build per-superpop sample vector for this ancestry
            samples_in_group = labels_df[labels_df["super_pop"] == anc]["sample"].tolist()
            if not samples_in_group:
                print(f"[warn] No labeled samples for ancestry {anc}; skipping.")
                continue

            # Compute PRS on reference samples
            scores, used_counts, per_chrom_found = compute_reference_scores(
                args.vcf_pattern, include_df, Wta_dedup, samples_in_group
            )

            # Sanity: used variant count
            n_used = used_counts["found"]
            if n_used < args.min_variants:
                print(f"[warn] Only used {n_used} variants (<{args.min_variants}) for {trait_name} {anc}; skipping calibration.")
                # record meta anyway
                trait_meta["ancestries"][anc] = {
                    "n_expected_variants": int(n_expected_keys),
                    "n_weights_rows": int(Wta_dedup.shape[0]),
                    "n_variants_used": int(n_used),
                    "skipped": True,
                    "reason": "too_few_variants"
                }
                continue

            # Fit ECDF (empirical quantiles) and Normal params
            q_probs = np.linspace(0.0, 1.0, num=int(args.quantiles), endpoint=True).astype(np.float32)
            q_values = np.quantile(scores, q_probs, method="linear").astype(np.float32)
            mu = float(np.mean(scores))
            sigma = float(np.std(scores, ddof=1)) if scores.size > 1 else 0.0

            # Save calib npz
            calib_path = os.path.join(out_trait_dir, f"calib__{anc}.npz")
            np.savez_compressed(
                calib_path,
                mu=np.array([mu], dtype=np.float32),
                sigma=np.array([sigma], dtype=np.float32),
                q_probs=q_probs,
                q_values=q_values
            )

            # Update meta
            trait_meta["ancestries"][anc] = {
                "n_expected_variants": int(n_expected_keys),
                "n_weights_rows": int(Wta_dedup.shape[0]),
                "n_variants_used": int(n_used),
                "allele_mismatch": int(used_counts["allele_mismatch"]),
                "not_found": int(used_counts["not_found"]),
                "n_samples": int(len(samples_in_group)),
                "quantiles": int(args.quantiles),
                "mu": mu,
                "sigma": sigma,
                "calib_file": os.path.basename(calib_path),
                "per_chrom_found": {k:int(v) for k,v in per_chrom_found.items()}
            }

            print(f"[✓] Calibrated {trait_name} {anc}: n_samples={len(samples_in_group)}, n_vars_used={n_used} → {calib_path}")

        # Write trait-level meta
        with open(os.path.join(out_trait_dir, "meta.json"), "w") as f:
            json.dump(trait_meta, f, indent=2)

        global_meta[trait_name] = trait_meta

    # Also write a top-level meta for convenience
    with open(os.path.join(args.out, "meta.json"), "w") as f:
        json.dump(global_meta, f, indent=2)

    print("\n[✓] Calibration complete.")

if __name__ == "__main__":
    main()
