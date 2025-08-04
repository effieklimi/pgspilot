#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
score_user_and_report.py

Purpose:
  Take a user’s imputed VCF, compute ancestry, compute trait scores, and
  emit percentiles & QC, plus a compact Markdown report per trait.

Inputs:
  --vcf users/<STEM>/<STEM>_imputed_all.vcf.gz
  --pca pca_model/                      # from fit_pca_1kg.py
  --traits_dir traits/                  # include lists (*.include.tsv)
  --weights_dir weights_hm/             # harmonized weights (*.tsv)
  --calib_dir calibration/              # per-trait subdirs with calib__<ANC>.npz
  --out_dir users/<STEM>/               # where to write outputs

Optional:
  --coverage-thresh 0.95
  --pcs <int>                           # default: read from loadings.npy
  --report-md                           # write per-trait report markdown (default on)
  --only-traits "Height,LDL"            # restrict to a subset for debugging

Outputs:
  users/<STEM>/ancestry.json
  users/<STEM>/scores_<trait>.json
  users/<STEM>/report_<trait>.md        (if --report-md)
  users/<STEM>/scores_all.json          (summary across traits)

Core logic:
  1) Project ancestry: DS at pca_sites → standardize with ref_means/stds → PCs
     → kNN classifier → ancestry posteriors w_k.
  2) For each trait × ancestry with weights & calibration:
     - Subset to frozen include list, align effect allele against ALT from VCF:
       use DS if effect==ALT else 2-DS.
     - Compute S_k = Σ β_k,i * DS_i(effect).
     - Coverage QC: n_used / n_expected, plus Σ INFO (R2) if present.
     - Map S_k → p_k via calibration ECDF.
  3) Mixture percentile p_mix = Σ_k w_k * p_k over ancestries available for that trait,
     renormalizing w_k to the set with valid p_k.
  4) QC gate: if coverage_min < threshold OR off_reference, suppress headline percentile.

Notes:
  - Assumes a single sample in the user VCF (your pipeline).
  - Reads DS; if missing, falls back to GP→expected dosage, else GT.
  - INFO R2/DR2/INFO is summed as transparency metric (optional QC).
"""

import argparse
import json
import os
import re
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
from cyvcf2 import VCF
from joblib import load as joblib_load

VALID_A = {"A","C","G","T"}
SUPERPOPS = ["AFR","AMR","EAS","EUR","SAS"]
AUTOSOMES = [f"chr{i}" for i in range(1,23)]

# ------------------------- small helpers ------------------------- #

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

def list_include_files(traits_dir: str) -> List[str]:
    return sorted([os.path.join(traits_dir, f) for f in os.listdir(traits_dir) if f.endswith(".include.tsv")])

def parse_trait_anc_from_include(filename: str) -> Tuple[str, str]:
    # traits/<trait>__<ANC>.include.tsv  -> (trait, ANC)
    # traits/<trait>.include.tsv         -> (trait, 'ALL')
    stem = os.path.basename(filename)
    m = re.match(r"(.+)\.include\.tsv$", stem)
    core = m.group(1)
    if "__" in core:
        trait, anc = core.split("__", 1)
        return trait, anc
    return core, "ALL"

def load_weights_dir(weights_dir: str) -> pd.DataFrame:
    files = [os.path.join(weights_dir, f) for f in os.listdir(weights_dir)
             if f.endswith(".tsv") or f.endswith(".tsv.gz")]
    if not files:
        raise RuntimeError(f"No weights files found in {weights_dir}")
    parts = []
    for f in files:
        df = pd.read_csv(f, sep="\t", dtype=str)
        need = ["trait","ancestry","chr","pos","ref","alt","effect_allele","beta"]
        for col in need:
            if col not in df.columns:
                raise RuntimeError(f"{f}: missing required column '{col}'")
        df["source_file"] = os.path.basename(f)
        parts.append(df[need + ["source_file"] + (["pgs_id"] if "pgs_id" in df.columns else [])])
    W = pd.concat(parts, ignore_index=True)
    # normalize
    W["chr"] = W["chr"].astype(str).apply(with_chr_prefix)
    W["pos"] = pd.to_numeric(W["pos"], errors="coerce").astype("Int64")
    W["ref"] = W["ref"].str.upper()
    W["alt"] = W["alt"].str.upper()
    W["effect_allele"] = W["effect_allele"].str.upper()
    W["beta"] = pd.to_numeric(W["beta"], errors="coerce")
    snv_mask = W["ref"].str.len().eq(1) & W["alt"].str.len().eq(1) & (W["ref"] != W["alt"])
    W = W[snv_mask & W["pos"].notna() & W["beta"].notna()].copy()
    W["ancestry"] = W["ancestry"].str.upper()
    W["key"] = W["chr"] + ":" + W["pos"].astype(int).astype(str) + ":" + W["ref"] + ":" + W["alt"]
    return W

def load_include_file(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = ["chr","pos","ref","alt"]
    for col in need:
        if col not in df.columns:
            raise RuntimeError(f"{path}: include file must have columns chr,pos,ref,alt[,effect_allele]")
    df = df[need + ([c for c in df.columns if c == "effect_allele"])]
    df["chr"] = df["chr"].astype(str).apply(with_chr_prefix)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    snv_mask = df["ref"].str.len().eq(1) & df["alt"].str.len().eq(1) & (df["ref"] != df["alt"])
    df = df[snv_mask & df["pos"].notna()].copy()
    df["key"] = df["chr"] + ":" + df["pos"].astype(int).astype(str) + ":" + df["ref"] + ":" + df["alt"]
    return df

def open_vcf(vcf_path: str) -> VCF:
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(vcf_path)
    return VCF(vcf_path)

def extract_ds_for_sample(var, sample_index: int = 0) -> Optional[float]:
    """
    Prefer DS; if absent, try GP→expected dosage; else GT (0/1/2).
    Returns float or None if not retrievable.
    """
    try:
        ds = var.format("DS")
        if ds is not None:
            # cyvcf2 returns shape (n_samples,) or (n_samples,1)
            val = float(ds[sample_index])
            if np.isfinite(val):
                return val
    except Exception:
        pass
    try:
        gp = var.format("GP")
        if gp is not None:
            # shape (n_samples, 3)
            p = gp[sample_index]
            if p is not None and len(p) >= 3:
                val = float(p[1] + 2.0*p[2])
                if np.isfinite(val):
                    return val
    except Exception:
        pass
    # Fall back to GT (hard call)
    try:
        g = var.genotypes  # (n_samples, 3)
        a1, a2 = g[sample_index][0], g[sample_index][1]
        if a1 >= 0 and a2 >= 0:
            return float(a1 + a2)
    except Exception:
        pass
    return None

def extract_info_r2(var) -> Optional[float]:
    # Beagle: R2; sometimes DR2 or INFO; grab the first numeric we find
    for key in ("R2","DR2","INFO","AR2","ER2"):
        try:
            val = var.INFO.get(key)
            if val is not None:
                f = float(val)
                if np.isfinite(f):
                    return f
        except Exception:
            continue
    return None

def safe_interp_percentile(x: float, q_values: np.ndarray, q_probs: np.ndarray) -> float:
    # Clamp to [0,1] with linear interpolation over monotone q_values
    if x <= q_values[0]:
        return float(q_probs[0])
    if x >= q_values[-1]:
        return float(q_probs[-1])
    return float(np.interp(x, q_values, q_probs))

def write_json(path: str, obj: dict):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)

def format_pct(x: float) -> str:
    return f"{100.0*x:.1f}%" if x is not None and np.isfinite(x) else "—"

# ------------------------- ancestry projection ------------------------- #

def load_pca_model(pca_dir: str):
    sites = pd.read_csv(os.path.join(pca_dir, "pca_sites.b38.tsv"), sep="\t", dtype=str)
    sites["chr"] = sites["chr"].astype(str).apply(with_chr_prefix)
    sites["pos"] = sites["pos"].astype(int)
    sites["ref"] = sites["ref"].str.upper()
    sites["alt"] = sites["alt"].str.upper()
    means = np.load(os.path.join(pca_dir, "ref_means.npy"))  # length M (2p)
    stds = np.load(os.path.join(pca_dir, "ref_stds.npy"))    # length M (sqrt(2p(1-p)))
    loadings = np.load(os.path.join(pca_dir, "loadings.npy")) # shape M x K
    clf_bundle = joblib_load(os.path.join(pca_dir, "classifier.pkl"))  # {"model": ..., "label_encoder": ...}
    # ref scores for off-reference check
    ref_scores = pd.read_csv(os.path.join(pca_dir, "ref_scores.csv"))
    return sites, means, stds, loadings, clf_bundle, ref_scores

def ancestry_from_user(vcf: VCF, sample_index: int, pca_sites: pd.DataFrame,
                       means: np.ndarray, stds: np.ndarray, loadings: np.ndarray,
                       clf_bundle, ref_scores_df: pd.DataFrame) -> Tuple[Dict[str,float], Dict[str,float], Dict[str,float]]:
    """
    Return:
      pcs_dict: {"PC1":...,}
      posteriors: {"AFR":..., "AMR":...}
      diag: {"rad2":..., "rad2_thresh":..., "off_reference": bool}
    """
    M, K = loadings.shape
    # Build lookup: per chr -> pos -> (ref,alt,index)
    per_chr = {}
    for i, r in pca_sites.iterrows():
        per_chr.setdefault(r["chr"], {})[int(r["pos"])] = (r["ref"], r["alt"], i)

    # Pre-allocate DS array
    ds = np.full(M, np.nan, dtype=np.float32)

    # Iterate VCF once across autosomes to fill DS for our positions
    for ch in sorted(per_chr.keys(), key=chrom_key):
        try:
            vcf.set_regions(ch)  # fast region filter in cyvcf2
        except Exception:
            # if indexing / region set fails, fall back to full scan (slower)
            pass
        wanted = per_chr[ch]
        pos_set = set(wanted.keys())
        for var in vcf(ch):
            pos = var.POS
            if pos not in pos_set:
                continue
            ref = var.REF
            if len(var.ALT) != 1:
                continue
            alt = var.ALT[0]
            want_ref, want_alt, idx = wanted[pos]
            if not (ref == want_ref and alt == want_alt):
                continue
            val = extract_ds_for_sample(var, sample_index)
            if val is None:
                continue
            ds[idx] = float(val)

    # Impute missing to means (2p)
    miss = np.isnan(ds)
    if np.any(miss):
        ds[miss] = means[miss]

    # Standardize and project
    z = (ds - means) / np.where(stds > 0, stds, 1.0)
    pcs = z.dot(loadings)  # shape (K,)

    pcs_dict = {f"PC{i+1}": float(pcs[i]) for i in range(pcs.shape[0])}

    # Predict ancestry posteriors
    model = clf_bundle["model"]
    le = clf_bundle["label_encoder"]
    proba = model.predict_proba(pcs.reshape(1, -1))[0]  # aligned with le.classes_
    labels = le.inverse_transform(np.arange(len(proba)))
    post = {lab: float(p) for lab, p in zip(labels, proba)}
    # Ensure we have a value for all 5 superpops (some may be zero/absent)
    for lab in SUPERPOPS:
        post.setdefault(lab, 0.0)

    # Off-reference heuristic: diagonal Mahalanobis radius in PC space
    # using reference PC means/SDs; mark off if radius^2 > ~30 (K≈6).
    ref_pc_cols = [c for c in ref_scores_df.columns if c.startswith("PC")]
    mu_ref = ref_scores_df[ref_pc_cols].mean().values
    sd_ref = ref_scores_df[ref_pc_cols].std(ddof=1).replace(0, 1.0).values
    rad2 = float(np.sum(((pcs - mu_ref) / sd_ref) ** 2))
    rad2_thresh = float(30.0 if len(ref_pc_cols) >= 6 else 20.0)
    off_reference = bool(rad2 > rad2_thresh)

    diag = {"rad2": rad2, "rad2_thresh": rad2_thresh, "off_reference": off_reference}
    return pcs_dict, post, diag

# ------------------------- DS extraction for all traits ------------------------- #

def collect_needed_keys(traits_dir: str, weights_df: pd.DataFrame) -> Tuple[Dict[str, set], set]:
    """
    Return:
      per_chr_positions: dict chr -> set(pos) to prefilter VCF streaming
      include_keys_all: set of keys "chr:pos:ref:alt" needed across all traits
    """
    include_files = list_include_files(traits_dir)
    include_keys_all = set()
    per_chr_positions: Dict[str, set] = {}
    for inc in include_files:
        inc_df = load_include_file(inc)
        # limit keys to those that also exist in weights (any ancestry)
        keys = set(inc_df["key"].tolist()) & set(weights_df["key"].tolist())
        include_keys_all |= keys
        # positions per chr
        for _, r in inc_df.iterrows():
            ch = r["chr"]; pos = int(r["pos"])
            per_chr_positions.setdefault(ch, set()).add(pos)
    # normalize only autosomes by default (traits constructed that way)
    per_chr_positions = {with_chr_prefix(k): v for k, v in per_chr_positions.items()
                         if with_chr_prefix(k) in AUTOSOMES}
    return per_chr_positions, include_keys_all

def scan_user_vcf_for_keys(vcf_path: str,
                           per_chr_positions: Dict[str, set],
                           include_keys_all: set) -> Tuple[Dict[str,float], Dict[str,float]]:
    """
    Single pass through the user's VCF to extract DS (ALT dosage) and R2 (if present)
    for all needed keys. Returns:
      ds_map: key -> DS_ALT
      r2_map: key -> R2/DR2/INFO (if present)
    """
    vcf = open_vcf(vcf_path)
    if len(vcf.samples) != 1:
        # choose the first sample; your pipeline should have a single sample
        sample_index = 0
    else:
        sample_index = 0

    ds_map: Dict[str, float] = {}
    r2_map: Dict[str, float] = {}

    # Iterate chromosome by chromosome to reduce scanning
    for ch in sorted(per_chr_positions.keys(), key=chrom_key):
        try:
            vcf.set_regions(ch)
        except Exception:
            pass
        pos_set = per_chr_positions[ch]
        for var in vcf(ch):
            pos = var.POS
            if pos not in pos_set:
                continue
            if len(var.ALT) != 1 or len(var.REF) != 1:
                continue
            key = f"{with_chr_prefix(var.CHROM)}:{pos}:{var.REF}:{var.ALT[0]}"
            if key not in include_keys_all:
                continue
            val = extract_ds_for_sample(var, sample_index)
            if val is not None and np.isfinite(val):
                ds_map[key] = float(val)
                r2 = extract_info_r2(var)
                if r2 is not None and np.isfinite(r2):
                    r2_map[key] = float(r2)

    vcf.close()
    return ds_map, r2_map

# ------------------------- calibration loading ------------------------- #

def load_calibration_for_trait(calib_dir: str, trait: str) -> Dict[str, dict]:
    """
    Returns dict anc -> {mu, sigma, q_probs, q_values}
    """
    tdir = os.path.join(calib_dir, trait)
    if not os.path.isdir(tdir):
        return {}
    out = {}
    for fn in os.listdir(tdir):
        if not fn.startswith("calib__") or not fn.endswith(".npz"):
            continue
        anc = fn[len("calib__"):-len(".npz")]
        data = np.load(os.path.join(tdir, fn))
        out[anc] = {
            "mu": float(data["mu"][0]),
            "sigma": float(data["sigma"][0]),
            "q_probs": data["q_probs"].astype(np.float32),
            "q_values": data["q_values"].astype(np.float32),
        }
    return out

# ------------------------- trait scoring ------------------------- #

def score_trait_for_ancestry(Wta: pd.DataFrame,
                             include_df: pd.DataFrame,
                             ds_map: Dict[str,float],
                             r2_map: Dict[str,float]) -> Tuple[Optional[float], dict]:
    """
    Compute S_k for this ancestry and trait.
    Returns (S_k or None if no variants used, qc_dict)
    qc_dict contains n_expected, n_used, coverage, sum_info.
    """
    # Restrict to keys in include
    inc_keys = set(include_df["key"].tolist())
    Wta = Wta[Wta["key"].isin(inc_keys)].copy()
    if Wta.empty:
        return None, {"n_expected": int(len(inc_keys)), "n_used": 0, "coverage": 0.0, "sum_info": 0.0}

    # Deduplicate by key (keep first)
    Wta = Wta.drop_duplicates(subset=["key"], keep="first").copy()

    # Build arrays aligned by rows
    keys = Wta["key"].tolist()
    betas = Wta["beta"].astype(float).values
    ealleles = Wta["effect_allele"].values
    alts = Wta["alt"].values
    # fetch DS_ALT for each key
    ds_alt = np.array([ds_map.get(k, np.nan) for k in keys], dtype=np.float32)
    present_mask = np.isfinite(ds_alt)
    n_expected = len(keys)
    n_used = int(present_mask.sum())
    coverage = float(n_used / n_expected) if n_expected > 0 else 0.0
    if n_used == 0:
        return None, {"n_expected": n_expected, "n_used": 0, "coverage": 0.0, "sum_info": 0.0}

    ds_alt_present = ds_alt[present_mask]
    betas_present = betas[present_mask]
    eal_present = ealleles[present_mask]
    alt_present = alts[present_mask]

    # effect dosage: DS if effect==ALT else 2-DS
    eff_ds = np.where(eal_present == alt_present, ds_alt_present, 2.0 - ds_alt_present)
    Sk = float(np.dot(betas_present, eff_ds))

    # sum INFO (optional)
    sum_info = 0.0
    for i, k in enumerate(np.array(keys)[present_mask]):
        r2 = r2_map.get(k)
        if r2 is not None and np.isfinite(r2):
            sum_info += float(r2)

    qc = {"n_expected": n_expected, "n_used": n_used, "coverage": coverage, "sum_info": float(sum_info)}
    return Sk, qc

# ------------------------- report rendering ------------------------- #

def render_markdown(trait: str,
                    mixture_percentile: Optional[float],
                    ancestry_posteriors: Dict[str,float],
                    per_anc: Dict[str, dict],
                    qc: dict) -> str:
    lines = []
    lines.append(f"# {trait}")
    lines.append("")
    if qc.get("low_confidence", False):
        lines.append(f"**Genetic percentile:** _withheld due to low confidence_")
    else:
        if mixture_percentile is not None:
            lines.append(f"**Genetic percentile:** **{100.0*mixture_percentile:.1f}**")
        else:
            lines.append(f"**Genetic percentile:** _not available_")
    lines.append("")
    # ancestry line
    top = sorted(ancestry_posteriors.items(), key=lambda kv: kv[1], reverse=True)
    top_str = " / ".join([f"{anc} {100.0*p:.0f}%" for anc, p in top if p > 0.01][:3])
    if top_str:
        lines.append(f"_Ancestry-aware calibration based on_: {top_str}")
    lines.append("")
    # ancestry-specific details
    lines.append("## Ancestry-specific details")
    for anc, entry in per_anc.items():
        p_k = entry.get("percentile")
        S_k = entry.get("S")
        cov = entry.get("coverage")
        n_exp = entry.get("n_expected")
        n_used = entry.get("n_used")
        sum_info = entry.get("sum_info")
        lines.append(f"- **{anc}**: percentile {format_pct(p_k)}, score S={S_k:.4f} "
                     f"(coverage {cov:.1%}, {n_used}/{n_exp} variants; ΣINFO≈{sum_info:.1f})")
    lines.append("")
    # QC
    lines.append("## Quality checks")
    lines.append(f"- Off-reference ancestry model: **{qc.get('off_reference', False)}**")
    lines.append(f"- Minimum coverage across ancestries: **{qc.get('coverage_min', 0):.1%}** "
                 f"(threshold {qc.get('coverage_thresh',0):.0%})")
    if qc.get("low_confidence", False):
        lines.append(f"- **Low confidence**: {'; '.join(qc.get('reasons', []))}")
    lines.append("")
    lines.append("_This percentile reflects ancestry-aware calibration. Not a medical diagnosis._")
    lines.append("")
    return "\n".join(lines)

# ------------------------- main ------------------------- #

def main():
    ap = argparse.ArgumentParser(description="Score a user, compute ancestry, and report mixture percentiles + QC.")
    ap.add_argument("--vcf", required=True, help="User's imputed VCF (single sample)")
    ap.add_argument("--pca", required=True, help="Directory with pca_model outputs")
    ap.add_argument("--traits_dir", required=True, help="Directory with include lists (*.include.tsv)")
    ap.add_argument("--weights_dir", required=True, help="Directory with harmonized weights (*.tsv)")
    ap.add_argument("--calib_dir", required=True, help="Calibration directory (per-trait subfolders)")
    ap.add_argument("--out_dir", required=True, help="Where to write outputs (users/<STEM>/)")
    ap.add_argument("--coverage-thresh", type=float, default=0.95, help="Coverage threshold (default 0.95)")
    ap.add_argument("--pcs", type=int, default=None, help="Number of PCs (default inferred from loadings)")
    ap.add_argument("--report-md", action="store_true", help="Write per-trait markdown reports")
    ap.add_argument("--only-traits", default=None, help="Comma-separated trait names to restrict scoring")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Load PCA model & ancestry classifier
    pca_sites, ref_means, ref_stds, loadings, clf_bundle, ref_scores_df = load_pca_model(args.pca)
    if args.pcs is not None:
        if args.pcs < loadings.shape[1]:
            loadings = loadings[:, :args.pcs]
            ref_scores_df = ref_scores_df[[c for c in ref_scores_df.columns if c.startswith("PC")][:args.pcs] + ["sample","super_pop"]]
    # Open user VCF once (for ancestry projection we’ll stream again in collect keys)
    vcf = open_vcf(args.vcf)
    sample_index = 0
    sample_id = vcf.samples[sample_index] if vcf.samples else "SAMPLE"

    # Compute ancestry PCs & posteriors
    pcs_dict, posteriors, anc_diag = ancestry_from_user(
        vcf, sample_index, pca_sites, ref_means, ref_stds, loadings, clf_bundle, ref_scores_df
    )
    vcf.close()

    # Normalize ancestry weights to sum to 1
    total_p = sum(posteriors.get(a, 0.0) for a in SUPERPOPS)
    if total_p > 0:
        post_norm = {a: posteriors.get(a, 0.0)/total_p for a in SUPERPOPS}
    else:
        post_norm = {a: 0.0 for a in SUPERPOPS}

    ancestry_json = {
        "sample": sample_id,
        "pcs": pcs_dict,
        "posteriors": post_norm,
        "top2": sorted(post_norm.items(), key=lambda kv: kv[1], reverse=True)[:2],
        "diag": anc_diag
    }
    write_json(os.path.join(args.out_dir, "ancestry.json"), ancestry_json)

    # Load weights and include lists
    W = load_weights_dir(args.weights_dir)

    # Determine which traits to score: by calibration subdirs, intersect with weights & includes
    trait_dirs = [d for d in os.listdir(args.calib_dir) if os.path.isdir(os.path.join(args.calib_dir, d))]
    if args.only_traits:
        keep = set([t.strip() for t in args.only_traits.split(",")])
        trait_dirs = [t for t in trait_dirs if t in keep]
    if not trait_dirs:
        raise RuntimeError("No traits found under --calib_dir (or filtered away).")

    # Collect keys we need from the VCF: union of include keys across traits
    per_chr_positions, include_keys_all = collect_needed_keys(args.traits_dir, W)
    # Extract DS & R2 for all needed keys in one scan
    ds_map, r2_map = scan_user_vcf_for_keys(args.vcf, per_chr_positions, include_keys_all)

    all_trait_summaries = {}

    for trait in sorted(trait_dirs):
        calib_by_anc = load_calibration_for_trait(args.calib_dir, trait)
        if not calib_by_anc:
            print(f"[skip] No calibration files for trait '{trait}'.")
            continue

        # Find includes for this trait
        include_candidates = [f for f in list_include_files(args.traits_dir)
                              if os.path.basename(f).startswith(trait)]
        if not include_candidates:
            print(f"[skip] No include lists for trait '{trait}'.")
            continue

        # Map ancestry -> include file to use (prefer trait__ANC.include.tsv, else trait.include.tsv)
        inc_all = [f for f in include_candidates if f.endswith(f"{trait}.include.tsv")]
        inc_by_anc = {}
        for f in include_candidates:
            if f.endswith(f"{trait}.include.tsv"):
                continue
            # pattern trait__ANC.include.tsv
            stem = os.path.basename(f)[:-len(".include.tsv")]
            if "__" in stem:
                _, anc = stem.split("__", 1)
                inc_by_anc[anc.upper()] = f

        # Which ancestries have weights and calibration?
        trait_W = W[W["trait"] == trait].copy()
        avail_anc_weights = sorted(set(trait_W["ancestry"].unique()))
        avail_anc_calib = sorted(calib_by_anc.keys())
        # Only consider ancestries present in both
        ancestries = [a for a in SUPERPOPS if a in avail_anc_weights and a in avail_anc_calib]
        if not ancestries:
            print(f"[skip] Trait '{trait}': no overlapping ancestries with weights+calibration.")
            continue

        per_anc_results = {}
        coverages = []

        for anc in ancestries:
            # Choose include file
            if anc in inc_by_anc:
                inc_path = inc_by_anc[anc]
            elif inc_all:
                inc_path = inc_all[0]
            else:
                print(f"[warn] Trait '{trait}' {anc}: no include list; skipping ancestry.")
                continue
            inc_df = load_include_file(inc_path)

            # Weights for this ancestry & trait
            Wta = trait_W[trait_W["ancestry"] == anc].copy()
            if Wta.empty:
                continue

            # Compute S_k and QC
            Sk, qc = score_trait_for_ancestry(Wta, inc_df, ds_map, r2_map)
            # Map to percentile
            p_k = None
            if Sk is not None:
                calib = calib_by_anc.get(anc)
                if calib:
                    p_k = safe_interp_percentile(Sk, calib["q_values"], calib["q_probs"])

            coverages.append(qc["coverage"])
            per_anc_results[anc] = {
                "S": float(Sk) if Sk is not None else None,
                "percentile": float(p_k) if p_k is not None else None,
                "coverage": float(qc["coverage"]),
                "n_expected": int(qc["n_expected"]),
                "n_used": int(qc["n_used"]),
                "sum_info": float(qc["sum_info"]),
            }

        # Mixture percentile using available ancestries with percentiles
        available_for_mix = [a for a in ancestries if per_anc_results.get(a, {}).get("percentile") is not None]
        w_raw = {a: post_norm.get(a, 0.0) for a in available_for_mix}
        w_sum = sum(w_raw.values())
        if w_sum > 0:
            w_mix = {a: w_raw[a]/w_sum for a in available_for_mix}
            p_mix = sum(w_mix[a] * per_anc_results[a]["percentile"] for a in available_for_mix)
        else:
            w_mix = {}
            p_mix = None

        # QC gating
        coverage_min = float(min(coverages)) if coverages else 0.0
        reasons = []
        low_conf = False
        if coverage_min < args.coverage_thresh:
            reasons.append(f"coverage {coverage_min:.1%} < {args.coverage_thresh:.0%}")
            low_conf = True
        if ancestry_json["diag"].get("off_reference", False):
            reasons.append("off-reference ancestry")
            low_conf = True

        # Build summary JSON for this trait
        trait_json = {
            "trait": trait,
            "ancestry_posteriors": post_norm,
            "mixture": {
                "available_ancestries": available_for_mix,
                "weights_used": w_mix,
                "p_mix": float(p_mix) if p_mix is not None else None
            },
            "ancestry_specific": per_anc_results,
            "qc": {
                "coverage_min": coverage_min,
                "coverage_thresh": args.coverage_thresh,
                "off_reference": ancestry_json["diag"].get("off_reference", False),
                "low_confidence": low_conf,
                "reasons": reasons
            }
        }

        # Headline: suppress if low confidence
        trait_json["headline_percentile"] = None if low_conf else (float(p_mix) if p_mix is not None else None)

        # Write per-trait JSON
        out_json = os.path.join(args.out_dir, f"scores_{trait}.json")
        write_json(out_json, trait_json)

        # Optional markdown
        if args.report_md:
            md = render_markdown(trait, None if low_conf else p_mix, post_norm, per_anc_results, trait_json["qc"])
            out_md = os.path.join(args.out_dir, f"report_{trait}.md")
            with open(out_md, "w") as f:
                f.write(md)

        all_trait_summaries[trait] = trait_json

    # Write a combined summary
    write_json(os.path.join(args.out_dir, "scores_all.json"), all_trait_summaries)

    print("\n[✓] Scoring complete.")
    print(f"    Ancestry: {os.path.join(args.out_dir, 'ancestry.json')}")
    print(f"    Per-trait JSON: users/<STEM>/scores_<trait>.json")
    if args.report_md:
        print(f"    Reports: users/<STEM>/report_<trait>.md")
    print(f"    Summary: {os.path.join(args.out_dir, 'scores_all.json')}")

if __name__ == "__main__":
    main()

