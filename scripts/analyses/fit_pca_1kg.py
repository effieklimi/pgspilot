#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import argparse
import json
import logging
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from cyvcf2 import VCF
from sklearn.covariance import LedoitWolf
from sklearn.decomposition import PCA
from scipy.stats import chi2

import subprocess
import shutil


# -------------------------------------------------------------
# Constants & helpers
# -------------------------------------------------------------

VALID_A = {"A", "C", "G", "T"}
AUTOSOMES = [str(i) for i in range(1, 23)]

HG38_MHC = ("chr6", 25_000_000, 34_000_000)

# A minimal set of commonly-problematic long-range LD regions on hg38
# (rough but effective; users may add more via --exclude-bed)
LONG_RANGE_REGIONS_HG38 = [
    ("chr8", 7_000_000, 13_000_000),    # 8p23 inversion neighborhood
    ("chr17", 40_000_000, 46_000_000),  # 17q21.31 (MAPT) inversion
]

# I/O utils
def with_chr_prefix(ch: str) -> str:
    ch = str(ch).strip()
    return ch if ch.startswith("chr") else "chr" + ch


def chrom_key(ch: str) -> int:
    ch = with_chr_prefix(ch)
    tok = ch[3:]
    try:
        return int(tok)
    except Exception:
        return 99


def is_biallelic_snp(ref: str, alt: Sequence[str]) -> bool:
    if len(alt) != 1:
        return False
    a = alt[0]
    return (
        len(ref) == 1 and len(a) == 1 and
        ref in VALID_A and a in VALID_A and ref != a
    )

@dataclass
class Region:
    chrom: str
    start: int
    end: int

    def contains(self, chrom: str, pos: int) -> bool:
        return self.chrom == with_chr_prefix(chrom) and self.start <= pos < self.end


def load_bed(path: str) -> List[Region]:
    regs: List[Region] = []
    with open(path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith(('#','track','browser')):
                continue
            toks = line.rstrip().split()[:3]
            if len(toks) < 3:
                continue
            chrom = with_chr_prefix(toks[0])
            try:
                start = int(toks[1])
                end = int(toks[2])
            except ValueError:
                continue
            regs.append(Region(chrom, start, end))
    return regs

def prefilter_with_bcftools(vcf_pattern: str, sites_df: pd.DataFrame, out_dir: str,
                            bcftools: str = "bcftools", tabix: str = "tabix",
                            logger: Optional[logging.Logger] = None) -> str:
    """
    Create per-chromosome VCFs subset to the provided sites using bcftools.
    Returns a new pattern string like f"{out_dir}/subset.chr{chr}.vcf.gz".
    """
    log = logger or logging.getLogger(__name__)
    if shutil.which(bcftools) is None:
        raise RuntimeError("bcftools not found in PATH; install or provide --bcftools path")
    if shutil.which(tabix) is None:
        raise RuntimeError("tabix not found in PATH; install or provide --tabix path")

    os.makedirs(out_dir, exist_ok=True)
    bed_path = os.path.join(out_dir, "panel.bed")

    with open(bed_path, "w") as bed:
        for ch, pos in sites_df[["chr", "pos"]].itertuples(index=False):
            chp = with_chr_prefix(ch)
            p = int(pos)
            bed.write(f"{chp}\t{p-1}\t{p}\n")

    for c in AUTOSOMES:
        in_path = vcf_pattern.replace("{chr}", str(c))
        if not os.path.exists(in_path):
            in_path = vcf_pattern.replace("{chr}", with_chr_prefix(str(c)))
        if not os.path.exists(in_path):
            log.warning("Prefilter: missing input for chr%s; skipping", c)
            continue
        out_path = os.path.join(out_dir, f"subset.chr{c}.vcf.gz")
        subprocess.run([bcftools, "view", "-R", bed_path, in_path, "-Oz", "-o", out_path], check=True)
        subprocess.run([tabix, "-f", out_path], check=True)
        log.info("Prefiltered chr%s → %s", c, out_path)

    return os.path.join(out_dir, "subset.chr{chr}.vcf.gz")


# VCF reading constrained to a panel of sites
def open_vcf_for_chr(pattern: str, ch: str) -> Optional[VCF]:
    path = pattern.replace("{chr}", ch)
    if os.path.exists(path):
        return VCF(path)
    alt = pattern.replace("{chr}", ch.replace("chr", ""))
    if os.path.exists(alt):
        return VCF(alt)
    return None


def ensure_sample_order_same(vcfs: List[VCF]) -> List[str]:
    for v in vcfs[1:]:
        if v.samples != vcfs[0].samples:
            raise RuntimeError("Sample order differs across chromosomes; harmonize first.")
    return vcfs[0].samples


@dataclass
class PanelIndex:
    by_chr_pos_to_j: Dict[str, Dict[int, int]]
    by_chr_pos_alleles: Dict[str, Dict[int, Tuple[str, str]]]
    order: pd.DataFrame  # columns: chr,pos,ref,alt in the exact used order


def build_panel_index(sites: pd.DataFrame) -> PanelIndex:
    s = sites.copy()
    s["chr"] = s["chr"].apply(with_chr_prefix)
    s["pos"] = s["pos"].astype(int)
    s["ref"] = s["ref"].str.upper()
    s["alt"] = s["alt"].str.upper()

    autos = {with_chr_prefix(c) for c in AUTOSOMES}
    s = s[s["chr"].isin(autos)].drop_duplicates(subset=["chr","pos","ref","alt"])\
                               .reset_index(drop=True)
    if s.empty:
        raise RuntimeError("No autosomal PCA sites after filtering input panel.")

    bad = ~s.apply(lambda r: is_biallelic_snp(r["ref"], [r["alt"]]), axis=1)
    if bad.any():
        example = s.loc[bad, ["chr","pos","ref","alt"]].head(5).to_dict(orient="records")
        raise RuntimeError(f"--sites contains non-biallelic/non-ACGT entries, e.g. {example}")

    by_chr_pos_to_j: Dict[str, Dict[int, int]] = {}
    by_chr_pos_alleles: Dict[str, Dict[int, Tuple[str,str]]] = {}
    for j, (ch, pos, ref, alt) in enumerate(s[["chr","pos","ref","alt"]].itertuples(index=False, name=None)):
        by_chr_pos_to_j.setdefault(ch, {})[pos] = j
        by_chr_pos_alleles.setdefault(ch, {})[pos] = (ref, alt)

    return PanelIndex(by_chr_pos_to_j, by_chr_pos_alleles, s)


# Genotype matrix builder
def variant_is_excluded(chrom: str, pos: int, regions: List[Region]) -> bool:
    c = with_chr_prefix(chrom)
    for r in regions:
        if r.chrom == c and r.start <= pos < r.end:
            return True
    return False


def try_match_alleles(ref_v: str, alt_v: str, want_ref: str, want_alt: str,
                      allow_swap: bool, allow_strand_flip: bool) -> Optional[Tuple[int, bool]]:
    """Return (swap, flipped) where swap∈{0,1} indicates REF/ALT swap, flipped indicates strand flip.
       None if not matchable under the allowed rules. Ambiguous A/T or C/G under strand flips are rejected.
    """
    if ref_v == want_ref and alt_v == want_alt:
        return (0, False)

    if not allow_swap and not allow_strand_flip:
        return None

    comp = {"A":"T","T":"A","C":"G","G":"C"}

    candidates: List[Tuple[str,str,int,bool]] = []
    if allow_swap:
        candidates.append((alt_v, ref_v, 1, False))
    if allow_strand_flip:
        if ref_v in comp and alt_v in comp:
            candidates.append((comp[ref_v], comp[alt_v], 0, True))
    if allow_swap and allow_strand_flip:
        if ref_v in comp and alt_v in comp:
            candidates.append((comp[alt_v], comp[ref_v], 1, True))

    ambig_pairs = {("A","T"), ("T","A"), ("C","G"), ("G","C")}

    for r,a,swap,flip in candidates:
        if flip and (r,a) in ambig_pairs:
            continue
        if r == want_ref and a == want_alt:
            return (swap, flip)
    return None


def load_genotype_matrix(vcfs: List[VCF], panel: PanelIndex, label_samples: List[str],
                         exclude_regions: List[Region],
                         impute_info_key: Optional[str] = None, impute_min: Optional[float] = None,
                         allow_swap: bool = False, allow_strand_flip: bool = False,
                         logger: Optional[logging.Logger] = None) -> Tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    """Build genotype dosage matrix X for labeled samples and panel sites in order.
       Returns (X, sites_used, filled_mask). X has NaNs for missing calls; dosage is 0/1/2.
    """
    log = logger or logging.getLogger(__name__)
    all_samples = ensure_sample_order_same(vcfs)
    samp_index = np.array([all_samples.index(s) for s in label_samples], dtype=np.int64)

    M = panel.order.shape[0]
    X = np.full((len(label_samples), M), np.nan, dtype=np.float32)
    filled = np.zeros(M, dtype=bool)

    for v in vcfs:
        for var in v:
            chrom = with_chr_prefix(var.CHROM)
            pos = var.POS
            if chrom not in panel.by_chr_pos_to_j:
                continue
            j = panel.by_chr_pos_to_j[chrom].get(pos)
            if j is None:
                continue
            if exclude_regions and variant_is_excluded(chrom, pos, exclude_regions):
                continue
            if not is_biallelic_snp(var.REF, var.ALT):
                continue
            want_ref, want_alt = panel.by_chr_pos_alleles[chrom][pos]
            match = try_match_alleles(var.REF, var.ALT[0], want_ref, want_alt,
                                      allow_swap=allow_swap, allow_strand_flip=allow_strand_flip)
            if match is None:
                continue
            swap, flipped = match

            # Optional: filter by imputation quality (DR2/R2/INFO)
            if impute_info_key is not None and impute_min is not None:
                val = var.INFO.get(impute_info_key)
                if val is None:
                    # If key absent, be conservative and skip
                    continue
                try:
                    if float(val) < float(impute_min):
                        continue
                except Exception:
                    continue

            g = np.asarray(var.genotypes, dtype=np.int16)  # shape (n_samples, 3)
            if g.ndim != 2 or g.shape[1] < 2:
                continue
            gt = (g[:, 0] + g[:, 1]).astype(np.float32)    # 0/1/2; -1 if missing
            miss = (g[:, 0] < 0) | (g[:, 1] < 0)

            # If REF/ALT are swapped relative to desired, invert dosage
            if swap == 1:
                # Swap 0<->2; 1 stays 1
                gt = 2.0 - gt

            gt[miss] = np.nan
            gt = gt[samp_index]

            X[:, j] = gt
            filled[j] = True

    # Subset sites by filled
    sites_used = panel.order.copy()
    if not filled.all():
        n_drop = int((~filled).sum())
        if n_drop > 0:
            log.warning("Dropping %d panel sites not found / mismatched / filtered", n_drop)
        sites_used = sites_used.loc[filled].reset_index(drop=True)
        X = X[:, filled]
    return X, sites_used, filled


# PCA + population params
def zscore_with_reference(X: np.ndarray, means_2p: np.ndarray, stds: np.ndarray) -> np.ndarray:
    Y = X.copy()
    inds = np.where(np.isnan(Y))
    if inds[0].size:
        Y[inds] = np.take(means_2p, inds[1])
    Y -= means_2p
    Y /= stds
    return Y


def fit_reference_pca(X_std: np.ndarray, k: int, random_state: int = 42) -> Tuple[PCA, np.ndarray, np.ndarray]:
    pca = PCA(n_components=k, svd_solver="randomized", random_state=random_state)
    scores = pca.fit_transform(X_std)             # (n_samples, k)
    loadings = pca.components_.T.astype(np.float32)
    return pca, scores, loadings


def compute_ref_stats(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    means = np.nanmean(X, axis=0).astype(np.float32)               # 2p
    ps = (means / 2.0).astype(np.float32)
    stds = np.sqrt(2.0 * ps * (1.0 - ps)).astype(np.float32)
    keep = stds > 1e-6
    return means, stds, keep


def pop_params_from_scores(ref_scores: pd.DataFrame, pc_cols: List[str]) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], List[str]]:
    pops = sorted(ref_scores["super_pop"].unique())
    centroids: Dict[str, np.ndarray] = {}
    inv_covs: Dict[str, np.ndarray] = {}
    for pop in pops:
        S = ref_scores.loc[ref_scores.super_pop == pop, pc_cols].to_numpy()
        mu = S.mean(axis=0)
        cov = LedoitWolf().fit(S).covariance_
        inv = np.linalg.inv(cov)
        centroids[pop] = mu.astype(np.float64)
        inv_covs[pop] = inv.astype(np.float64)
    return centroids, inv_covs, pops


def mahalanobis2(x: np.ndarray, mu: np.ndarray, inv_cov: np.ndarray) -> float:
    d = x - mu
    return float(d @ inv_cov @ d)


# CLI subcommands
def cmd_fit_ref(args: argparse.Namespace) -> None:
    log = logging.getLogger("fit-ref")

    df_lab = pd.read_csv(args.labels, sep=None, engine="python", dtype=str)
    df_lab.columns = [c.strip().lower() for c in df_lab.columns]
    sample_col = "sample" if "sample" in df_lab.columns else ("iid" if "iid" in df_lab.columns else None)
    if sample_col is None or "super_pop" not in df_lab.columns:
        raise RuntimeError("labels file must have columns: sample (or IID), super_pop")
    df_lab = df_lab[[sample_col, "super_pop"]].rename(columns={sample_col: "sample"})
    df_lab["super_pop"] = df_lab["super_pop"].str.upper()
    allowed = {"AFR","AMR","EAS","EUR","SAS"}
    df_lab = df_lab[df_lab.super_pop.isin(allowed)].copy()

    df_sites = pd.read_csv(args.sites, sep=None, engine="python", dtype=str)
    cols = {c.lower(): c for c in df_sites.columns}
    def need(*names):
        return next((cols[n] for n in names if n in cols), None)
    c_ch = need("chr","chrom","chromosome"); c_pos = need("pos","position")
    c_ref = need("ref","reference","ref_allele"); c_alt = need("alt","alternate","alt_allele")
    if not all([c_ch,c_pos,c_ref,c_alt]):
        raise RuntimeError("--sites must have columns chr,pos,ref,alt (common synonyms ok)")
    df_sites = df_sites[[c_ch,c_pos,c_ref,c_alt]].rename(columns={c_ch:"chr",c_pos:"pos",c_ref:"ref",c_alt:"alt"})

    panel = build_panel_index(df_sites)

    # Determine VCF pattern to use (optionally prefiltered to panel sites)
    vcf_pattern_to_use = args.vcf_pattern
    if args.prefilter_with_bcftools:
        tmp_dir = args.tmp_dir or os.path.join(args.out, "tmp_prefilt")
        vcf_pattern_to_use = prefilter_with_bcftools(
            vcf_pattern=args.vcf_pattern,
            sites_df=panel.order[["chr","pos"]],
            out_dir=tmp_dir,
            bcftools=args.bcftools,
            tabix=args.tabix,
            logger=log,
        )
        log.info("Using prefiltered VCFs: %s", vcf_pattern_to_use)


    exclude_regions: List[Region] = []
    if args.exclude_mhc:
        exclude_regions.append(Region(*HG38_MHC))
    if args.exclude_long_range:
        for ch, s, e in LONG_RANGE_REGIONS_HG38:
            exclude_regions.append(Region(ch, s, e))
    if args.exclude_bed:
        exclude_regions.extend(load_bed(args.exclude_bed))

    # Open per-chromosome VCF/BCF files AFTER deciding pattern
    vcfs: List[VCF] = []
    for ch in AUTOSOMES:
        v = open_vcf_for_chr(vcf_pattern_to_use, ch)
        if v is None:
            v = open_vcf_for_chr(vcf_pattern_to_use, with_chr_prefix(ch))
        if v is None:
            log.warning("Missing VCF for chr%s; continuing", ch)
            continue
        vcfs.append(v)
    if not vcfs:
        raise RuntimeError("No VCF/BCF files could be opened for autosomes 1..22.")

    # Keep only labeled samples that are present in the VCFs
    v0_samples = set(vcfs[0].samples)
    df_lab = df_lab[df_lab["sample"].isin(v0_samples)].copy()
    if df_lab.empty:
        raise RuntimeError("No overlap between labeled samples and VCF samples.")
    df_lab.sort_values("sample", inplace=True)

    X, sites_used, filled_mask = load_genotype_matrix(
        vcfs=vcfs,
        panel=panel,
        label_samples=df_lab["sample"].tolist(),
        exclude_regions=exclude_regions,
        impute_info_key=args.impute_info_key,
        impute_min=args.impute_min,
        allow_swap=args.allow_allele_swap,
        allow_strand_flip=args.allow_strand_flip,
        logger=log,
    )

    if X.shape[1] == 0:
        raise RuntimeError("No panel sites available after filtering/matching.")

    means_2p, stds, keep = compute_ref_stats(X)
    if keep.sum() < X.shape[1]:
        log.warning("Dropping %d monomorphic/near-constant sites before PCA", int(X.shape[1]-keep.sum()))
        X = X[:, keep]
        means_2p = means_2p[keep]
        stds = stds[keep]
        sites_used = sites_used.loc[keep].reset_index(drop=True)

    X_std = zscore_with_reference(X, means_2p, stds)

    k = int(args.pcs)
    pca, scores, loadings = fit_reference_pca(X_std, k=k, random_state=args.random_seed)

    os.makedirs(args.out, exist_ok=True)
    sites_path = os.path.join(args.out, "pca_sites.b38.tsv")
    sites_used.to_csv(sites_path, sep="\t", index=False)

    np.save(os.path.join(args.out, "ref_means.npy"), means_2p.astype(np.float32))
    np.save(os.path.join(args.out, "ref_stds.npy"), stds.astype(np.float32))
    np.save(os.path.join(args.out, "loadings.npy"), loadings.astype(np.float32))
    np.save(os.path.join(args.out, "explained_variance.npy"), pca.explained_variance_.astype(np.float32))
    np.save(os.path.join(args.out, "explained_variance_ratio.npy"), pca.explained_variance_ratio_.astype(np.float32))

    pc_cols = [f"PC{i}" for i in range(1, k+1)]
    ref_scores = pd.DataFrame(scores, columns=pc_cols)
    ref_scores.insert(0, "super_pop", df_lab.super_pop.values)
    ref_scores.insert(0, "sample", df_lab["sample"].values)
    ref_scores.to_csv(os.path.join(args.out, "ref_scores.csv"), index=False)

    centroids, inv_covs, pops = pop_params_from_scores(ref_scores, pc_cols)
    np.savez(os.path.join(args.out, "pop_params.npz"),
             pops=np.array(pops, dtype=object),
             centroids=np.array([centroids[p] for p in pops], dtype=object),
             inv_covs=np.array([inv_covs[p] for p in pops], dtype=object),
             k=np.array(k))

    meta = {
        "n_samples": int(ref_scores.shape[0]),
        "n_sites": int(sites_used.shape[0]),
        "pcs": int(k),
        "labels_file": os.path.abspath(args.labels),
        "sites_file": os.path.abspath(args.sites),
        "vcf_pattern": args.vcf_pattern,
        "exclude_mhc": bool(args.exclude_mhc),
        "exclude_long_range": bool(args.exclude_long_range),
        "exclude_bed": os.path.abspath(args.exclude_bed) if args.exclude_bed else None,
        "impute_info_key": args.impute_info_key,
        "impute_min": args.impute_min,
        "allow_allele_swap": bool(args.allow_allele_swap),
        "allow_strand_flip": bool(args.allow_strand_flip),
    }
    with open(os.path.join(args.out, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    log.info("[✓] Reference PCA written to %s", args.out)


def cmd_project(args: argparse.Namespace) -> None:
    log = logging.getLogger("project")

    means_2p = np.load(os.path.join(args.model, "ref_means.npy"))
    stds = np.load(os.path.join(args.model, "ref_stds.npy"))
    loadings = np.load(os.path.join(args.model, "loadings.npy"))

    df_sites = pd.read_csv(os.path.join(args.model, "pca_sites.b38.tsv"), sep="\t", dtype=str)

    vcfs: List[VCF] = []
    for ch in AUTOSOMES:
        v = open_vcf_for_chr(args.vcf_pattern, ch)
        if v is None:
            v = open_vcf_for_chr(args.vcf_pattern, with_chr_prefix(ch))
        if v is None:
            log.warning("Missing VCF for chr%s; continuing", ch)
            continue
        vcfs.append(v)
    if not vcfs:
        raise RuntimeError("No VCF/BCF files could be opened for autosomes 1..22.")

    samples_all = ensure_sample_order_same(vcfs)
    if args.samples_keep:
        keep = [s.strip() for s in open(args.samples_keep) if s.strip()]
        keep = [s for s in keep if s in samples_all]
        if not keep:
            raise RuntimeError("--samples-keep had no overlap with VCF samples")
        label_samples = keep
    else:
        label_samples = samples_all

    panel = build_panel_index(df_sites)

    # Build genotype matrix for these samples (no extra region exclusion; model already did)
    X, sites_used, filled = load_genotype_matrix(
        vcfs=vcfs,
        panel=panel,
        label_samples=label_samples,
        exclude_regions=[],
        impute_info_key=None,
        impute_min=None,
        allow_swap=False,
        allow_strand_flip=False,
        logger=log,
    )

    if X.shape[1] != loadings.shape[0]:
        raise RuntimeError("Projected sites count mismatch with model loadings; ensure same panel.")

    X_std = zscore_with_reference(X, means_2p, stds)
    scores = X_std @ loadings

    pc_cols = [f"PC{i}" for i in range(1, scores.shape[1] + 1)]
    df_scores = pd.DataFrame(scores, columns=pc_cols)
    df_scores.insert(0, "sample", label_samples)

    os.makedirs(args.out, exist_ok=True)
    out_path = os.path.join(args.out, "query_scores.csv")
    df_scores.to_csv(out_path, index=False)
    log.info("[✓] Wrote projected scores: %s", out_path)


def cmd_assign(args: argparse.Namespace) -> None:
    log = logging.getLogger("assign")

    df_scores = pd.read_csv(args.scores)
    d = np.load(args.pop_params, allow_pickle=True)
    pops = list(d["pops"].tolist())
    centroids = list(d["centroids"].tolist())
    inv_covs = list(d["inv_covs"].tolist())
    K_model = int(d["k"]) if "k" in d else len([c for c in df_scores.columns if c.startswith("PC")])

    K_use = args.k or K_model
    pc_cols = [f"PC{i}" for i in range(1, K_use+1)]
    for c in pc_cols:
        if c not in df_scores.columns:
            raise RuntimeError(f"Scores missing column {c}")

    thr95 = chi2.ppf(args.p95, df=K_use)
    thr99 = chi2.ppf(args.p99, df=K_use)

    pop_to_mu = {p: np.asarray(mu)[:K_use] for p, mu in zip(pops, centroids)}
    pop_to_inv = {p: np.asarray(ic)[:K_use, :K_use] for p, ic in zip(pops, inv_covs)}

    rows = []
    X = df_scores[pc_cols].to_numpy()
    for i in range(X.shape[0]):
        x = X[i]
        d2 = {p: mahalanobis2(x, pop_to_mu[p], pop_to_inv[p]) for p in pops}
        ordered = sorted(d2.items(), key=lambda kv: kv[1])
        (best_pop, best_d2), (second_pop, second_d2) = ordered[0], ordered[1]
        if best_d2 > thr99:
            label = "unassigned"
        elif best_d2 > thr95:
            label = "borderline"
        else:
            label = best_pop
        rows.append({
            "sample": df_scores.loc[i, "sample"],
            "label": label,
            "best_pop": best_pop,
            "best_d2": best_d2,
            "second_pop": second_pop,
            "second_d2": second_d2,
            "margin": float(second_d2 - best_d2),
            "thr95": float(thr95),
            "thr99": float(thr99),
            "k_used": int(K_use),
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out, index=False)
    log.info("[✓] Wrote assignments: %s", args.out)


# Argument parsing
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Reference PCA, projection, and Mahalanobis-based ancestry assignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"], help="Logging level")

    sub = p.add_subparsers(dest="command", required=True)

    q = sub.add_parser("fit-ref", help="Fit reference PCA model")
    q.add_argument("--vcf-pattern", required=True, help="Path with {chr} placeholder, e.g. /data/ALL.chr{chr}.vcf.gz")
    q.add_argument("--labels", required=True, help="TSV/CSV with columns: sample, super_pop")
    q.add_argument("--sites", required=True, help="Panel sites (TSV: chr,pos,ref,alt). Order is preserved.")
    q.add_argument("--out", required=True, help="Output directory for PCA model")
    q.add_argument("--pcs", type=int, default=12, help="Number of PCs to keep")
    q.add_argument("--random-state", type=int, default=42, dest="random_seed")

    q.add_argument("--exclude-mhc", action="store_true", help="Exclude xMHC (hg38 chr6:25–34Mb)")
    q.add_argument("--exclude-long-range", action="store_true", help="Exclude built-in long-range LD regions")
    q.add_argument("--exclude-bed", default=None, help="Additional BED of regions to exclude")

    q.add_argument("--allow-allele-swap", action="store_true", help="Allow REF/ALT swap if necessary")
    q.add_argument("--allow-strand-flip", action="store_true", help="Allow strand flip (disallows ambiguous AT/CG)")

    q.add_argument("--impute-info-key", default=None, help="INFO key with imputation r^2 (e.g., DR2,R2,INFO)")
    q.add_argument("--impute-min", type=float, default=None, help="Minimum imputation r^2 to keep variant")
    q.add_argument("--prefilter-with-bcftools", action="store_true", help="Prefilter each chr VCF to the panel sites using bcftools for speed")
    q.add_argument("--bcftools", default="bcftools", help="Path to bcftools executable")
    q.add_argument("--tabix", default="tabix", help="Path to tabix executable")
    q.add_argument("--tmp-dir", default=None, help="Temporary directory for prefiltered VCFs")

    r = sub.add_parser("project", help="Project samples onto reference PCs")
    r.add_argument("--vcf-pattern", required=True, help="Path with {chr} placeholder for query VCFs")
    r.add_argument("--model", required=True, help="Reference PCA model directory from fit-ref")
    r.add_argument("--out", required=True, help="Output directory for scores")
    r.add_argument("--samples-keep", default=None, help="Optional list of sample IDs to project (one per line)")

    s = sub.add_parser("assign", help="Assign ancestry by Mahalanobis on PCs")
    s.add_argument("--scores", required=True, help="CSV with columns sample,PC1..PCK")
    s.add_argument("--pop-params", required=True, help="pop_params.npz from fit-ref")
    s.add_argument("--out", required=True, help="Output CSV for assignments")
    s.add_argument("--k", type=int, default=None, help="#PCs to use (≤ model K); default = model K")
    s.add_argument("--p95", type=float, default=0.95, help="Lower inlier quantile")
    s.add_argument("--p99", type=float, default=0.99, help="Upper inlier quantile")

    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    argv = argv if argv is not None else sys.argv[1:]
    parser = build_parser()
    if not argv:
        parser.print_help(sys.stderr)
        sys.exit(2)

    log_level = "INFO"
    for i, tok in enumerate(argv):
        if tok == "--log-level" and i+1 < len(argv):
            log_level = argv[i+1]
            break
    logging.basicConfig(level=getattr(logging, log_level), format="%(levelname)s: %(message)s")

    args = parser.parse_args(argv)
    if args.command == "fit-ref":
        cmd_fit_ref(args)
    elif args.command == "project":
        cmd_project(args)
    elif args.command == "assign":
        parser.error("The 'assign' subcommand is disabled in this project. Use scripts/analyses/assign_ancestry.py instead.")
    else:
        parser.error("Unknown command")


if __name__ == "__main__":
    main()
