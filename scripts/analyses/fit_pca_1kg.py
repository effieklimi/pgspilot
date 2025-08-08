#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fit_pca_1kg.py  (panel-first, simplified)

Build a global-ancestry PCA model from a reference cohort (e.g., 1000G)
USING A PRE-SELECTED, LD-PRUNED SNP PANEL.

Required inputs:
  --vcf-pattern "<path>/ALL.chr{chr}.vcf.gz"   # {chr} replaced by 1..22 or chr1..chr22
  --labels labels.tsv                          # columns: sample, super_pop (AFR/AMR/EAS/EUR/SAS)
  --sites pca_sites.b38.tsv                    # columns: chr, pos, ref, alt (order = order used)
  --out pca_model/                             # output directory

Options:
  [--pcs 6]                                    # number of PCs to keep (default 6)
  [--random-seed 42]

Outputs (in --out):
  pca_sites.b38.tsv     (chr pos ref alt; EXACT order used after any drops)
  ref_means.npy         (float32; length M; 2p)
  ref_stds.npy          (float32; length M; sqrt(2p(1-p)))
  loadings.npy          (float32; shape M x K; variant loadings)
  ref_scores.csv        (sample, super_pop, PC1..PCK)
  classifier.pkl        (sklearn kNN model with classes AFR/AMR/EAS/EUR/SAS)
  meta.json             (counts, files, classes)
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional

from cyvcf2 import VCF
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder
from joblib import dump

# -------------------------- utils -------------------------- #

AUTOSOMES = [str(i) for i in range(1, 23)]
VALID_A = {"A", "C", "G", "T"}

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

def is_biallelic_snp(ref: str, alt: List[str]) -> bool:
    if len(alt) != 1:
        return False
    a = alt[0]
    return (len(ref) == 1 and len(a) == 1 and ref in VALID_A and a in VALID_A and ref != a)

def load_labels(path: str) -> pd.DataFrame:
    print(f"[*] Loading labels from {path}...")
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    sample_col = "sample" if "sample" in df.columns else ("iid" if "iid" in df.columns else None)
    if sample_col is None or "super_pop" not in df.columns:
        raise RuntimeError("labels file must have columns: sample (or IID), super_pop")
    df = df[[sample_col, "super_pop"]].rename(columns={sample_col: "sample"})
    df["super_pop"] = df["super_pop"].str.upper()
    allowed = {"AFR", "AMR", "EAS", "EUR", "SAS"}
    df = df[df["super_pop"].isin(allowed)].copy()
    return df

def open_vcf_for_chr(pattern: str, ch: str) -> Optional[VCF]:
    path = pattern.replace("{chr}", ch)
    if not os.path.exists(path):
        # try flipping chr prefix presence
        alt_path = pattern.replace("{chr}", ch.replace("chr", ""))
        if os.path.exists(alt_path):
            path = alt_path
        else:
            return None
    return VCF(path)

def ensure_sample_order_same(vcfs: List[VCF]) -> List[str]:
    # Use the first VCF as canonical
    for v in vcfs[1:]:
        if v.samples != vcfs[0].samples:
            raise RuntimeError("Sample order differs across chromosomes; harmonize first.")
    return vcfs[0].samples

# ---------------------- PCA build ---------------------- #

def build_pca(vcf_pattern: str,
              vcfs: List[VCF],
              pca_sites: pd.DataFrame,
              label_df: pd.DataFrame,
              pcs: int,
              random_seed: int):
    """
    Returns:
      means_out (M,), stds_out (M,), loadings (M,K), ref_scores_df, pca_sites_final
    Ensures that genotype matrix columns follow the EXACT order of pca_sites rows.
    """
    print(f"[*] Building PCA model...")
    # Confirm sample order and subset to labeled samples
    all_samples = ensure_sample_order_same(vcfs)
    label_df = label_df[label_df["sample"].isin(all_samples)].copy()
    if label_df.empty:
        raise RuntimeError("None of the labeled samples are present in the VCF(s).")
    label_samples = label_df["sample"].tolist()
    sample_index = np.array([all_samples.index(s) for s in label_samples], dtype=np.int64)
    n = len(label_samples)

    # Canonicalize sites but preserve ROW ORDER
    sites = pca_sites.copy()
    sites["chr"] = sites["chr"].apply(with_chr_prefix)
    sites["pos"] = sites["pos"].astype(int)
    sites["ref"] = sites["ref"].str.upper()
    sites["alt"] = sites["alt"].str.upper()

    # Validate that all sites are valid bi-allelic SNPs
    bad_mask = ~sites.apply(lambda r: is_biallelic_snp(r["ref"], [r["alt"]]), axis=1)
    if bad_mask.any():
        bad_rows = sites.loc[bad_mask, ["chr", "pos", "ref", "alt"]].head(5).to_dict(orient="records")
        raise RuntimeError(f"--sites contains non-biallelic or non-ACGT SNPs (e.g., {bad_rows} …)")

    # Restrict to autosomes but KEEP original order among those
    autosomes_set = {with_chr_prefix(c) for c in AUTOSOMES}
    sites = sites[sites["chr"].isin(autosomes_set)].drop_duplicates(subset=["chr", "pos", "ref", "alt"])
    M = sites.shape[0]
    if M == 0:
        raise RuntimeError("No autosomal PCA sites after filtering input panel.")

    print(f"[*] Panel rows (autosomes): {M}. Preserving given order.")

    # Build (per-chr) position → absolute column index maps following the given order
    by_chr_idx: Dict[str, Dict[int, int]] = {}
    for j, (ch, pos) in enumerate(zip(sites["chr"].tolist(), sites["pos"].tolist())):
        by_chr_idx.setdefault(ch, {})[pos] = j

    # Also map expected alleles to verify identity
    by_chr_alleles: Dict[str, Dict[int, Tuple[str, str]]] = {}
    for ch, pos, ref, alt in zip(sites["chr"].tolist(), sites["pos"].tolist(), sites["ref"].tolist(), sites["alt"].tolist()):
        by_chr_alleles.setdefault(ch, {})[pos] = (ref, alt)

    # Allocate genotype matrix (float32), fill with NaN
    X = np.full((n, M), np.nan, dtype=np.float32)
    filled = np.zeros(M, dtype=bool)

    # Iterate chromosomes; open VCFs via pattern (support 1 vs chr1)
    for ch in AUTOSOMES:
        ch_key = with_chr_prefix(ch)
        print(f"[*] Processing chromosome {ch_key} ...", flush=True)
        if ch_key not in by_chr_idx:
            continue

        v = open_vcf_for_chr(vcf_pattern, ch)
        if v is None:
            v = open_vcf_for_chr(vcf_pattern, ch_key)
        if v is None:
            print(f"[warn] Missing VCF for {ch_key}; skipping its sites.")
            continue

        idx_map = by_chr_idx[ch_key]             # pos -> j
        allele_map = by_chr_alleles[ch_key]      # pos -> (ref, alt)

        pos_set = set(idx_map.keys())

        for var in v:
            pos = var.POS
            if pos not in pos_set:
                continue
            if not is_biallelic_snp(var.REF, var.ALT):
                continue
            want_ref, want_alt = allele_map[pos]
            ref = var.REF
            alt = var.ALT[0]
            if not (ref == want_ref and alt == want_alt):
                # strict allele identity to avoid strand issues
                continue

            g = np.asarray(var.genotypes, dtype=np.int16)  # (n_samples, 3)
            gt = (g[:, 0] + g[:, 1]).astype(np.float32)    # 0/1/2
            miss = (g[:, 0] < 0) | (g[:, 1] < 0)
            gt[miss] = np.nan
            gt = gt[sample_index]                          # reorder to label order

            j = idx_map[pos]
            X[:, j] = gt
            filled[j] = True

        v.close()

    # Drop sites never filled (not found / allele-mismatch)
    if not filled.all():
        n_drop = int((~filled).sum())
        if n_drop > 0:
            dropped = sites.loc[~filled, ["chr", "pos"]].head(5).to_dict(orient="records")
            print(f"[warn] Dropping {n_drop} sites not found or allele-mismatched (e.g., {dropped} …)")
        sites = sites.loc[filled].reset_index(drop=True)
        X = X[:, filled]
        M = X.shape[1]

    if M == 0:
        raise RuntimeError("No panel sites could be filled from the VCFs.")

    # Compute means/stds per site among labeled samples
    means = np.nanmean(X, axis=0)          # mean dosage = 2p
    ps = means / 2.0
    stds = np.sqrt(2.0 * ps * (1.0 - ps))

    # Drop monomorphic/near-zero variance sites
    keep_mask = stds > 1e-6
    if keep_mask.sum() < M:
        n_drop = int(M - keep_mask.sum())
        print(f"[warn] Dropping {n_drop} monomorphic/low-variance sites before PCA.")
        X = X[:, keep_mask]
        sites = sites.loc[keep_mask].reset_index(drop=True)
        means = means[keep_mask]
        stds = stds[keep_mask]
        ps = ps[keep_mask]
        M = X.shape[1]

    # Impute missing to mean (2p) and standardize
    inds = np.where(np.isnan(X))
    if inds[0].size:
        X[inds] = np.take(means, inds[1])
    X -= means
    X /= stds
    X_std = X

    # Fit PCA
    k = int(pcs)
    print(f"[*] Fitting PCA for K={k} on matrix {X_std.shape[0]}×{X_std.shape[1]} …")
    pca = PCA(n_components=k, svd_solver="randomized", random_state=random_seed)
    scores = pca.fit_transform(X_std)                 # (n_samples, k)
    loadings = pca.components_.T.astype(np.float32)   # (M, k)

    # Reference scores dataframe
    cols = [f"PC{i}" for i in range(1, k + 1)]
    ref_scores = pd.DataFrame(scores, columns=cols)
    ref_scores.insert(0, "super_pop", label_df["super_pop"].values)
    ref_scores.insert(0, "sample", label_df["sample"].values)

    means_out = (ps * 2.0).astype(np.float32)         # store as 2p
    stds_out = stds.astype(np.float32)

    return means_out, stds_out, loadings, ref_scores, sites

# ---------------------- classifier ---------------------- #

def train_classifier(ref_scores: pd.DataFrame, n_neighbors: int = 100):
    X = ref_scores[[c for c in ref_scores.columns if c.startswith("PC")]].values
    y = ref_scores["super_pop"].values
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    k = min(n_neighbors, X.shape[0])
    clf = KNeighborsClassifier(n_neighbors=k, weights="distance")
    clf.fit(X, y_enc)
    return clf, le

# ---------------------- main ---------------------- #

def main():
    ap = argparse.ArgumentParser(description="Fit a 1000G-based PCA ancestry model using a preselected LD-pruned panel.")
    ap.add_argument("--vcf-pattern", required=True, help="Pattern with '{chr}', e.g., /data/ALL.chr{chr}.vcf.gz")
    ap.add_argument("--labels", required=True, help="TSV/CSV with columns: sample, super_pop (AFR/AMR/EAS/EUR/SAS)")
    ap.add_argument("--sites", required=True, help="Preselected sites (TSV with chr pos ref alt); order is preserved")
    ap.add_argument("--out", required=True, help="Output directory (pca_model/)")
    ap.add_argument("--pcs", type=int, default=6, help="Number of principal components (default 6)")
    ap.add_argument("--random-seed", type=int, default=42, help="Random seed")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    print(f"[*] Panel-first PCA setup")

    # Open autosome VCFs (best-effort — some datasets might lack certain chromosomes)
    vcfs: List[VCF] = []
    for ch in AUTOSOMES:
        v = open_vcf_for_chr(args.vcf_pattern, ch)
        if v is None:
            v = open_vcf_for_chr(args.vcf_pattern, "chr" + ch)
        if v is None:
            print(f"[warn] Missing VCF for chr{ch}; continuing.")
            continue
        vcfs.append(v)
    if not vcfs:
        raise RuntimeError("No VCF/BCF files could be opened for autosomes 1..22.")

    # Load labels and intersect with VCF0 samples
    label_df = load_labels(args.labels)
    v0_samples = set(vcfs[0].samples)
    label_df = label_df[label_df["sample"].isin(v0_samples)].copy()
    if label_df.empty:
        raise RuntimeError("No overlap between label samples and VCF samples.")
    label_df.sort_values("sample", inplace=True)
    print(f"[*] Using {label_df.shape[0]} labeled samples across {len(vcfs)} chromosomes.")

    # Load preselected sites (PRESERVE ORDER)
    df_sites = pd.read_csv(args.sites, sep=None, engine="python", dtype=str)
    cols = {c.lower(): c for c in df_sites.columns}
    def need(*names):
        return next((cols[n] for n in names if n in cols), None)
    c_ch = need("chr", "chrom", "chromosome"); c_pos = need("pos", "position")
    c_ref = need("ref", "reference", "ref_allele"); c_alt = need("alt", "alternate", "alt_allele")
    if not all([c_ch, c_pos, c_ref, c_alt]):
        raise RuntimeError("--sites file must have columns chr,pos,ref,alt (any common synonyms accepted)")
    pca_sites = df_sites[[c_ch, c_pos, c_ref, c_alt]].rename(
        columns={c_ch: "chr", c_pos: "pos", c_ref: "ref", c_alt: "alt"}
    )
    print(f"[*] Loaded preselected sites: {pca_sites.shape[0]} rows (order preserved).")

    # Build PCA
    means, stds, loadings, ref_scores, pca_sites_final = build_pca(
        vcf_pattern=args.vcf_pattern,
        vcfs=vcfs,
        pca_sites=pca_sites,
        label_df=label_df,
        pcs=args.pcs,
        random_seed=args.random_seed
    )

    # Train classifier
    clf, le = train_classifier(ref_scores, n_neighbors=100)

    # Save outputs
    # 1) PCA sites (final order == order used)
    sites_path = os.path.join(args.out, "pca_sites.b38.tsv")
    pca_sites_final[["chr", "pos", "ref", "alt"]].to_csv(sites_path, sep="\t", index=False)

    # 2) Means & stds
    np.save(os.path.join(args.out, "ref_means.npy"), means)
    np.save(os.path.join(args.out, "ref_stds.npy"), stds)

    # 3) Loadings (M x K)
    np.save(os.path.join(args.out, "loadings.npy"), loadings)

    # 4) Reference scores
    ref_scores.to_csv(os.path.join(args.out, "ref_scores.csv"), index=False)

    # 5) Classifier (+ label encoder)
    dump({"model": clf, "label_encoder": le}, os.path.join(args.out, "classifier.pkl"))

    # 6) Meta
    meta = {
        "n_samples": int(ref_scores.shape[0]),
        "n_sites": int(pca_sites_final.shape[0]),
        "pcs": int(args.pcs),
        "vcf_pattern": args.vcf_pattern,
        "labels_file": os.path.abspath(args.labels),
        "sites_file": os.path.abspath(args.sites),
        "classes": list(sorted(ref_scores["super_pop"].unique())),
    }
    with open(os.path.join(args.out, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"\n[✓] PCA model written to: {args.out}")
    print(f"    Sites:        {sites_path}  (order == model order)")
    print(f"    Means/stds:   ref_means.npy, ref_stds.npy")
    print(f"    Loadings:     loadings.npy")
    print(f"    Ref scores:   ref_scores.csv")
    print(f"    Classifier:   classifier.pkl")
    print(f"    Meta:         meta.json")

    # Close VCFs
    for v in vcfs:
        try:
            v.close()
        except Exception:
            pass

if __name__ == "__main__":
    main()
