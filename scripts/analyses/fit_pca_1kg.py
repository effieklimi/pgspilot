#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fit_pca_1kg.py

Build a global-ancestry PCA model from a reference cohort (e.g., 1000G).

Inputs:
  --vcf-pattern "<path>/ALL.chr{chr}.vcf.gz"   # {chr} replaced by 1..22
  --labels labels.tsv                          # columns: sample, super_pop (AFR/AMR/EAS/EUR/SAS)
  --out pca_model/                             # output directory
  [--pcs 6]                                    # number of PCs to keep
  [--maf 0.05]                                 # minimum minor-allele frequency
  [--max-missing 0.05]                         # max fraction missing per variant
  [--thin-kb 50]                               # distance-based thinning (kb)
  [--max-snps 50000]                           # target max number of SNPs
  [--sites pca_sites.preset.tsv]               # optional: preselected sites (chr pos ref alt)
  [--random-seed 42]

Outputs (in --out):
  pca_sites.b38.tsv     (chr pos ref alt; exactly the order used)
  ref_means.npy         (float32; length M; 2p)
  ref_stds.npy          (float32; length M; sqrt(2p(1-p)))
  loadings.npy          (float32; shape M x K; variant loadings)
  ref_scores.csv        (sample, super_pop, PC1..PCK)
  classifier.pkl        (sklearn kNN model with classes AFR/AMR/EAS/EUR/SAS)
  meta.json             (counts, thresholds, file list)

Notes:
  - Autosomes only (chr1..chr22). VCFs can use "1" or "chr1"; we normalize to "chrN".
  - Genotypes are taken as hard calls (0/1/2 ALT alleles). Missing is imputed to 2p.
  - Standardization per site uses (x - 2p) / sqrt(2p(1-p)).
  - PCA uses randomized SVD (efficient for M >> N).
  - Classifier is k-NN (k=100 or fewer if ref N < 100) with distance weights.
"""

import argparse
import os
import sys
import json
import gzip
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional

from cyvcf2 import VCF
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder
from joblib import dump

# -------------------------- utils -------------------------- #

AUTOSOMES = [str(i) for i in range(1,23)]
VALID_A = {"A","C","G","T"}

def with_chr_prefix(ch: str) -> str:
    ch = str(ch).strip()
    if ch.startswith("chr"): return ch
    return "chr" + ch

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
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    # Accept sample or iid
    sample_col = "sample" if "sample" in df.columns else ("iid" if "iid" in df.columns else None)
    if sample_col is None or "super_pop" not in df.columns:
        raise RuntimeError("labels file must have columns: sample (or IID), super_pop")
    df = df[[sample_col, "super_pop"]].rename(columns={sample_col: "sample"})
    df["super_pop"] = df["super_pop"].str.upper()
    # restrict to 5 superpops
    allowed = {"AFR","AMR","EAS","EUR","SAS"}
    df = df[df["super_pop"].isin(allowed)].copy()
    return df

def open_vcf_for_chr(pattern: str, ch: str) -> Optional[VCF]:
    path = pattern.replace("{chr}", ch)
    if not os.path.exists(path):
        # try without 'chr' prefix
        alt_path = pattern.replace("{chr}", ch.replace("chr",""))
        if os.path.exists(alt_path):
            path = alt_path
        else:
            return None
    v = VCF(path)
    return v

def ensure_sample_order_same(vcfs: List[VCF]) -> List[str]:
    # Use the first VCF as canonical
    base = [with_chr_prefix(s) for s in vcfs[0].samples]  # just to avoid accidental numbers? keep raw.
    for v in vcfs[1:]:
        if v.samples != vcfs[0].samples:
            raise RuntimeError("Sample order differs across chromosomes; please harmonize or use PLINK to merge first.")
    return vcfs[0].samples

def iter_selected_variants(v: VCF,
                           keep_positions: Dict[int, Tuple[str,str]],
                           n_samples: int,
                           sample_index: Optional[np.ndarray]=None):
    """
    Iterate variants of a chromosome and yield selected ones as tuples:
    (pos, ref, alt, gt_float_array[n_samples], missing_frac)
    keep_positions: dict pos -> (ref, alt)
    sample_index: optional array of indices to reorder samples (if label subset)
    """
    for var in v:
        if var.CHROM not in (str(v.seqnames),):  # cyvcf2 iteration already by this VCF file
            pass
        pos = var.POS
        if pos not in keep_positions:
            continue
        ref = var.REF
        alt_list = var.ALT
        if not is_biallelic_snp(ref, alt_list):
            continue
        alt = alt_list[0]
        # exact allele match
        want_ref, want_alt = keep_positions[pos]
        if not (ref == want_ref and alt == want_alt):
            # could be swapped due to REF/ALT in source; for PCA we require exact identity to avoid strand issues
            continue
        g = np.asarray(var.genotypes, dtype=np.int16)  # shape (n, 3) [a1,a2,phased]
        gt = (g[:,0] + g[:,1]).astype(np.float32)      # 0/1/2, -2 when missing
        miss = (g[:,0] < 0) | (g[:,1] < 0)
        gt[miss] = np.nan
        if sample_index is not None:
            gt = gt[sample_index]
        missing_frac = np.mean(np.isnan(gt))
        yield pos, ref, alt, gt, missing_frac

# ---------------------- site selection ---------------------- #

def select_pca_sites(args, vcfs: List[VCF], label_samples: List[str]) -> pd.DataFrame:
    """
    Scan autosomes, filter by MAF/missingness, and thin by distance until max_snps.
    Returns dataframe with columns: chr pos ref alt maf missing
    """
    np.random.seed(args.random_seed)

    # Use sample subset in label order
    all_samples = ensure_sample_order_same(vcfs)
    # The PCA will only use samples that have superpop labels
    sample_index = np.array([all_samples.index(s) for s in label_samples], dtype=np.int64)

    selected_rows = []
    total_kept = 0

    for ch in AUTOSOMES:
        v = open_vcf_for_chr(args.vcf_pattern, ch if args.use_chr_prefix else ch.replace("chr",""))
        if v is None:
            # try with/without chr prefix automatically
            v = open_vcf_for_chr(args.vcf_pattern, "chr"+ch if not ch.startswith("chr") else ch)
        if v is None:
            print(f"[warn] VCF for chr{ch} not found; skipping.")
            continue

        # distance-based thinning
        last_kept_pos = -10**12
        min_dist = int(args.thin_kb * 1000)

        for var in v:
            # only consider biallelic SNPs
            if not is_biallelic_snp(var.REF, var.ALT):
                continue
            pos = var.POS
            if pos - last_kept_pos < min_dist:
                continue
            # fetch genotype vector for labeled samples
            g = np.asarray(var.genotypes, dtype=np.int16)
            gt = (g[:,0] + g[:,1]).astype(np.float32)
            miss = (g[:,0] < 0) | (g[:,1] < 0)
            gt[miss] = np.nan
            gt = gt[sample_index]
            missing_frac = np.mean(np.isnan(gt))
            if missing_frac > args.max_missing:
                continue
            # compute p from alt dosage (ignoring nans)
            if np.all(np.isnan(gt)):
                continue
            p = np.nanmean(gt) / 2.0
            maf = min(p, 1.0 - p)
            if maf < args.maf or maf <= 0.0 or maf >= 0.5:
                continue

            # keep
            selected_rows.append((with_chr_prefix(var.CHROM), pos, var.REF, var.ALT[0], maf, missing_frac))
            last_kept_pos = pos
            total_kept += 1

            if total_kept >= args.max_snps:
                break

        v.close()
        if total_kept >= args.max_snps:
            break

    if total_kept == 0:
        raise RuntimeError("No PCA SNPs selected. Check thresholds or inputs.")

    df = pd.DataFrame(selected_rows, columns=["chr","pos","ref","alt","maf","missing"])
    # Make sure chrom labels are unified to "chrN"
    df["chr"] = df["chr"].apply(with_chr_prefix)
    # Order by chr/pos
    df.sort_values(by=["chr","pos"], key=lambda s: s.map(lambda x: (chrom_key(x), x)), inplace=True)
    return df

# ---------------------- PCA build ---------------------- #

def build_pca(args, vcfs: List[VCF], pca_sites: pd.DataFrame, label_df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Returns:
      means (M,), stds (M,), loadings (M,K), ref_scores_df (samples x PCs)
    """
    # Confirm sample order and subset to labeled samples
    all_samples = ensure_sample_order_same(vcfs)
    label_df = label_df[label_df["sample"].isin(all_samples)].copy()
    if label_df.empty:
        raise RuntimeError("None of the labeled samples are present in the VCF(s).")
    label_samples = label_df["sample"].tolist()
    sample_index = np.array([all_samples.index(s) for s in label_samples], dtype=np.int64)
    n = len(label_samples)

    # Map selected sites per chr → dict[pos] = (ref, alt)
    site_order = []
    by_chr: Dict[str, Dict[int, Tuple[str,str]]] = {}
    for _, r in pca_sites.iterrows():
        ch = with_chr_prefix(str(r["chr"]))
        pos = int(r["pos"])
        ref = str(r["ref"])
        alt = str(r["alt"])
        by_chr.setdefault(ch, {})[pos] = (ref, alt)
        site_order.append((ch, pos))

    M = len(site_order)
    print(f"[*] Building genotype matrix for {n} samples × {M} SNPs…")

    # Allocate genotype matrix (float32)
    X = np.empty((n, M), dtype=np.float32)
    X[:] = np.nan

    # Fill X by iterating chromosomes
    col_idx = 0
    for ch in AUTOSOMES:
        ch_key = with_chr_prefix(ch)
        if ch_key not in by_chr:
            continue
        keep_positions = by_chr[ch_key]  # pos -> (ref,alt)

        # open VCF
        v = open_vcf_for_chr(args.vcf_pattern, ch if args.use_chr_prefix else ch.replace("chr",""))
        if v is None:
            v = open_vcf_for_chr(args.vcf_pattern, ch_key)
        if v is None:
            raise RuntimeError(f"VCF for {ch_key} not found during PCA build.")

        # Iterate only selected positions
        # Build a quick set for membership test
        pos_set = set(keep_positions.keys())

        for var in v:
            pos = var.POS
            if pos not in pos_set:
                continue
            if not is_biallelic_snp(var.REF, var.ALT):
                continue
            ref = var.REF
            alt = var.ALT[0]
            want_ref, want_alt = keep_positions[pos]
            if not (ref == want_ref and alt == want_alt):
                # skip mismatched (rare)
                continue
            g = np.asarray(var.genotypes, dtype=np.int16)
            gt = (g[:,0] + g[:,1]).astype(np.float32)
            miss = (g[:,0] < 0) | (g[:,1] < 0)
            gt[miss] = np.nan
            gt = gt[sample_index]  # reorder to label order
            X[:, col_idx] = gt
            col_idx += 1
            if col_idx % 5000 == 0:
                print(f"    filled {col_idx}/{M} variants…", flush=True)

        v.close()

    if col_idx == 0:
        raise RuntimeError("Failed to fill any PCA variants — allele mismatches or inputs corrupted?")
    if col_idx < M:
        print(f"[warn] Filled {col_idx} of {M} requested SNPs (some mismatched). Truncating to filled set.")
        # Truncate the site list accordingly (keep first col_idx order)
        pca_sites = pca_sites.iloc[:col_idx].copy()
        X = X[:, :col_idx]
        M = col_idx

    # Compute means and stds per site using labeled samples
    means = np.nanmean(X, axis=0)  # mean dosage
    ps = means / 2.0
    stds = np.sqrt(2.0 * ps * (1.0 - ps))
    # Guard against near-zero stds (rare if maf≥0.05)
    keep_mask = stds > 1e-6
    if keep_mask.sum() < M:
        print(f"[warn] Dropping {M - keep_mask.sum()} monomorphic/low-var sites before PCA.")
        X = X[:, keep_mask]
        pca_sites = pca_sites.loc[keep_mask].copy()
        means = means[keep_mask]
        stds = stds[keep_mask]
        ps = ps[keep_mask]
        M = X.shape[1]

    # Impute missing genotypes to mean (2p)
    # (X - means)/std
    # Vectorized imputation
    inds = np.where(np.isnan(X))
    if inds[0].size:
        X[inds] = np.take(means, inds[1])
    X_std = (X - means) / stds

    # Fit PCA (features = variants; samples = rows)
    k = int(args.pcs)
    print(f"[*] Fitting PCA for K={k} PCs on matrix {X_std.shape[0]}×{X_std.shape[1]} …")
    pca = PCA(n_components=k, svd_solver="randomized", random_state=args.random_seed)
    scores = pca.fit_transform(X_std)  # shape (n_samples, k)
    # sklearn PCA components_ shape: (k, n_features). We want loadings (n_features, k)
    loadings = pca.components_.T.astype(np.float32)

    # Prepare reference scores dataframe
    cols = [f"PC{i}" for i in range(1, k+1)]
    ref_scores = pd.DataFrame(scores, columns=cols)
    ref_scores.insert(0, "super_pop", label_df["super_pop"].values)
    ref_scores.insert(0, "sample", label_df["sample"].values)

    # Means/stdevs to save are based on dosage mean (2p) and std sqrt(2p(1-p))
    means_out = (ps * 2.0).astype(np.float32)
    stds_out = stds.astype(np.float32)

    return means_out, stds_out, loadings, ref_scores, pca_sites

# ---------------------- classifier ---------------------- #

def train_classifier(ref_scores: pd.DataFrame, n_neighbors: int = 100):
    X = ref_scores[[c for c in ref_scores.columns if c.startswith("PC")]].values
    y = ref_scores["super_pop"].values
    le = LabelEncoder()
    y_enc = le.fit_transform(y)  # classes sorted alphabetically
    # Keep k <= number of samples
    k = min(n_neighbors, X.shape[0])
    clf = KNeighborsClassifier(n_neighbors=k, weights="distance")
    clf.fit(X, y_enc)
    return clf, le

# ---------------------- main ---------------------- #

def main():
    ap = argparse.ArgumentParser(description="Fit a 1000G-based PCA ancestry model and kNN classifier.")
    ap.add_argument("--vcf-pattern", required=True, help="Path pattern with '{chr}' placeholder, e.g., /data/ALL.chr{chr}.vcf.gz")
    ap.add_argument("--labels", required=True, help="TSV/CSV with columns: sample, super_pop (AFR/AMR/EAS/EUR/SAS)")
    ap.add_argument("--out", required=True, help="Output directory (pca_model/)")
    ap.add_argument("--pcs", type=int, default=6, help="Number of principal components to keep (default 6)")
    ap.add_argument("--maf", type=float, default=0.05, help="MAF threshold for PCA site selection (default 0.05)")
    ap.add_argument("--max-missing", type=float, default=0.05, help="Max missingness per site (default 0.05)")
    ap.add_argument("--thin-kb", type=float, default=50.0, help="Minimum distance in kb between kept SNPs (default 50)")
    ap.add_argument("--max-snps", type=int, default=50000, help="Target maximum SNPs for PCA (default 50,000)")
    ap.add_argument("--sites", default=None, help="Optional preselected sites (TSV with chr pos ref alt), to skip auto-selection")
    ap.add_argument("--random-seed", type=int, default=42, help="Random seed for reproducibility")
    # hidden arg: sometimes inputs are 'chrN' or 'N'. We'll try both automatically; keep flag for clarity
    ap.add_argument("--use-chr-prefix", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Open one VCF to sniff if uses 'chr' prefix
    probe = args.vcf_pattern.replace("{chr}", "1")
    args.use_chr_prefix = os.path.exists(probe) or os.path.exists(args.vcf_pattern.replace("{chr}", "chr1"))

    # Open all available autosome VCFs (some datasets might be missing a chr file; we handle it)
    vcfs = []
    for ch in AUTOSOMES:
        v = open_vcf_for_chr(args.vcf_pattern, ch)
        if v is None:
            v = open_vcf_for_chr(args.vcf_pattern, "chr"+ch)
        if v is None:
            print(f"[warn] Missing VCF for chr{ch}; continuing.")
            continue
        vcfs.append(v)
    if not vcfs:
        raise RuntimeError("No VCF/BCF files could be opened for autosomes 1..22.")

    # Load labels and intersect with VCF samples
    label_df = load_labels(args.labels)
    # Ensure overlap
    v0_samples = set(vcfs[0].samples)
    label_df = label_df[label_df["sample"].isin(v0_samples)].copy()
    if label_df.empty:
        raise RuntimeError("No overlap between label samples and VCF samples.")
    label_df.sort_values("sample", inplace=True)
    print(f"[*] Using {label_df.shape[0]} labeled samples across {len(vcfs)} chromosomes.")

    # Select sites (or use provided)
    if args.sites:
        df_sites = pd.read_csv(args.sites, sep=None, engine="python", dtype=str)
        cols = {c.lower(): c for c in df_sites.columns}
        need = lambda *names: next((cols[n] for n in names if n in cols), None)
        c_ch = need("chr","chrom","chromosome"); c_pos = need("pos","position")
        c_ref = need("ref","reference","ref_allele"); c_alt = need("alt","alternate","alt_allele")
        if not all([c_ch,c_pos,c_ref,c_alt]):
            raise RuntimeError("--sites file must have columns chr,pos,ref,alt")
        df_sites = df_sites[[c_ch,c_pos,c_ref,c_alt]].rename(columns={c_ch:"chr", c_pos:"pos", c_ref:"ref", c_alt:"alt"})
        df_sites["chr"] = df_sites["chr"].apply(with_chr_prefix)
        df_sites["pos"] = df_sites["pos"].astype(int)
        df_sites["ref"] = df_sites["ref"].str.upper()
        df_sites["alt"] = df_sites["alt"].str.upper()
        # autosomes only
        df_sites = df_sites[df_sites["chr"].isin([with_chr_prefix(c) for c in AUTOSOMES])].copy()
        df_sites = df_sites.drop_duplicates(subset=["chr","pos","ref","alt"])
        # we won't enforce maf/missing thresholds in this mode; they'll be handled implicitly by standardization (but monomorphic will drop)
        print(f"[*] Using preselected sites: {df_sites.shape[0]} rows.")
        pca_sites = df_sites[["chr","pos","ref","alt"]].copy()
        pca_sites.sort_values(by=["chr","pos"], key=lambda s: s.map(lambda x: (chrom_key(x), x)), inplace=True)
    else:
        print("[*] Selecting PCA sites by MAF, missingness, and distance thinning…")
        pca_sites = select_pca_sites(args, vcfs, label_df["sample"].tolist())
        print(f"[✓] Selected {pca_sites.shape[0]} PCA SNPs.")

    # Build PCA
    means, stds, loadings, ref_scores, pca_sites_final = build_pca(args, vcfs, pca_sites, label_df)

    # Train classifier
    clf, le = train_classifier(ref_scores, n_neighbors=100)

    # Save outputs
    # 1) PCA sites
    sites_path = os.path.join(args.out, "pca_sites.b38.tsv")
    pca_sites_final[["chr","pos","ref","alt"]].to_csv(sites_path, sep="\t", index=False)

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
        "maf": float(args.maf),
        "max_missing": float(args.max_missing),
        "thin_kb": float(args.thin_kb) if not args.sites else None,
        "max_snps": int(args.max_snps) if not args.sites else None,
        "used_preselected_sites": bool(args.sites),
        "vcf_pattern": args.vcf_pattern,
        "labels_file": os.path.abspath(args.labels),
        "classes": list(sorted(ref_scores["super_pop"].unique())),
    }
    with open(os.path.join(args.out, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"\n[✓] PCA model written to: {args.out}")
    print(f"    Sites:        {sites_path}")
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
