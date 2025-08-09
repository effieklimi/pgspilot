#!/usr/bin/env python3
# /app/scripts/analyses/call_ancestry.py
import argparse, os, sys, tempfile, shutil, subprocess
import numpy as np
import pandas as pd

PCA_DIR = "/app/pca_model"

def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.stderr.write(p.stderr)
        raise SystemExit(f"[ERROR] {' '.join(cmd)}\n{p.stderr}")
    return p

def detect_score_col(df):
    cols = list(map(str, df.columns))
    # Prefer *_SUM, but ignore dosage sum
    sum_cols = [c for c in cols if c.endswith("_SUM") and c != "NAMED_ALLELE_DOSAGE_SUM"]
    if sum_cols:
        return sum_cols[0]
    # Then *_AVG
    avg_cols = [c for c in cols if c.endswith("_AVG")]
    if avg_cols:
        return avg_cols[0]
    # Legacy fallbacks
    for c in ("SCORE1_SUM","SCORE1_AVG","SCORE","SCORE1","BETA_SUM","BETA_AVG"):
        if c in cols:
            return c
    # Last resort: anything starting with SCORE or BETA
    for c in cols:
        if c.startswith(("SCORE","BETA")):
            return c
    raise SystemExit("[ERROR] SCORE column not found in .sscore; got columns: " + ", ".join(cols))



def main():
    ap = argparse.ArgumentParser(description="Project user onto prefit PCA and assign subpopulation.")
    ap.add_argument("--user-pfile", required=True, help="User PLINK2 PFILE prefix (no extension)")
    ap.add_argument("--out-ancestry", required=True, help="TSV: sample_id,assigned_subpop,distance_over_thr")
    ap.add_argument("--out-pcs", required=True, help="TSV: sample_id,PC1..PCk")
    ap.add_argument("--pcs", type=int, default=6)
    ap.add_argument("--reject-quantile", type=float, default=0.99)
    args = ap.parse_args()

    # Load PCA artifacts (trained by your fit_pca_1kg.py)
    sites = pd.read_csv(os.path.join(PCA_DIR, "pca_sites.b38.tsv"), sep="\t", dtype=str).rename(columns=str.lower)
    for c in ("chr","pos","ref","alt"):
        if c not in sites.columns:
            raise SystemExit("[ERROR] pca_sites.b38.tsv missing chr/pos/ref/alt")

    means = np.load(os.path.join(PCA_DIR, "ref_means.npy"))  # = 2p, length M
    stds  = np.load(os.path.join(PCA_DIR, "ref_stds.npy"))   # sqrt(2p(1-p)), length M
    L     = np.load(os.path.join(PCA_DIR, "loadings.npy"))   # shape M x K
    ref   = pd.read_csv(os.path.join(PCA_DIR, "ref_scores.csv"))  # sample, super_pop, PC1..PCK
    ref = ref.rename(columns=str.lower).rename(columns={"super_pop":"subpop"})
    pcs_cols = [c for c in ref.columns if c.startswith("pc")]
    K = min(args.pcs, len(pcs_cols), L.shape[1])
    M = L.shape[0]
    if len(means)!=M or len(stds)!=M:
        raise SystemExit("[ERROR] means/stds length mismatch with loadings")

    # Build scoring betas and intercepts so PLINK2 can do the dot products on raw dosages
    betas = L[:, :K] / stds[:, None]                      # M x K
    intercept = -(betas * means[:, None]).sum(axis=0)     # length K

    # Standardize IDs on the user set (chr:pos:ref:alt)
    work = tempfile.mkdtemp(prefix="ancestry_")
    try:
        std_prefix = os.path.join(work, "user_stdids")
        run(["plink2", "--pfile", args.user_pfile,
             "--set-all-var-ids", "chr@:#:$r:$a", "--new-id-max-allele-len", "200", "truncate",
             "--make-pgen", "--out", std_prefix])

        # Build per-PC score files and score with PLINK2
        ids = (sites["chr"].where(sites["chr"].str.startswith("chr"), "chr"+sites["chr"])
               + ":" + sites["pos"] + ":" + sites["ref"] + ":" + sites["alt"]).tolist()
        a1s = sites["alt"].str.upper().tolist()

        pcs_vals = np.zeros(K, dtype=float)
        iid = None

        # Precompute an index for ids -> row
        id_to_idx = {vid: i for i, vid in enumerate(ids)}
        used_count = None  # filled after first PC

        for k in range(K):
            score_path = os.path.join(work, f"pc{k+1}.score.tsv")
            with open(score_path, "w") as f:
                f.write("ID\tA1\tBETA\n")
                for i in range(M):
                    f.write(f"{ids[i]}\t{a1s[i]}\t{betas[i, k]:.10g}\n")

            out_pref = os.path.join(work, f"pc{k+1}")
            run([
                "plink2", "--pfile", std_prefix,
                "--score", score_path, "1","2","3","header-read",
                "no-mean-imputation",          # <-- add this
                "cols=+scoresums",
                "list-variants",
                "--out", out_pref
            ])
            ssc = pd.read_csv(out_pref + ".sscore", sep=r"\s+")
            if iid is None:
                iid = ssc["IID"].iloc[0]

            if used_count is None:
                used_ids = pd.read_csv(out_pref + ".sscore.vars", header=None)[0].tolist()
                idx = [id_to_idx[v] for v in used_ids if v in id_to_idx]
                used_count = len(idx)
                used_frac = used_count / float(M)

            # Subset-aware intercept on the same idx
            score_sum = float(ssc[detect_score_col(ssc)].iloc[0])
            pcs_vals[k] = score_sum - float((betas[idx, k] * means[idx]).sum())


        # Save user PCs
        with open(args.out_pcs, "w") as g:
            g.write("sample_id" + "".join([f"\tPC{i+1}" for i in range(K)]) + "\n")
            g.write(iid + "".join([f"\t{pcs_vals[i]:.10g}" for i in range(K)]) + "\n")

        # Ancestry assignment (nearest centroid; reject at quantile)
        ref_use = ref[["subpop"] + [f"pc{i+1}" for i in range(K)]].copy()
        centroids = ref_use.groupby("subpop").mean()

        # rejection thresholds per subpop
        thr = {}
        for sp, grp in ref_use.groupby("subpop"):
            c = centroids.loc[sp].values
            d = np.linalg.norm(grp.iloc[:, 1:].values - c, axis=1)
            thr[sp] = np.quantile(d, args.reject_quantile)

        # assign
        x = pcs_vals
        dists = np.linalg.norm(centroids.values - x, axis=1)
        idx = dists.argmin()
        sp  = centroids.index[idx]
        ratio = dists[idx] / (thr[sp] if thr[sp] > 0 else np.inf)
        label = sp if ratio <= 1.0 else "Uncertain"

        pd.DataFrame([{
            "sample_id": iid,
            "assigned_subpop": label,
            "distance_over_thr": float(ratio),
            "used_markers": int(used_count),
            "used_fraction": float(used_frac)
        }]).to_csv(args.out_ancestry, sep="\t", index=False)

        print(f"[OK] Wrote {args.out_pcs} and {args.out_ancestry} (assigned={label})")

    finally:
        shutil.rmtree(work, ignore_errors=True)

if __name__ == "__main__":
    main()
