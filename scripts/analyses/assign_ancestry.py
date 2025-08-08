#!/usr/bin/env python3
# assign_ancestry.py
import argparse, pandas as pd, numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-pcs", required=True)   # sample_id, subpop, PC1..PCn (reference)
    ap.add_argument("--user-pcs", required=True)  # sample_id, PC1..PCn (users, projected on same basis)
    ap.add_argument("--out", required=True)
    ap.add_argument("--npc", type=int, default=6)
    ap.add_argument("--reject-quantile", type=float, default=0.99,
                    help="Reject if distance to assigned centroid exceeds this within-pop quantile")
    args = ap.parse_args()

    ref = pd.read_csv(args.ref_pcs, sep="\t")
    usr = pd.read_csv(args.user_pcs, sep="\t")

    pcs = [f"PC{i}" for i in range(1, args.npc+1)]
    for c in pcs:
        if c not in ref.columns or c not in usr.columns:
            raise SystemExit(f"Missing {c} in input")

    # Centroids per subpop
    centroids = ref.groupby("subpop")[pcs].mean()

    # Within-pop distance distribution (for rejection threshold)
    dists = []
    for sp, grp in ref.groupby("subpop"):
        c = centroids.loc[sp].values
        d = np.linalg.norm(grp[pcs].values - c, axis=1)
        thr = np.quantile(d, args.reject_quantile)
        dists.append({"subpop": sp, "thr": thr})
    thr_df = pd.DataFrame(dists).set_index("subpop")

    out_rows = []
    for _, r in usr.iterrows():
        x = r[pcs].values.astype(float)
        # nearest centroid
        diffs = centroids[pcs].values - x
        d = np.linalg.norm(diffs, axis=1)
        idx = d.argmin()
        sp = centroids.index[idx]
        dist = d[idx]
        pct = (dist / thr_df.loc[sp, "thr"]) if thr_df.loc[sp, "thr"] > 0 else np.inf
        label = sp if pct <= 1.0 else "Uncertain"
        out_rows.append({"sample_id": r["sample_id"], "assigned_subpop": label,
                         "distance_over_thr": float(pct)})

    pd.DataFrame(out_rows).to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
