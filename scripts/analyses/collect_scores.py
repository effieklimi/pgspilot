#!/usr/bin/env python3
# collect_scores.py
import argparse, json, os, glob, pandas as pd

def read_sscore(path, user_iid):
    df = pd.read_csv(path, sep=r"\s+")
    row = df.loc[df["IID"] == user_iid]
    if row.empty:
        return None
    # plink2 writes SCORE1_SUM by default; fall back to SCORE1_AVG if needed
    for col in ["SCORE1_SUM", "SCORE1_AVG", "SCORE"]:
        if col in row.columns:
            val = float(row.iloc[0][col])
            break
    else:
        raise ValueError(f"No score column found in {path}")
    # Matched variant count is in .sscore.vars lines (we’ll count in caller)
    return val

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True, help="Directory with *.sscore and *.sscore.vars")
    ap.add_argument("--registry", default="/app/weights_hm/pgs_harmonization_registry.tsv")
    ap.add_argument("--standardization", required=True, help="TSV: pgs_id subpop mean sd")
    ap.add_argument("--user-iid", required=True)
    ap.add_argument("--subpop", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    reg = pd.read_csv(args.registry, sep="\t")
    std = pd.read_csv(args.standardization, sep="\t").rename(columns=str.lower)
    std = std.rename(columns={"pgs":"pgs_id","population":"subpop"})
    std["subpop"] = std["subpop"].str.upper()

    payload = {"user_iid": args.user_iid, "subpop": args.subpop, "scores": []}

    for sscore in glob.glob(os.path.join(args.results_dir, f"*.{args.subpop}.sscore")):
        base = os.path.basename(sscore).replace(".sscore", "")
        # expect name like PGS000300.EUR
        pgs_id = base.split(".")[0]

        raw = read_sscore(sscore, args.user_iid)
        if raw is None:
            continue

        vars_path = sscore + ".vars"
        matched = 0
        if os.path.exists(vars_path):
            with open(vars_path) as fh:
                matched = sum(1 for _ in fh) - 1  # minus header line

        # Registry enrich
        r = reg[(reg["pgs_id"] == pgs_id) & (reg["ancestry"] == args.subpop)]
        trait = r["trait"].iloc[0] if not r.empty else ""
        pgs_url = r["pgs_url"].iloc[0] if not r.empty else ""
        out_path = r["out_path"].iloc[0] if not r.empty else ""
        maf_file = r["maf_file"].iloc[0] if not r.empty else ""

        # Standardization (zscore)
        zscore = None
        mean = sd = None
        srow = std[(std["pgs_id"] == pgs_id) & (std["subpop"] == args.subpop)]
        if not srow.empty:
            mean = float(srow["mean"].iloc[0])
            sd = float(srow["sd"].iloc[0])
            if sd and sd > 0:
                zscore = (raw - mean) / sd

        payload["scores"].append({
            "pgs_id": pgs_id,
            "trait": trait,
            "pgs_url": pgs_url,
            "subpop_weights_used": args.subpop,
            "weights_path": out_path,
            "maf_file_used": maf_file,
            "matched_variant_count": matched,
            "raw_score": raw,
            "zscore": zscore,
            "standardization_mean": mean,
            "standardization_sd": sd
        })

    with open(args.out_json, "w") as fh:
        json.dump(payload, fh, indent=2)

if __name__ == "__main__":
    main()
