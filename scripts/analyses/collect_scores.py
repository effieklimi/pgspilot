#!/usr/bin/env python3
# collect_scores.py
import argparse, json, os, glob, subprocess, shutil, sys, re
import pandas as pd

def read_sscore(path, user_iid):
    df = pd.read_csv(path, sep=r"\s+", dtype=str)
    if "IID" not in df.columns:
        raise ValueError(f"IID column not found in {path}")
    row = df.loc[df["IID"].astype(str) == str(user_iid)]
    if row.empty:
        return None
    for col in ["SCORE1_SUM", "SCORE1_AVG", "SCORE"]:
        if col in row.columns:
            val = float(row.iloc[0][col])
            break
    else:
        raise ValueError(f"No score column found in {path}")
    return val

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True, help="Directory containing PGS results; searched recursively for *.sscore and *.sscore.vars")
    ap.add_argument("--registry", default="/app/pgs/weights/harmonized/registry.csv")
    ap.add_argument("--standardization", required=True, help="TSV: pgs_id subpop mean sd")
    ap.add_argument("--user-iid", required=True)
    ap.add_argument("--subpop", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    if not (args.registry and os.path.exists(args.registry) and os.path.getsize(args.registry) > 0):
        raise SystemExit(f"[ERROR] Registry CSV not found or empty: {args.registry}")
    reg = pd.read_csv(args.registry, dtype=str)
    reg.columns = [str(c).strip().lower() for c in reg.columns]
    if "subpopulation" in reg.columns:
        reg["subpopulation"] = reg["subpopulation"].astype(str).str.upper()
    reg["ancestry"] = reg.get("subpopulation", "").astype(str)

    def ensure_standardization_table(std_path: str):
        candidate_paths = [
            std_path,
            "/app/pgs/weights/standardization/standardization.tsv",
        ]

        if os.path.exists(std_path) and os.path.getsize(std_path) > 0:
            return std_path

        for candidate in candidate_paths[1:]:
            if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                # Mirror into requested path for consistency if different
                try:
                    out_dir = os.path.dirname(std_path) or "."
                    os.makedirs(out_dir, exist_ok=True)
                    shutil.copyfile(candidate, std_path)
                    return std_path
                except Exception:
                    return candidate

        build_script = "/app/scripts/analyses/build_standardization.py"
        if os.path.exists(build_script):
            try:
                subprocess.run(["python3", build_script], check=True)
            except Exception as e:
                print(f"[WARN] Failed to build standardization via {build_script}: {e}", file=sys.stderr)
                return None
            for candidate in candidate_paths:
                if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                    if candidate != std_path:
                        try:
                            out_dir = os.path.dirname(std_path) or "."
                            os.makedirs(out_dir, exist_ok=True)
                            shutil.copyfile(candidate, std_path)
                            return std_path
                        except Exception:
                            return candidate
                    return candidate

        print("[WARN] Standardization TSV not found and could not be built. Proceeding without z-scores.", file=sys.stderr)
        return None

    std_path = ensure_standardization_table(args.standardization)
    if std_path and os.path.exists(std_path) and os.path.getsize(std_path) > 0:
        std = pd.read_csv(std_path, sep="\t").rename(columns=str.lower)
    else:
        std = pd.DataFrame(columns=["pgs_id","subpop","mean","sd"])  # empty → zscore stays None
    std = std.rename(columns={"pgs":"pgs_id","population":"subpop"})
    std["subpop"] = std["subpop"].str.upper()

    subpop_norm = str(args.subpop).strip().upper()
    payload = {"user_iid": args.user_iid, "subpop": subpop_norm, "scores": []}

    # Recurse into subdirectories to support per-PGS layout results_dir/PGS_ID/PGS_ID.SUBPOP.*
    sscore_paths = glob.glob(os.path.join(args.results_dir, "**", "*.sscore"), recursive=True)
    for sscore in sorted(sscore_paths):
        base = os.path.basename(sscore).replace(".sscore", "")
        parts = base.split(".")
        file_subpop = parts[-1].upper() if len(parts) >= 2 else ""
        if file_subpop != subpop_norm:
            continue
        prefix = base[: -(len(file_subpop) + 1)] if base.endswith(f".{file_subpop}") else base
        m = re.search(r"(PGS\d+)", prefix, flags=re.IGNORECASE)
        pgs_id = (m.group(1).upper() if m else parts[0].upper())

        raw = read_sscore(sscore, args.user_iid)
        if raw is None:
            continue

        vars_path = sscore + ".vars"
        matched = 0
        if os.path.exists(vars_path):
            with open(vars_path) as fh:
                matched = sum(1 for _ in fh) - 1  # minus header line

        r = reg[(reg["pgs_id"].astype(str) == pgs_id) & (reg["subpopulation"].astype(str) == subpop_norm)]

        trait = (r["trait_reported"].iloc[0] if ("trait_reported" in r.columns and not r.empty) else "")
        out_path = (r["weights_path"].iloc[0] if ("weights_path" in r.columns and not r.empty) else "")
        maf_file = (r["maf_file"].iloc[0] if ("maf_file" in r.columns and not r.empty) else "")
        pgs_url = f"https://www.pgscatalog.org/score/{pgs_id}"

        zscore = None
        mean = sd = None
        srow = std[(std["pgs_id"] == pgs_id) & (std["subpop"] == subpop_norm)]
        if not srow.empty:
            mean = float(srow["mean"].iloc[0])
            sd = float(srow["sd"].iloc[0])
            if sd and sd > 0:
                zscore = (raw - mean) / sd

        payload["scores"].append({
            "pgs_id": pgs_id,
            "trait": trait,
            "pgs_url": pgs_url,
            "subpop_weights_used": subpop_norm,
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
