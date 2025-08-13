#!/usr/bin/env python3
import os, sys, subprocess, tempfile, shutil, glob, argparse
import pandas as pd

# ---------- Defaults (can be overridden by CLI) ----------
VCF_PATTERN = "/app/genome_data/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
LABELS_PATH = "/app/scripts/helpers/integrated_call_samples_v3.20130502.ALL.panel"
REGISTRY    = "/app/pgs/weights/harmonized/registry.tsv"
OUT_TSV     = "/app/pgs/weights/standardization/standardization.tsv"
WORK_DIR    = "/app/pgs/weights/standardization/tmp_std"
SUBPOPS     = ["AFR","AMR","EAS","EUR","SAS"]

def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.stderr.write(p.stderr)
        raise SystemExit(f"[ERROR] Command failed: {' '.join(cmd)}")
    return p

def ensure_dir(path):
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def pfile_exists(prefix: str) -> bool:
    return all(os.path.exists(prefix + ext) for ext in (".pgen", ".pvar", ".psam"))

def build_cohort_pfile(merged_prefix: str, reuse: bool) -> str:
    """Create (or reuse) one merged PFILE with chr:pos:ref:alt IDs for chr1..22."""
    ensure_dir(os.path.dirname(merged_prefix))
    if reuse and pfile_exists(merged_prefix):
        print(f"[HIT] Reusing cohort PFILE at {merged_prefix}", file=sys.stderr)
        return merged_prefix

    work_dir = os.path.dirname(merged_prefix)
    ensure_dir(work_dir)
    per_chr = []
    for c in range(1, 23):
        vcf = VCF_PATTERN.replace("{chr}", str(c))
        if not os.path.exists(vcf):
            raise SystemExit(f"[ERROR] Missing VCF: {vcf}")
        out = os.path.join(work_dir, f"cohort_chr{c}")
        per_chr.append(out)
        if pfile_exists(out):
            continue
        run([
            "plink2", "--vcf", vcf,
            "--snps-only", "just-acgt", "--max-alleles", "2",
            "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "200", "truncate",
            "--make-pgen", "--out", out
        ])

    if not pfile_exists(merged_prefix):
        merge_list = os.path.join(work_dir, "pmerge_list.txt")
        with open(merge_list, "w") as fh:
            for p in per_chr[1:]:
                fh.write(p + "\n")
        run([
            "plink2",
            "--pfile", per_chr[0],
            "--pmerge-list", merge_list,
            "--make-pgen", "--out", merged_prefix
        ])
    return merged_prefix

def detect_score_col(df):
    for c in ("SCORE1_SUM","SCORE1_AVG","SCORE"):
        if c in df.columns: return c
    raise SystemExit("[ERROR] Could not find SCORE column in .sscore")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pgs-id", nargs="+", help="Only build for these PGS IDs (e.g., PGS000300 PGS001234)")
    ap.add_argument("--append", action="store_true", help="Append/dedupe into OUT_TSV instead of overwriting")
    ap.add_argument("--reuse-cohort", action="store_true", help="Reuse existing merged cohort PFILE if present")
    ap.add_argument("--labels", default=LABELS_PATH)
    ap.add_argument("--registry", default=REGISTRY)
    ap.add_argument("--out-tsv", default=OUT_TSV)
    ap.add_argument("--work-dir", default=WORK_DIR)
    ap.add_argument("--cohort-prefix", default=os.path.join(WORK_DIR, "cohort_merged"))
    args = ap.parse_args()

    ensure_dir(os.path.dirname(args.out_tsv))
    ensure_dir(args.work_dir)

    # 0) Read labels
    lab = pd.read_csv(args.labels, sep="\t", dtype=str)
    lab.columns = [c.strip().lower() for c in lab.columns]
    iid_col = "sample" if "sample" in lab.columns else ("iid" if "iid" in lab.columns else None)
    sp_col  = "super_pop" if "super_pop" in lab.columns else ("subpop" if "subpop" in lab.columns else None)
    if not iid_col or not sp_col:
        raise SystemExit("[ERROR] LABELS must include sample/IID and super_pop/subpop columns.")
    lab = lab.rename(columns={iid_col: "iid", sp_col: "subpop"})
    lab["subpop"] = lab["subpop"].str.upper()

    # 1) Build/reuse cohort PFILE
    merged_prefix = build_cohort_pfile(args.cohort_prefix, reuse=args.reuse_cohort)

    # 2) Load registry and optionally restrict to requested PGS IDs
    reg = pd.read_csv(args.registry, sep="\t", dtype=str)
    reg.columns = [c.strip().lower() for c in reg.columns]
    needed = {"pgs_id","ancestry","out_path","trait"}
    if not needed.issubset(reg.columns):
        raise SystemExit(f"[ERROR] Registry missing columns: {needed - set(reg.columns)}")
    reg = reg[reg["out_path"].apply(os.path.exists)]
    reg["ancestry"] = reg["ancestry"].str.upper()
    reg = reg[reg["ancestry"].isin(SUBPOPS)]
    if args.pgs_id:
        want = {x.upper() for x in args.pgs_id}
        reg = reg[reg["pgs_id"].str.upper().isin(want)]
    if reg.empty:
        raise SystemExit("[ERROR] No valid (PGS, subpop) rows to process.")

    # 2b) If appending, skip keys already present
    skip_keys = set()
    if args.append and os.path.exists(args.out_tsv) and os.path.getsize(args.out_tsv) > 0:
        old = pd.read_csv(args.out_tsv, sep="\t", dtype=str)
        old.columns = [c.strip().lower() for c in old.columns]
        if {"pgs_id","subpop"}.issubset(old.columns):
            skip_keys = {(r["pgs_id"], r["subpop"].upper()) for _, r in old.iterrows()}

    rows_out = []
    scores_dir = os.path.join(args.work_dir, "scores")
    ensure_dir(scores_dir)

    reg = reg.sort_values(["pgs_id","ancestry"])
    for _, rr in reg.iterrows():
        pgs_id  = rr["pgs_id"]
        subpop  = rr["ancestry"].upper()
        if (pgs_id, subpop) in skip_keys:
            print(f"[SKIP] {pgs_id} {subpop} already in {args.out_tsv}", file=sys.stderr)
            continue

        weights = rr["out_path"]
        trait   = rr.get("trait","")

        # Build minimal score file (ID A1 BETA)
        score_file = os.path.join(scores_dir, f"{pgs_id}.{subpop}.score.tsv")
        with open(weights) as fin, open(score_file, "w") as fout:
            header = fin.readline().rstrip("\n").split("\t")
            idx = {h:i for i,h in enumerate(header)}
            for req in ("chr","pos","ref","alt","effect_allele","beta"):
                if req not in idx:
                    raise SystemExit(f"[ERROR] Missing '{req}' in {weights}")
            fout.write("ID\tA1\tBETA\n")
            for line in fin:
                t = line.rstrip("\n").split("\t")
                vid = f"{t[idx['chr']]}:{t[idx['pos']]}:{t[idx['ref']]}:{t[idx['alt']]}"
                fout.write(f"{vid}\t{t[idx['effect_allele']]}\t{t[idx['beta']]}\n")

        out_pref = os.path.join(scores_dir, f"{pgs_id}.{subpop}")
        run([
            "plink2",
            "--pfile", merged_prefix,
            "--score", score_file, "1", "2", "3", "header-read", "no-mean-imputation", "list-variants",
            "--out", out_pref
        ])

        sscore = out_pref + ".sscore"
        svars  = out_pref + ".sscore.vars"
        if not os.path.exists(sscore):
            raise SystemExit(f"[ERROR] Missing {sscore}")

        df = pd.read_csv(sscore, sep=r"\s+")
        df.columns = [c.lstrip("#") for c in df.columns]  # robust to '#IID'
        score_col = detect_score_col(df)

        j = df.merge(lab, left_on="IID", right_on="iid", how="left")
        sub = j[j["subpop"] == subpop]
        n_samples = sub.shape[0]
        if n_samples == 0:
            print(f"[WARN] No samples for {subpop}; skipping {pgs_id}.", file=sys.stderr)
            continue

        mean = float(sub[score_col].mean())
        sd   = float(sub[score_col].std(ddof=0))  # population SD

        n_variants = 0
        if os.path.exists(svars):
            with open(svars) as fh:
                n_variants = max(0, sum(1 for _ in fh) - 1)

        rows_out.append({
            "pgs_id": pgs_id,
            "subpop": subpop,
            "mean": mean,
            "sd": sd,
            "n_samples": int(n_samples),
            "n_variants_used": int(n_variants),
            "trait": trait
        })
        print(f"[OK] {pgs_id} {subpop}: mean={mean:.5f} sd={sd:.5f} n={n_samples} vars={n_variants}")

    if not rows_out:
        print("[NOTE] Nothing to write (all rows existed or no data).", file=sys.stderr)
        return

    out_df_new = pd.DataFrame(rows_out, columns=["pgs_id","subpop","mean","sd","n_samples","n_variants_used","trait"])

    # Atomic write (append/dedupe if requested)
    if args.append and os.path.exists(args.out_tsv) and os.path.getsize(args.out_tsv) > 0:
        old = pd.read_csv(args.out_tsv, sep="\t", dtype=str)
        old.columns = [c.strip().lower() for c in old.columns]
        combo = pd.concat([old, out_df_new], ignore_index=True)
        # Keep the last occurrence (new rows win)
        combo = combo.drop_duplicates(subset=["pgs_id","subpop"], keep="last")
        tmp = args.out_tsv + ".tmp"
        combo.to_csv(tmp, sep="\t", index=False)
        os.replace(tmp, args.out_tsv)
        print(f"[DONE] Updated {args.out_tsv} (now {combo.shape[0]} rows).")
    else:
        tmp = args.out_tsv + ".tmp"
        out_df_new.to_csv(tmp, sep="\t", index=False)
        os.replace(tmp, args.out_tsv)
        print(f"[DONE] Wrote {args.out_tsv} with {out_df_new.shape[0]} rows.")

if __name__ == "__main__":
    main()
