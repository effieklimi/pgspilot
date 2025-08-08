#!/usr/bin/env python3
"""
Build ancestry-aware standardization stats (mean, sd) for each PGS.

Tailored to your layout:
- 1kG high-coverage VCFs per chromosome:
    /app/genome_data/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
- 1kG panel file with subpop labels:
    /app/scripts/helpers/integrated_call_samples_v3.20130502.ALL.panel
- Harmonized weights registry (fixed):
    /app/weights_hm/pgs_harmonization_registry.tsv
- Outputs:
    /app/weights_hm/pgs_standardization.tsv

What it does
------------
1) Converts chr1..22 1kG VCFs to per-chr PGEN with IDs chr:pos:ref:alt
2) Merges them into ONE cohort PFILE for scoring
3) For each (PGS, subpop) row in the registry:
   - Builds minimalist score file (ID A1 BETA)
   - plink2 --score across all samples
   - Joins with panel labels (super-population) and computes mean/sd for that subpop
   - Counts matched variants from .sscore.vars
4) Writes a tidy TSV:
   pgs_id subpop mean sd n_samples n_variants_used trait

Requirements
------------
- plink2 on PATH
- Python: pandas

Run
---
python build_standardization_from_1kg.py
"""

import os, sys, subprocess, tempfile, shutil, glob
import pandas as pd

# ---------- Fixed paths from your setup ----------
VCF_PATTERN = "/app/genome_data/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
LABELS_PATH = "/app/scripts/helpers/integrated_call_samples_v3.20130502.ALL.panel"
REGISTRY    = "/app/weights_hm/pgs_harmonization_registry.tsv"
OUT_TSV     = "/app/weights_hm/pgs_standardization.tsv"
WORK_DIR    = "/app/pca_model/tmp_std"   # will be created; safe to delete/rebuild

SUBPOPS = ["AFR","AMR","EAS","EUR","SAS"]

# -------------------------------------------------

def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.stderr.write(p.stderr)
        raise SystemExit(f"[ERROR] Command failed: {' '.join(cmd)}")
    return p

def ensure_dir(path):
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def build_cohort_pfile():
    """Create one merged PFILE with chr:pos:ref:alt IDs for chr1..22."""
    ensure_dir(WORK_DIR)
    per_chr = []
    # Step 1: per-chr PGEN
    for c in range(1, 23):
        vcf = VCF_PATTERN.replace("{chr}", str(c))
        if not os.path.exists(vcf):
            raise SystemExit(f"[ERROR] Missing VCF: {vcf}")
        out = os.path.join(WORK_DIR, f"cohort_chr{c}")
        per_chr.append(out)
        # Convert to PGEN with standardized IDs
        run([
            "plink2", "--vcf", vcf,
            "--snps-only", "just-acgt", "--max-alleles", "2",
            "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "200", "truncate",
            "--make-pgen", "--out", out
        ])

    # Step 2: merge into one PFILE
    merged_prefix = os.path.join(WORK_DIR, "cohort_merged")
    merge_list = os.path.join(WORK_DIR, "pmerge_list.txt")
    with open(merge_list, "w") as fh:
        for p in per_chr[1:]:
            fh.write(p + "\n")
    # Start with chr1, then merge the rest
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
    ensure_dir(os.path.dirname(OUT_TSV))
    ensure_dir(WORK_DIR)

    # 0) Read labels (panel). Expect columns: sample  pop  super_pop  sex
    lab = pd.read_csv(LABELS_PATH, sep="\t", dtype=str)
    lab.columns = [c.strip().lower() for c in lab.columns]
    # Be robust to header names
    iid_col = "sample" if "sample" in lab.columns else "iid"
    sp_col  = "super_pop" if "super_pop" in lab.columns else "subpop"
    if iid_col not in lab.columns or sp_col not in lab.columns:
        raise SystemExit("[ERROR] LABELS_PATH must have sample (IID) and super_pop/subpop columns.")
    lab = lab.rename(columns={iid_col: "iid", sp_col: "subpop"})
    lab["subpop"] = lab["subpop"].str.upper()

    # 1) Build/refresh cohort PFILE (merged)
    merged_prefix = build_cohort_pfile()

    # 2) Load registry
    reg = pd.read_csv(REGISTRY, sep="\t", dtype=str)
    reg.columns = [c.strip().lower() for c in reg.columns]
    needed = {"pgs_id","ancestry","out_path","trait"}
    if not needed.issubset(reg.columns):
        raise SystemExit(f"[ERROR] Registry missing columns: {needed - set(reg.columns)}")
    # Keep rows that exist on disk, and ancestry in allowed set
    reg = reg[reg["out_path"].apply(os.path.exists)]
    reg = reg[reg["ancestry"].str.upper().isin(SUBPOPS)]
    if reg.empty:
        raise SystemExit("[ERROR] No valid (PGS, subpop) rows found in registry with existing files.")

    # 3) Work dirs for scoring
    scores_dir = os.path.join(WORK_DIR, "scores")
    ensure_dir(scores_dir)

    rows_out = []

    # Deterministic order
    reg = reg.sort_values(["pgs_id","ancestry"])
    for _, rr in reg.iterrows():
        pgs_id  = rr["pgs_id"]
        subpop  = rr["ancestry"].upper()
        weights = rr["out_path"]
        trait   = rr.get("trait","")

        # Build minimal score file (ID A1 BETA) from harmonized weights
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

        # Score the entire cohort
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
        score_col = detect_score_col(df)

        # Join with labels and filter to this subpop
        j = df.merge(lab, left_on="IID", right_on="iid", how="left")
        sub = j[j["subpop"] == subpop]
        n_samples = sub.shape[0]
        if n_samples == 0:
            # no samples for this subpop in labels; skip
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
        raise SystemExit("[ERROR] No rows produced; check registry/panel alignment.")

    out_df = pd.DataFrame(rows_out, columns=["pgs_id","subpop","mean","sd","n_samples","n_variants_used","trait"])
    out_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[DONE] Wrote {OUT_TSV} with {out_df.shape[0]} rows.")

if __name__ == "__main__":
    main()
