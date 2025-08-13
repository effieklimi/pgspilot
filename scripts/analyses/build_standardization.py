#!/usr/bin/env python3
import os, sys, subprocess, tempfile, shutil, glob, argparse, shlex, gzip
import pandas as pd

# Defaults (can be overridden by CLI)
VCF_PATTERN = "/app/genome_data/1000G/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
LABELS_PATH = "/app/scripts/helpers/integrated_call_samples_v3.20130502.ALL.panel"
REGISTRY    = "/app/pgs/weights/harmonized/registry.csv"
OUT_TSV     = "/app/pgs/weights/standardization/standardization.tsv"
WORK_DIR    = "/app/pgs/weights/standardization/tmp_std"
SUBPOPS     = ["AFR","AMR","EAS","EUR","SAS"]

def run(cmd):
    print("+", " ".join(shlex.quote(c) for c in cmd), file=sys.stderr, flush=True)
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"[ERROR] Command failed with exit code {e.returncode}: {' '.join(cmd)}")

def ensure_dir(path):
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def pfile_exists(prefix: str) -> bool:
    return all(os.path.exists(prefix + ext) for ext in (".pgen", ".pvar", ".psam"))

def count_lines_excluding_header(path: str) -> int:
    """Fast line count via `wc -l`, minus 1 for header. Returns 0 if file missing."""
    if not os.path.exists(path):
        return 0
    out = subprocess.check_output(["wc", "-l", path], text=True)
    n = int(out.split()[0])
    return max(0, n - 1)

def build_cohort_pfile(merged_prefix: str, reuse: bool, threads: int) -> str:
    """Create (or reuse) one merged PFILE with chr:pos:ref:alt IDs for chr1..22.

    Note: Merging large per-chromosome PGENs can occasionally fail with
    'Invalid variant count in .pgen file' if any partial/old artifact exists.
    We keep this path for callers that explicitly want a merged cohort, but
    the main flow below avoids merging and instead scores per chromosome and
    accumulates results.
    """
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
        def vcf_has_chr_prefix(path: str) -> bool:
            try:
                opener = gzip.open if path.endswith((".gz", ".bgz")) else open
                with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        if line.startswith('#'):
                            continue
                        chrom = line.split('\t', 1)[0]
                        return chrom.startswith('chr')
            except Exception:
                return False
            return False
        has_chr = vcf_has_chr_prefix(vcf)
        set_fmt = "@:#:$r:$a" if has_chr else "chr@:#:$r:$a"
        run([
            "plink2", "--vcf", vcf,
            "--snps-only", "just-acgt", "--max-alleles", "2",
            "--set-all-var-ids", set_fmt, "--new-id-max-allele-len", "200", "truncate",
            "--threads", str(threads),
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
            "--threads", str(threads),
            "--make-pgen", "--out", merged_prefix
        ])
    return merged_prefix

def pvar_has_chr_prefix(prefix: str) -> bool:
    try:
        with open(prefix + ".pvar", "r") as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                chrom = line.split('\t', 1)[0]
                return chrom.startswith('chr')
    except Exception:
        return False
    return False

def build_per_chr_pfiles(work_dir: str, threads: int, expect_chr_prefix: bool = True) -> list:
    """Create per-chromosome PFILEs with standardized variant IDs and return prefixes.
    If existing PFILEs are present but CHROM prefix style mismatches expectation, rebuild them.
    """
    ensure_dir(work_dir)
    prefixes = []
    for c in range(1, 23):
        vcf = VCF_PATTERN.replace("{chr}", str(c))
        if not os.path.exists(vcf):
            raise SystemExit(f"[ERROR] Missing VCF: {vcf}")
        pref = os.path.join(work_dir, f"cohort_chr{c}")
        prefixes.append(pref)
        if pfile_exists(pref):
            # Validate style: rebuild if mismatch with expected chr prefix
            has_chr_now = pvar_has_chr_prefix(pref)
            if has_chr_now != expect_chr_prefix:
                for ext in (".pgen", ".pvar", ".psam"):
                    try:
                        os.remove(pref + ext)
                    except Exception:
                        pass
            else:
                continue
        def vcf_has_chr_prefix(path: str) -> bool:
            try:
                opener = gzip.open if path.endswith((".gz", ".bgz")) else open
                with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        if line.startswith('#'):
                            continue
                        chrom = line.split('\t', 1)[0]
                        return chrom.startswith('chr')
            except Exception:
                return False
            return False
        has_chr = vcf_has_chr_prefix(vcf)
        set_fmt = "@:#:$r:$a" if has_chr else "chr@:#:$r:$a"
        run([
            "plink2", "--vcf", vcf,
            "--snps-only", "just-acgt", "--max-alleles", "2",
            "--set-all-var-ids", set_fmt, "--new-id-max-allele-len", "200", "truncate",
            "--threads", str(threads),
            "--make-pgen", "--out", pref
        ])
    return prefixes

def detect_score_col(df):
    upper_cols = [str(c).strip() for c in df.columns]
    # Prefer *_SUM over *_AVG
    sum_candidates = [c for c in upper_cols if c.endswith('_SUM') and 'ALLELE' not in c.upper() and c.upper() != 'NAMED_ALLELE_DOSAGE_SUM']
    if sum_candidates:
        # If SCORE1_SUM exists, prefer it; else take the first
        if 'SCORE1_SUM' in sum_candidates:
            return 'SCORE1_SUM'
        return sum_candidates[0]
    avg_candidates = [c for c in upper_cols if c.endswith('_AVG') and 'ALLELE' not in c.upper()]
    if avg_candidates:
        if 'SCORE1_AVG' in avg_candidates:
            return 'SCORE1_AVG'
        return avg_candidates[0]
    for c in ("SCORE1_SUM","SCORE1_AVG","SCORE"):
        if c in upper_cols:
            return c
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
    ap.add_argument("--threads", type=int, default=int(os.environ.get("PLINK_THREADS", os.cpu_count() or 1)),
                    help="Number of threads to pass to PLINK2 (default: $PLINK_THREADS or CPU count).")
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
    lab["subpop"] = lab["subpop"].astype(str).str.strip().str.upper()
    lab = lab.drop_duplicates(subset=["iid"])

    # 1) Build per-chromosome PFILEs (avoid a global merge step)
    per_chr_prefixes = build_per_chr_pfiles(args.work_dir, threads=args.threads)

    # 2) Load registry (CSV schema you described)
    reg = pd.read_csv(args.registry, dtype=str)
    reg.columns = [c.strip().lower() for c in reg.columns]
    need = {"pgs_id","subpopulation","weights_path"}
    if not need.issubset(reg.columns):
        missing = need - set(reg.columns)
        raise SystemExit(f"[ERROR] Registry missing required columns: {missing}")
    reg["subpopulation"] = reg["subpopulation"].astype(str).str.upper()
    reg = reg[reg["subpopulation"].isin(SUBPOPS)]
    reg = reg[reg["weights_path"].apply(lambda p: isinstance(p, str) and os.path.exists(p))]
    if args.pgs_id:
        want = {x.upper() for x in args.pgs_id}
        reg = reg[reg["pgs_id"].str.upper().isin(want)]
    if reg.empty:
        raise SystemExit("[ERROR] No valid (PGS, subpop) rows to process.")

    # 2b) Skip keys already present if appending
    skip_keys = set()
    if args.append and os.path.exists(args.out_tsv) and os.path.getsize(args.out_tsv) > 0:
        old = pd.read_csv(args.out_tsv, sep="\t", dtype=str)
        old.columns = [c.strip().lower() for c in old.columns]
        if {"pgs_id","subpop"}.issubset(old.columns):
            skip_keys = {(r["pgs_id"], r["subpop"].upper()) for _, r in old.iterrows()}

    rows_out = []
    scores_dir = os.path.join(args.work_dir, "scores")
    ensure_dir(scores_dir)

    reg = reg.sort_values(["pgs_id","subpopulation"])
    for _, rr in reg.iterrows():
        pgs_id  = rr["pgs_id"]
        subpop  = rr["subpopulation"].upper()
        if (pgs_id, subpop) in skip_keys:
            print(f"[SKIP] {pgs_id} {subpop} already in {args.out_tsv}", file=sys.stderr)
            continue

        weights = rr["weights_path"]
        trait   = rr.get("trait_reported", "")

        acc = None
        total_vars = 0
        for pref in per_chr_prefixes:
            basepref = os.path.basename(pref)
            chr_num = ''.join(ch for ch in basepref if ch.isdigit())
            if not chr_num:
                try:
                    with open(pref + ".pvar", "r") as pv:
                        for _line in pv:
                            if _line.startswith('#'): continue
                            chr_num = _line.split('\t',1)[0].lstrip('chr')
                            break
                except Exception:
                    chr_num = ""
            has_chr_prefix = False
            try:
                with open(pref + ".pvar", "r") as pv:
                    for _line in pv:
                        if _line.startswith('#'): continue
                        has_chr_prefix = _line.split('\t',1)[0].startswith('chr')
                        break
            except Exception:
                has_chr_prefix = True  

            score_file = os.path.join(scores_dir, f"{pgs_id}.{subpop}.chr{chr_num}.score.tsv")
            opener = gzip.open if weights.endswith((".gz", ".gs")) else open
            with opener(weights, "rt") as fin, open(score_file, "w") as fout:
                header = fin.readline().rstrip("\n").split("\t")
                idx = {h:i for i,h in enumerate(header)}
                for req in ("chr","pos","ref","alt","effect_allele","beta"):
                    if req not in idx:
                        raise SystemExit(f"[ERROR] Missing '{req}' in {weights}")
                fout.write("ID\tA1\tBETA\n")
                for line in fin:
                    t = line.rstrip("\n").split("\t")
                    w_chr = str(t[idx['chr']]).lstrip('c').lstrip('h').lstrip('r') if t[idx['chr']].startswith('chr') else str(t[idx['chr']])
                    if chr_num and w_chr != chr_num:
                        continue
                    id_chr = ("chr" + chr_num) if has_chr_prefix else chr_num
                    vid = f"{id_chr}:{t[idx['pos']]}:{t[idx['ref']]}:{t[idx['alt']]}"
                    fout.write(f"{vid}\t{t[idx['effect_allele']]}\t{t[idx['beta']]}\n")

            out_pref = os.path.join(scores_dir, f"{pgs_id}.{subpop}." + os.path.basename(pref))
            run([
                "plink2",
                "--pfile", pref,
                "--score", score_file, "1", "2", "3", "header-read", "no-mean-imputation", "list-variants",
                "--threads", str(args.threads),
                "--out", out_pref
            ])
            sscore = out_pref + ".sscore"
            svars  = out_pref + ".sscore.vars"
            if not os.path.exists(sscore):
                continue
            df = pd.read_csv(sscore, sep=r"\s+")
            df.columns = [c.lstrip("#") for c in df.columns]
            score_col = detect_score_col(df)
            df = df[["IID", score_col]].rename(columns={score_col: "score_sum"})
            if acc is None:
                acc = df.copy()
            else:
                acc = acc.merge(df, on="IID", how="outer").fillna(0)
                if "score_sum_x" in acc.columns and "score_sum_y" in acc.columns:
                    acc["score_sum"] = acc["score_sum_x"].astype(float) + acc["score_sum_y"].astype(float)
                    acc = acc.drop(columns=["score_sum_x","score_sum_y"])        
            total_vars += count_lines_excluding_header(svars)

        if acc is None or acc.empty:
            print(f"[WARN] No per-chromosome scores for {pgs_id} {subpop}; skipping.", file=sys.stderr)
            continue

        j = acc.merge(lab, left_on="IID", right_on="iid", how="left")
        sub = j[j["subpop"] == subpop]
        n_samples = sub.shape[0]
        if n_samples == 0:
            print(f"[WARN] No samples for {subpop}; skipping {pgs_id}.", file=sys.stderr)
            continue

        mean = float(sub["score_sum"].astype(float).mean())
        sd   = float(sub["score_sum"].astype(float).std(ddof=0))
        n_variants = int(total_vars)

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

    if args.append and os.path.exists(args.out_tsv) and os.path.getsize(args.out_tsv) > 0:
        old = pd.read_csv(args.out_tsv, sep="\t", dtype=str)
        old.columns = [c.strip().lower() for c in old.columns]
        combo = pd.concat([old, out_df_new], ignore_index=True)
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
