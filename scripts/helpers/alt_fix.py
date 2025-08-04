#!/usr/bin/env python3
"""
alt_fix.py  – replace missing ALT alleles ('.') in a VCF.gz using a pre-built
SQLite lookup table.

Project layout assumed
├─ static_files/
│   ├─ alt_alleles.db            (the lookup DB you built once)
│   └─ scripts/
│       └─ alt_fix.py            (this file)
└─ <sample>_results/
    └─ <sample>.vcf.gz           (VCF you want to fix)
        <sample>.vcf.gz.tbi
"""

import argparse, gzip, os, sqlite3, sys, time

def replace_alts(vcf_in, vcf_out, db_path, preview=5):
    if not os.path.exists(db_path):
        sys.exit(f"[ERROR] DB not found: {db_path}")

    conn   = sqlite3.connect(db_path)
    cur    = conn.cursor()
    t0     = time.time()
    total  = rep = 0

    with gzip.open(vcf_in,  "rt") as fin, gzip.open(vcf_out, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line); continue
            total += 1
            f      = line.rstrip("\n").split("\t")
            if f[4] == ".":                                    # ALT missing
                cur.execute("SELECT alt FROM alt_alleles WHERE chrom=? AND pos=? AND ref=?",
                            (f[0], int(f[1]), f[3]))
                row = cur.fetchone()
                if row:
                    f[4] = row[0]; rep += 1
                    if rep <= preview:
                        print(f"  • {f[0]}:{f[1]}  {f[3]}→{f[4]}")
            fout.write("\t".join(f) + "\n")

    conn.close()
    dt = time.time() - t0
    print(f"\n✓ {rep:,} of {total:,} records patched  ({dt:.1f}s)")
    print(f"✓ Output written: {vcf_out}")

def main():
    p = argparse.ArgumentParser(
        description="Fill '.' ALT alleles in VCF.gz from SQLite lookup")
    p.add_argument("vcf",  help="input VCF.gz (must be indexed)")
    p.add_argument("--db", default="genome_data/alt_alleles.db",
                   help="SQLite DB (default genome_data/alt_alleles.db)")
    args = p.parse_args()

    if not args.vcf.endswith(".vcf.gz"):
        sys.exit("Input must be bgzipped *.vcf.gz")

    vcf_out = args.vcf.replace(".vcf.gz", ".alt.vcf.gz")
    replace_alts(args.vcf, vcf_out, args.db)

if __name__ == "__main__":
    main()
