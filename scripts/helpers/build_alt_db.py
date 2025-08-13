#!/usr/bin/env python3
"""
build_alt_db.py - Build SQLite database of ALT alleles from 1000 Genomes reference files

Usage: python3 build_alt_db.py [--ref-dir PATH]
"""

import sqlite3
import gzip
import os
import glob
import argparse
import subprocess

def extract_variants_from_bcf(bcf_file):
    """Extract variants from BCF file using bcftools"""
    try:
        # Use bcftools to extract variant info
        cmd = ['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', bcf_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip().split('\n')
    except subprocess.CalledProcessError as e:
        print(f"❌ Error reading {bcf_file}: {e}")
        return []

def build_alt_database(ref_dir=None):
    """Build SQLite database from 1000G reference files"""
    
    if ref_dir is None:
        # Use the same path structure as the main imputation script
        root_dir = os.getcwd()
        ref_dir = os.path.join(root_dir, "genome_data", "ref_bcfs_b38")  # Eagle BCF files
    
    if not os.path.exists(ref_dir):
        print(f"❌ Reference directory not found: {ref_dir}")
        print("💡 Make sure reference files are downloaded with your setup script")
        return
    
    db_path = "genome_data/alt_alleles.db"
    os.makedirs("genome_data", exist_ok=True)
    
    # Create database
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    
    # Create table
    cur.execute("""
        CREATE TABLE IF NOT EXISTS alt_alleles (
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            PRIMARY KEY (chrom, pos, ref)
        )
    """)
    
    # Look for BCF files (same pattern as your imputation script)
    file_pattern = f"{ref_dir}/1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.bcf"
    ref_files = glob.glob(file_pattern)
    
    if not ref_files:
        print(f"❌ No BCF files found in {ref_dir}")
        print(f"💡 Expected pattern: 1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.bcf")
        return
    
    print(f"📊 Processing {len(ref_files)} reference files...")
    
    total_variants = 0
    for ref_file in sorted(ref_files):
        # Extract chromosome from filename (1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.bcf -> 1)
        basename = os.path.basename(ref_file)
        chr_part = basename.replace('1kGP_high_coverage_Illumina.chr', '').replace('.filtered.SNV_INDEL_SV_phased_panel.bcf', '')
            
        print(f"  Processing chr{chr_part}...")
        
        # Extract variants from BCF using bcftools
        variant_lines = extract_variants_from_bcf(ref_file)
        
        # Process variants
        for line in variant_lines:
            if not line.strip():
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
                
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[2]
            alt = fields[3]
            
            # Normalize chromosome names (remove chr prefix if present)
            chrom = chrom.replace('chr', '')
            
            # Skip complex variants and multi-allelic sites
            if ',' in alt or '*' in alt or len(ref) > 50 or len(alt) > 50:
                continue
            
            # Insert with both chr-prefixed and non-prefixed versions
            for chr_format in [chrom, f"chr{chrom}"]:
                cur.execute(
                    "INSERT OR IGNORE INTO alt_alleles (chrom, pos, ref, alt) VALUES (?, ?, ?, ?)",
                    (chr_format, pos, ref, alt)
                )
            
            total_variants += 1
            
            if total_variants % 100000 == 0:
                print(f"    Processed {total_variants:,} variants...")
                conn.commit()  # Periodic commits
    
    # Create index for faster lookups
    print("📊 Creating database index...")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_lookup ON alt_alleles (chrom, pos, ref)")
    
    conn.commit()
    conn.close()
    
    print(f"✅ ALT alleles database created: {db_path}")
    print(f"📊 Total variants: {total_variants:,}")
    
    # Test the database
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM alt_alleles")
    count = cur.fetchone()[0]
    print(f"🔍 Database contains {count:,} entries")
    conn.close()

def main():
    parser = argparse.ArgumentParser(description="Build ALT alleles database")
    parser.add_argument("--ref-dir", help="Directory containing reference VCF/BCF files")
    args = parser.parse_args()
    
    build_alt_database(args.ref_dir)

if __name__ == "__main__":
    main()