#!/usr/bin/env node

import fs from "fs";
import path from "path";
import { execSync } from "child_process";
import { GENOME_DATA_DIR, BEAGLE_JAR, BREF3_JAR } from "./paths.js";
import { ensureDir, rel, downloadFile } from "./helpers.js";

// ────────────── create bcf+bref3 for b38 1000G panels ──────────────
console.log("\n📦 1000 Genomes reference panels (create bcf+bref3 for b38)");
const tempDir = path.join(GENOME_DATA_DIR, "temp_downloads");
ensureDir(tempDir);
const chrList = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), "X"];

const bcftools = (cmd) => execSync(`bcftools ${cmd}`, { stdio: "inherit" });
const makeBref3 = (vcf, outPrefix) => {
  const jar = path.join(GENOME_DATA_DIR, "jars", BREF3_JAR);
  const bref = `${outPrefix}.bref3`;
  let fd;
  try {
    ensureDir(path.dirname(bref));
    fd = fs.openSync(bref, "w");

    execSync(`java -jar "${jar}" "${vcf}"`, {
      stdio: ["inherit", fd, "inherit"],
    });

    return fs.existsSync(bref) && fs.statSync(bref).size > 100_000;
  } catch (e) {
    console.error("✖ ERROR: bref3 gen failed.", e.message);
    return false;
  } finally {
    if (fd !== undefined)
      try {
        fs.closeSync(fd);
      } catch {}
  }
};

for (const chr of chrList) {
  console.log(`\n🧬 chr${chr}`);

  const vcfUrl = `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`;
  const vcfGz = path.join(
    tempDir,
    `1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`
  );
  const vcfTbi = `${vcfGz}.tbi`;
  const bcfOut = path.join(
    GENOME_DATA_DIR,
    "ref_bcfs_b38",
    `1000GP_chr${chr}.bcf`
  );
  const bref3Out = path.join(
    GENOME_DATA_DIR,
    "ref_brefs_b38",
    `1000GP_chr${chr}.bref3`
  );

  const needBcf = !(
    fs.existsSync(bcfOut) &&
    fs.existsSync(`${bcfOut}.csi`) &&
    fs.statSync(bcfOut).size > 1_000_000
  );
  const needBref = !(
    fs.existsSync(bref3Out) && fs.statSync(bref3Out).size > 100_000
  );

  if (!needBcf && !needBref) {
    console.log("✓ Already processed");
    continue;
  }

  // Download if needed
  // Download if needed
  if (!fs.existsSync(vcfGz)) {
    if (!downloadFile(vcfGz, vcfUrl, `1000G chr${chr} VCF`)) continue;
  }

  // ensure .tbi index is present
  if (!fs.existsSync(vcfTbi)) {
    downloadFile(vcfTbi, `${vcfUrl}.tbi`, `1000G chr${chr} index`);
  }

  // check for corruption
  try {
    execSync(`gzip -t "${vcfGz}"`);
  } catch {
    console.error("✖ ERROR: corrupted VCF, re-downloading");
    try {
      fs.unlinkSync(vcfGz);
    } catch {}
    try {
      fs.unlinkSync(vcfTbi);
    } catch {}

    if (!downloadFile(vcfGz, vcfUrl, `1000G chr${chr} VCF (retry)`)) continue;
    downloadFile(vcfTbi, `${vcfUrl}.tbi`, `1000G chr${chr} index (retry)`);
  }

  // Build a SNP+INDEL-only VCF once if we need either output
  if (needBcf || needBref) {
    const snvindelVcf = path.join(tempDir, `snvindel.chr${chr}.vcf.gz`);

    // Keep only SNPs and INDELs (filters out symbolic SVs like <DEL>, <INS>, etc.)
    bcftools(
      `view -i 'TYPE="snp" || TYPE="indel"' -Oz -o "${snvindelVcf}" "${vcfGz}"`
    );
    bcftools(`index -f "${snvindelVcf}"`);

    // Convert to BCF (from SNP/INDEL-only VCF)
    if (needBcf) {
      console.log("🔄 VCF -> BCF (SNVs+INDELs only)");
      ensureDir(path.dirname(bcfOut));
      bcftools(`view -Ob -o "${bcfOut}" "${snvindelVcf}"`);
      bcftools(`index "${bcfOut}"`);
    }

    // Convert to BREF3 (from SNP/INDEL-only VCF)
    if (needBref) {
      console.log("🔄 VCF -> bref3 (SNVs+INDELs only)");

      // makeBref3 expects an outPrefix (no .bref3 extension)
      const bref3Prefix = bref3Out.replace(/\.bref3$/, "");

      if (!makeBref3(snvindelVcf, bref3Prefix)) {
        console.error(`✖ ERROR: failed to create bref3 for chr${chr}`);
      }
    }

    // Clean temp files (both the original and filtered VCFs)
    // try {
    //   fs.unlinkSync(vcfGz);
    //   fs.unlinkSync(vcfTbi);
    //   fs.unlinkSync(snvindelVcf);
    //   try {
    //     fs.unlinkSync(`${snvindelVcf}.tbi`);
    //   } catch {}
    // } catch {}
  } else {
    // Nothing to do for this chr; just clean the original downloads
    // try {
    //   fs.unlinkSync(vcfGz);
    //   fs.unlinkSync(vcfTbi);
    // } catch {}
  }
}

// ───────────────────── final check + clean up ─────────────────────

// Per-chr refs
const missingChr = chrList.filter((chr) => {
  const bcf = path.join(
    GENOME_DATA_DIR,
    "ref_bcfs_b38",
    `1000GP_chr${chr}.bcf`
  );
  const bref = path.join(
    GENOME_DATA_DIR,
    "ref_brefs_b38",
    `1000GP_chr${chr}.bref3`
  );
  return !fs.existsSync(bcf) || !fs.existsSync(bref);
});
if (missingChr.length) {
  console.error("✖ ERROR: Incomplete 1000G set.", missingChr.join(", "));
  process.exit(1);
}

console.log("\nCOMPLETE! Reference stack ready for run_imputation.sh");
console.log(`📁 Path: ${GENOME_DATA_DIR}`);
// console.log(
//   "You can now run: ./run_imputation.sh target_genomes/your_sample.txt"
// );
