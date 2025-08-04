import fs from "fs";
import path from "path";
import { GENOME_DATA_DIR, BEAGLE_JAR, BREF3_JAR } from "./paths.js";
import { rel } from "./helpers.js";

const mustFiles = [
  path.join(GENOME_DATA_DIR, "fasta", "Homo_sapiens_assembly38.fasta"),
  path.join(GENOME_DATA_DIR, "fasta", "Homo_sapiens_assembly19.fasta"),
  path.join(GENOME_DATA_DIR, "chain", "hg19ToHg38.over.chain.gz"),
  path.join(GENOME_DATA_DIR, "jars", BREF3_JAR),
  path.join(GENOME_DATA_DIR, "jars", BEAGLE_JAR),
  path.join(GENOME_DATA_DIR, "bin", "eagle"),
  path.join(GENOME_DATA_DIR, "scripts", "addchr.txt"),
  path.join(GENOME_DATA_DIR, "scripts", "dropchr.txt"),
];

let ok = true;
mustFiles.forEach((f) => {
  if (!fs.existsSync(f)) {
    console.error("✖ ERROR: Missing", rel(f));
    ok = false;
  }
});
if (!ok) {
  console.error("✖ ERROR: Critical files missing, abort");
  process.exit(1);
}

console.log("\n🧹 cleaning temp_downloads");
try {
  fs.rmSync(tempDir, { recursive: true, force: true });
} catch {}
