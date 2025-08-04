#!/usr/bin/env node

import { promises as fs } from "fs";
import fsSync from "fs";
import path from "path";
import { execSync } from "child_process";
import {
  IODINE_DATA_DIR,
  USER_DATA_DIR,
  GENOME_DATA_DIR,
  TRAITS_DIR,
  UTILS_DIR,
  BEAGLE_JAR,
  BREF3_JAR,
  EAGLE_BIN_GZ,
  APP_NAME,
} from "./paths.js";
import { ensureDir, rel, downloadFile, gunzipTo } from "./helpers.js";

// add a preflight check for bcftools, java, tabix/gzip, unzip, and gunzip

async function setupDirectories() {
  try {
    console.log(`Setting up ${APP_NAME} data directories...`);

    // Create directories
    await fs.mkdir(IODINE_DATA_DIR, { recursive: true });
    await fs.mkdir(USER_DATA_DIR, { recursive: true });
    await fs.mkdir(GENOME_DATA_DIR, { recursive: true });
    await fs.mkdir(TRAITS_DIR, { recursive: true });
    await fs.mkdir(UTILS_DIR, { recursive: true });

    const userDataPath = path.join(USER_DATA_DIR, "user-data.json");
    try {
      await fs.access(userDataPath);
      console.log("\n✔ User data file already exists");
    } catch {
      const initialData = {
        profile: {
          name: "",
          age: "",
          race: "",
          ethnicity: "",
          healthConditions: "",
        },
        genomeFile: null,
      };
      await fs.writeFile(userDataPath, JSON.stringify(initialData, null, 2));
      console.log("✓ Created initial user data file");
    }

    console.log(
      `\n${APP_NAME} data directories set up successfully at: ${IODINE_DATA_DIR}`
    );
  } catch (error) {
    console.error("✖ ERROR: failed to set up directories:", error);
    process.exit(1);
  }
}
// Run the setup
await setupDirectories();

//────────────────── small stuff - chain, tools ──────────────────
console.log("\n📦 CHAIN & TOOLS");
ensureDir(path.join(GENOME_DATA_DIR, "eagle_maps_b38"));
ensureDir(path.join(GENOME_DATA_DIR, "beagle_maps_b38"));
ensureDir(path.join(GENOME_DATA_DIR, "temp_downloads"));
ensureDir(path.join(GENOME_DATA_DIR, "bin"));

const eagleTar = path.join(GENOME_DATA_DIR, "temp_downloads", EAGLE_BIN_GZ);

downloadFile(
  path.join(GENOME_DATA_DIR, "chain", "hg19ToHg38.over.chain.gz"),
  "https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz",
  "hg19→hg38 chain"
);
downloadFile(
  path.join(GENOME_DATA_DIR, "jars", BEAGLE_JAR),
  `https://faculty.washington.edu/browning/beagle/${BEAGLE_JAR}`,
  "Beagle 5.5 jar"
);
downloadFile(
  path.join(GENOME_DATA_DIR, "jars", BREF3_JAR),
  `https://faculty.washington.edu/browning/beagle/${BREF3_JAR}`,
  "Bref3 jar"
);
downloadFile(
  path.join(GENOME_DATA_DIR, "temp_downloads", EAGLE_BIN_GZ),
  `https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/${EAGLE_BIN_GZ}`,
  "Eagle binary"
);

try {
  execSync(
    `tar -xzf "${eagleTar}" -C "${path.join(
      GENOME_DATA_DIR,
      "bin"
    )}" --strip-components=1 ${EAGLE_BIN_GZ.replace(/\.tar\.gz$/, "")}/eagle`,
    { stdio: "inherit" }
  );
  fsSync.chmodSync(path.join(GENOME_DATA_DIR, "bin", "eagle"), 0o755);
  console.log("✓ Extracted Eagle binary");
} catch (e) {
  console.error("✖ ERROR: extracting Eagle binary:", e.message);
}

//────────────────── FASTAs ──────────────────
console.log("\n📦 FASTA (GRCh37 & GRCh38)");

const fastaDir = path.join(GENOME_DATA_DIR, "fasta");
const f37gz = path.join(fastaDir, "Homo_sapiens_assembly19.fasta.gz");
const f37 = path.join(fastaDir, "Homo_sapiens_assembly19.fasta");
const f38gz = path.join(fastaDir, "Homo_sapiens_assembly38.fasta.gz");
const f38 = path.join(fastaDir, "Homo_sapiens_assembly38.fasta");

if (fsSync.existsSync(f37)) {
  console.log("✔ GRCh37 FASTA already unzipped");
} else {
  if (fsSync.existsSync(f37gz)) {
    console.log("✔ GRCh37 FASTA .gz found, skipping download");
  } else {
    downloadFile(
      f37gz,
      "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz",
      "GRCh37 FASTA"
    );
  }
  gunzipTo(f37gz, f37);
}

if (fsSync.existsSync(f38)) {
  console.log("✔ GRCh38 FASTA already unzipped");
} else {
  if (fsSync.existsSync(f38gz)) {
    console.log("✔ GRCh38 FASTA .gz found, skipping download");
  } else {
    downloadFile(
      f38gz,
      "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
      "GRCh38 FASTA"
    );
  }
  gunzipTo(f38gz, f38);
}

//────────────────── maps ──────────────────
console.log("\n📦 GENETIC MAPS (b38)");
const eagleMap38 = path.join(
  GENOME_DATA_DIR,
  "eagle_maps_b38",
  "genetic_map_hg38_withX.txt.gz"
);
downloadFile(
  eagleMap38,
  "https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz",
  "Eagle map hg38"
);
const fullMap = path.join(
  GENOME_DATA_DIR,
  "eagle_maps_b38",
  "genetic_map_hg38_withX.txt.gz"
);

// Create per-chr symlinks (best), or copies (fallback)
for (const c of [...Array.from({ length: 22 }, (_, i) => String(i + 1)), "X"]) {
  const perChr = path.join(
    GENOME_DATA_DIR,
    "eagle_maps_b38",
    `eagle_chr${c}_b38.map.gz`
  );
  try {
    // clean any previous stub
    if (!fsSync.existsSync(perChr)) {
      try {
        fsSync.symlinkSync(fullMap, perChr);
      } catch {
        fsSync.copyFileSync(fullMap, perChr);
      }
    }
  } catch {
    // fallback to a real copy if symlinks are disallowed
    fsSync.copyFileSync(fullMap, perChr);
  }
}

const beagleZip38 = path.join(
  GENOME_DATA_DIR,
  "temp_downloads",
  "beagle_maps_b38.zip"
);
const beagleDir = path.join(GENOME_DATA_DIR, "beagle_maps_b38");
downloadFile(
  beagleZip38,
  "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip",
  "Beagle maps hg38 zip"
);
ensureDir(beagleDir);
try {
  const expected = path.join(beagleDir, "plink.chr22.GRCh38.map"); // top-level
  if (!fsSync.existsSync(expected)) {
    execSync(`unzip -o "${beagleZip38}" -d "${beagleDir}"`, {
      stdio: "inherit",
    });
    console.log("✓ Extracted Beagle maps (b38)");
  } else {
    console.log("✔ Beagle maps already extracted");
  }
} catch (e) {
  console.error("✖ ERROR: unzip beagle maps:", e.message);
}

//──────────── 1000G subpopulations pannel ────────────
console.log("\n📦 1000G SUBPOPULATIONS PANEL");
const subpopDir = path.join(GENOME_DATA_DIR, "1000G_subpopulations_panel");
ensureDir(subpopDir);
downloadFile(
  path.join(subpopDir, "integrated_call_samples_v3.20130502.ALL.panel"),
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
  "1000G subpopulations panel"
);

//────────────────── chr rename maps ──────────────────
console.log("\n📄 CHROMOSOME RENAME MAPS");
const scriptsDir = path.join(GENOME_DATA_DIR, "scripts");
ensureDir(scriptsDir);

// Helper: write only if missing or content differs
const writeIfChanged = (p, data) => {
  if (fsSync.existsSync(p)) {
    const cur = fsSync.readFileSync(p, "utf8");
    if (cur === data) {
      console.log("✔", rel(p), "up-to-date");
      return;
    }
    fsSync.writeFileSync(p, data);
    console.log("🔄 Updated", rel(p));
    return;
  }
  fsSync.writeFileSync(p, data);
  console.log("✓ Wrote", rel(p));
};

const autosomes = Array.from({ length: 22 }, (_, i) => String(i + 1));

// addchr.txt: no-chr → chr (e.g., 1 → chr1, X → chrX, MT/M → chrM)
const addChrContent =
  [
    ...autosomes.map((n) => `${n}\tchr${n}`),
    "X\tchrX",
    "Y\tchrY",
    "MT\tchrM",
    "M\tchrM",
  ].join("\n") + "\n";

// dropchr.txt: chr → no-chr (e.g., chr1 → 1, chrX → X, chrM/chrMT → MT)
const dropChrContent =
  [
    ...autosomes.map((n) => `chr${n}\t${n}`),
    "chrX\tX",
    "chrY\tY",
    "chrM\tMT",
    "chrMT\tMT",
  ].join("\n") + "\n";

writeIfChanged(path.join(scriptsDir, "addchr.txt"), addChrContent);
writeIfChanged(path.join(scriptsDir, "dropchr.txt"), dropChrContent);

for (const c of [...Array.from({ length: 22 }, (_, i) => String(i + 1)), "X"]) {
  const src = path.join(beagleDir, `plink.chr${c}.GRCh38.map`);
  const dst = path.join(beagleDir, `beagle_chr${c}_b38.map`);
  if (!fsSync.existsSync(src)) {
    console.warn(`! Missing ${path.basename(src)} after unzip`);
    continue;
  }
  if (!fsSync.existsSync(dst)) fsSync.copyFileSync(src, dst);
}
