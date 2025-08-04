import fs from "fs";
import path from "path";
import { execSync } from "child_process";
import { USER_DATA_DIR, APP_ROOT } from "./paths.js";

export const ensureDir = (p) => fs.mkdirSync(p, { recursive: true });

export const rel = (p) =>
  path.relative(USER_DATA_DIR, p).startsWith("..")
    ? path.relative(APP_ROOT, p)
    : path.relative(USER_DATA_DIR, p);

export const downloadFile = (dest, url, desc = "") => {
  if (fs.existsSync(dest) && fs.statSync(dest).size > 0) {
    console.log("✔", rel(dest), "exists");
    return true;
  }
  ensureDir(path.dirname(dest));
  for (let attempt = 1; attempt <= 2; attempt++) {
    console.log(
      `⬇︎ DOWNLOADING${attempt > 1 ? " (retry)" : ""}: ${
        desc || path.basename(dest)
      }: ${url}`
    );
    try {
      execSync(`wget -q -O "${dest}" "${url}"`);
      if (fs.statSync(dest).size === 0) throw new Error("empty file");
      if (path.basename(dest) === "eagle") fs.chmodSync(dest, 0o755);
      console.log("✓ Downloaded", rel(dest), "\n");
      return true;
    } catch (e) {
      try {
        fs.unlinkSync(dest);
      } catch {}
      if (attempt === 2) {
        console.error("✖ ERROR:", url, e.message);
        return false;
      }
    }
  }
};

export const gunzipTo = (srcGz, finalName) => {
  if (fs.existsSync(finalName)) {
    console.log("✔", rel(finalName), "exists");
    return true;
  }
  if (!fs.existsSync(srcGz)) return false;
  console.log("→ gunzip:", rel(srcGz));
  try {
    execSync(`gunzip -f "${srcGz}"`, { stdio: "inherit" });
  } catch (err) {
    console.warn("! gunzip warning:", err.message);
  }
  const out = srcGz.replace(/\.gz$/, "");
  if (!fs.existsSync(out)) {
    console.error("✖ ERROR: gunzip failed:", rel(srcGz));
    return false;
  }
  if (out !== finalName) fs.renameSync(out, finalName);
  console.log("✓ Decompressed", rel(finalName));
  return true;
};
