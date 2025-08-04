#!/usr/bin/env node

import { existsSync } from "fs";
import path from "path";
import os from "os";
import { fileURLToPath } from "url";

export const APP_NAME = "Iodine";

export const BEAGLE_JAR = "beagle.27Feb25.75f.jar";
export const BREF3_JAR = "bref3.27Feb25.75f.jar";
export const EAGLE_BIN_GZ = "Eagle_v2.4.1.tar.gz";

export const getUserDataDir = () => {
  if (process.env.DOCKER_CONTAINER || existsSync("/.dockerenv")) {
    return "/app/genome_data";
  }
  const platform = os.platform();
  const home = os.homedir();
  switch (platform) {
    case "darwin":
      return path.join(home, "Library", "Application Support", "iodine");
    case "win32":
      return path.join(
        process.env.APPDATA || path.join(home, "AppData", "Roaming"),
        "iodine"
      );
    case "linux":
      return path.join(
        process.env.XDG_DATA_HOME || path.join(home, ".local", "share"),
        "iodine"
      );
    default:
      return path.join(home, ".iodine");
  }
};

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

export const IODINE_DATA_DIR = getUserDataDir();
export const USER_DATA_DIR = path.join(IODINE_DATA_DIR, "user_data");
export const GENOME_DATA_DIR = path.join(IODINE_DATA_DIR, "genome_data");
export const TRAITS_DIR = path.join(IODINE_DATA_DIR, "traits");
export const UTILS_DIR = path.join(IODINE_DATA_DIR, "utils");
export const APP_ROOT = process.env.APP_ROOT || path.join(__dirname, "../..");
