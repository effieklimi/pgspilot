#!/usr/bin/env bash
set -euo pipefail

if [ $# -eq 0 ]; then
  echo "Error: Please provide the path to your input file."
  echo "Usage: bash run_pipeline.sh /path/to/your/file.txt"
  exit 1
fi

# Find the directory of this script, so it can be run from anywhere.
SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)

# 1) Run the PCA setup (in Docker via setup.sh)
bash "${SCRIPT_DIR}/pipeline/qc_genome.sh" "$@"