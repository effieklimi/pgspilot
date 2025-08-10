#!/bin/sh
set -e

# Usage: bash run_pipeline.sh PGS000123 [PGS000456 ...]
if [ $# -lt 1 ]; then
  echo "Error: Please provide one or more PGS IDs."
  echo "Usage: bash run_pipeline.sh PGS000123 [PGS000456 ...]"
  exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)

# Pass all IDs through to the real runner
bash "${SCRIPT_DIR}/pipeline/add_pgs.sh" "$@"
