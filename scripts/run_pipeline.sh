#!/bin/sh
set -e # Exit immediately if a command exits with a non-zero status.

# This is the main entry point for the client.
# It avoids the need for them to run `chmod` on any files.

if [ $# -eq 0 ]; then
  echo "Error: Please provide the path to your input file."
  echo "Usage: bash run_pipeline.sh /path/to/your/file.txt"
  exit 1
fi

# We find the directory of this script, so it can be run from anywhere.
SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)

# We call the *real* wrapper script using bash, passing along all arguments.
bash "${SCRIPT_DIR}/pipeline/user.sh" "$@"