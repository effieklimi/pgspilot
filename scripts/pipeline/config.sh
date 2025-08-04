# Shared settings for all scripts. 

THREADS=${THREADS:-4}
ROOT_DIR=$(pwd)

BASE=$(basename "$IN_TXT")
STEM=${BASE%%.*}

INPUT_BUILD="${INPUT_BUILD:-auto}"   # override with INPUT_BUILD=38 ./impute_23andme.sh file.txt
TARGET_BUILD=38

# References
FASTA_37="${ROOT_DIR}/genome_data/fasta/Homo_sapiens_assembly19.fasta"
FASTA_38="${ROOT_DIR}/genome_data/fasta/Homo_sapiens_assembly38.fasta"
CHAIN="${ROOT_DIR}/genome_data/chain/hg19ToHg38.over.chain.gz"

# 1000G / reference cohort (for calibration, not for PCA)
ONEKG_VCF_PATTERN="${ROOT_DIR}/genome_data/1000g/1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"  # {chr} replaced by 1..22
ONEKG_LABELS="${ROOT_DIR}/genome_data/1000g/integrated_call_samples_v3.20130502.ALL.panel"         # i don't know what this is! help!

# PCA model output directory (built once via fit_pca_1kg.py)
PCA_DIR="${ROOT_DIR}/pca_model"

# Project directories
WEIGHTS_RAW_DIR="${ROOT_DIR}/weights_raw"   # source PGS files (from PGS Catalog)
WEIGHTS_HM_DIR="${ROOT_DIR}/weights_hm"     # harmonized weights (output)
TRAITS_DIR="${ROOT_DIR}/traits"             # frozen include lists (output)
CALIB_DIR="${ROOT_DIR}/calibration"         # trait calibrations (output)
USERS_DIR="${ROOT_DIR}/users"               # per-user outputs

# Scorable & optional mask
SCORABLE_SITES="scorable_sites.b38.tsv"   # produced via your export script
SITE_MASK=""                              # optional; leave empty if none
MASK_MIN_DR2="0.50"                       # used if SITE_MASK is provided

# Build knobs
INCLUDE_CHROMS="1-22"   # pilot: autosomes only
DROP_MHC="--drop-mhc"   # leave empty "" to keep MHC
QUANTILES=1001
MIN_VARIANTS=20

# Python entry points (if you want to pin a venv)
PYTHON="python"

# Scripts (paths relative to repo root)
WEIGHTS_HARMONIZER="${ROOT_DIR}/scripts/analyses/weights_harmonize_b38.py"
INCLUDE_BUILDER="${ROOT_DIR}/scripts/analyses/build_include_lists.py"
CALIB_BUILDER="${ROOT_DIR}/scripts/analyses/build_calibration.py"
PCA_FITTER="${ROOT_DIR}/scripts/analyses/fit_pca_1kg.py"
SCORER="${ROOT_DIR}/scripts/analyses/score_user_and_report.py"

# Imputation script (your existing pipeline)
IMPUTE_SH="${ROOT_DIR}/scripts/analyses/impute.sh"   # path to your imputation script
THREADS=${THREADS:-4}
