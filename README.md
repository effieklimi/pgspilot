# Pilot pipelines for imputation and PGS scoring

## Run once at the start:

1. **Generate a docker image:**

```
cd pgspilot
docker build --no-cache -t pgspilot .
```

2. **Perform initial setup:**

```
bash scripts/run_setup.sh
```

**This will (idempotently):**

- Download all required files (using the 1000 Genome project expansion)
- Build bcf and bref3 genome references
- Generate the alt_alleles.db needed for imputation
- Run a PCA for subpopulation assignment

## Run every time you need to add a new trait to your PGS registry:

```
# An example for "insomnia"
bash scripts/run_add_pgs.sh PGS002149
```

**Key output:** a folder called `/pgs/weights` in the project root with everything you need to score users later.

**This will:**

- Fetch the already harmonized PGS weights file from the PGS catalogue
- Create weights files per subpopulation
- Append the PGS ID, along with other metadata, to a global registry of all PGSs added so far, stored in `/pgs/weights/harmonized`
- Standardize it per subpopulation (needed for subpopulation-based ranking)
- Keep a stats summary of the subpopulation standardization, per PGS ID

## Run every time you need to score a user:

```
bash scripts/run_qc_genome.sh <path/to/23andMe/file>.txt
bash scripts/run_user.sh <path/to/23andMe/file>.txt
```

**The key output:** a file called `PGSXXXX.SUB.normalized.tsv` per PGS ID.

- **PGSXXXX:** PGS ID
- **SUB:** the subpopulation the user was assigned using the PCA (one of EUR, AMR, SAS, EAS, AFR)

**These will:**

- Do quality control (QC) on the given genome file. `run_qc_genome` writes:

  1. A single text report with QC/metadata and a final status line (PASS/FAIL):

  ```
  # ... QC details
  ---------------------------------
  Flags:          (none)
  QC STATUS:      PASS
  =================================


  # example of a failed QC run:
  Flags:
  ---------------------------------
  - Variant count (1415218) outside expected 23andMe v3/v4/v5 windows; file may be imputed, merged, or from a different vendor
  QC STATUS:      FAIL
  =================================
  ```

  2. Filename: {STEM}\_initial_qc.txt where {STEM} is your input filename without .txt/.gz.

- Run imputation on the user's genome, and generate a directory called `/users/{STEM}` inside where all subsequent results will be written for the user
- Imputed genome: `/users/{STEM}/{STEM}_imputed_all.vcf.gz` and its corresponding indexed file `.tbi`
- Imputation quality control files: in `/users/{STEM}/imputation_qc_reports`
- User's plink2 files: in `/users/{STEM}/pfiles`
- At this stage, the user's imputed genome has been filtered to keep only variants with high imputation quality (i.e. DR2/INFO >=0.8). Important for PGS scoring.
- Subpopulation assigned using PCA: `/users/{STEM}/ancestry.tsv`
- User's PGS scores: `/users/{STEM}/pgs_scores`. They are organised in directories, each corresponding to one PGS score:

```
pgs_scores/
├── PGS000300/ # example of a PGS weights file that has NO OVERLAP with the user's variants
│   ├── PGS000300.AMR.normalized.tsv
│   ├── PGS000300.AMR.skipped.tsv # `.skipped` folder is written when there is no overlap
│   ├── PGS000300.AMR.sscore
│   └── PGS000300.AMR.sscore.vars # contain's the user's exact variants that were included in the PGS weights file (should be empty here)
└── PGS002149/ # example of a PGS weights file that has overlap with the user's variants
    ├── PGS002149.AMR.log
    ├── PGS002149.AMR.normalized.tsv # this contains both RAW score and SUBPOPULATION-NORMALISED score
    ├── PGS002149.AMR.sscore
    └── PGS002149.AMR.sscore.vars # this one should not be empty
```

- SUBPOPULATION-NORMALISED score is a standard deviation - where the user is placed within their assigned subpopulation (i.e. how much higher or lower risk VS their subpopulation's average)

## Redo PCA locally

**Important!** PCA is used to assign each user to a subpopulation out of EUR, AMR, SAS, EAS, AFR, for PGS scoring purposes. It is **not** an ancestry, and it should not be reported as such to the user. Instead, it assigns the user into the subpopulation he/she is most similar to, specifically within the 1000 Genomes expanded panel used here. This is expected to positivelly correlate with the ancestry of the user, but it is not a robust prediction of their ancestry.

In the current pilot implementation, the PCA was run with 6 principal components for speed. The analyses will run using those PCA results, but the scores obtained with it will not be optimal. Six PCs is not enough here.

**Things to keep in mind and steps that need to be taken:**

- Then, delete everything inside the `pca_model` folder
- Re-run the setup script with `bash scripts/run_setup.sh`. This should re-run the PCA
- This might take some time to run
- I have already configured `scripts/pipeline/setup.sh` (the "project setup" script) to run with 12 principal components next time by default by passing it the flag: `--pcs 12`
- PGS IDs will need to be re-run after the PCA is re-run. Delete all old PGS-related directories in `/pgs/weights`, and re-run the `bash scripts/run_add_pgs.sh <PGS_ID>`

The new PCA will need to be inspected after it finishes running. Let me know at this stage, I will inspect it for you. If any changes would be needed, they would most likely be using 10 principal components instead of 12.
