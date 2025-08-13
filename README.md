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

**These will (idempotently):**

- Download all required files
- Build bcf and bref3 genome references
- Generate the alt_alleles.db needed for imputation
- Run a PCA for subpopulation assignment

## Run every time you need to add a new trait to your PGS registry:

```
# An example for "insomnia"
bash scripts/run_add_pgs.sh PGS002149
```

**These will:**

- Fetch the already harmonized PGS weights file from the PGS catalogue
- Create weights files per subpopulation
- Append the PGS ID, along with other metadata, to a global registry of all PGSs added so far, stored in `/pgs/weights/harmonized`

## Run every time you need to score a user:

```
bash run_qc_genome.sh <path/to/23andMe/file>.txt
bash run_user.sh <path/to/23andMe/file>.txt
```

**These will:**

- Do quality control on the given genome file
- Run imputation on the user's genome, and generate a directory inside `/users` where all subsequent results will be saved for the user
- Imputated genome: `/user/userPath/<name/of/genome/file>_imputed_all.vcf.gz` and its corresponding indexed file `.tbi`
- Imputation quality control files: in `/user/userPath/imputation_qc_reports`
- User's plink2 files: in `/user/userPath/pfiles`
- Subpopulation assigned using PCA: `/user/userPath/ancestry.tsv`
- User's PGS scores: `/user/userPath/pgs_scores`

## Redo PCA

**Important!** PCA is used to assign each user to a subpopulation out of EUR, AMR, SAS, EAS, AFR, for PGS scoring purposes. It is **not** an ancestry, and it should not be reported as such to the user. Instead, it assigns the user into the subpopulation he/she is most similar to, specifically within the 1000 Genomes expanded panel used here. This is expected to positivelly correlate with the ancestry of the user, but it is not a robust prediction of their ancestry.

In the current pilot implementation, the PCA was run with 6 principal components for speed. The analyses will run using those PCA results, but the scores obtained with it will not be optimal. Six PCs is not enough here.

Steps that need to be taken:

- Make sure that the analyses work on your local machine with the files generated using 6 prinicipal components
- Then, delete everything inside the `pca_model` folder
- Re-run the setup script with `bash run_setup.sh`. This should re-run the PCA
- This might take some time to run
- I have already configured the `scripts/pipeline/setup.sh` script (the "setup" pipeline orchestrator) to run with 12 principal components next time by configuring the flag: `--pcs 12`

The new PCA will need to be inspected after it finishes running. Let me know at this stage, I will inspect it for you. If any changes would be needed, they would most likely be using 10 principal components instead of 12.
