# Pilot pipelines for imputation and PGS scoring

### Run once at the start:

1. Generate a docker image:

```
cd pgspilot
docker build --no-cache -t pgspilot .
```

2. Perform initial setup:

```
bash scripts/run_setup.sh
```

**These will (idempotently):**

- Download all required files
- Build bcf and bref3 genome references
- Generate the alt_alleles.db needed for imputation
- Run a PCA for subpopulation assignment

### Run every time you need to add a new trait to your PGS registry:

```
# An example for "insomnia"
bash scripts/run_add_pgs.sh PGS002149
```

**These will:**

- Fetch the already harmonized PGS weights file from the PGS catalogue
- Create weights files per subpopulation
- Append the PGS ID, along with other metadata, to a global registry of all PGSs added so far, stored in `/pgs/weights/harmonized`

### Run every time you need to score a user:

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
