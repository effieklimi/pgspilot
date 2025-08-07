docker build -t pgspilot .

docker run -it --rm \
 -v "$(pwd)/genome_data:/app/genome_data:ro" \
 -v "$(pwd)/input_data:/app/input_data:ro" \
 -v "$(pwd)/pca_model:/app/pca_model" \
 -v "$(pwd)/users:/app/users" \
 -v "$(pwd)/traits:/app/traits" \
 -v "$(pwd)/weights_hm:/app/weights_hm" \
 -v "$(pwd)/weights_raw:/app/weights_raw" \
 -v "$(pwd)/calibration:/app/calibration" \
 pgspilot





 # End-to-end checklist

1. **Run once (per reference update):**

   - Build PCA model:

   ```
   $PYTHON fit_pca_1kg.py \
     --vcf-pattern "$ONEKG_VCF_PATTERN" \
     --labels "$ONEKG_LABELS" \
     --out "$PCA_DIR" \
     --pcs 6 --maf 0.05 --max-missing 0.05 --thin-kb 50 --max-snps 50000
   ```

   - Build/export scorable_sites.b38.tsv (your “production-lite #3”).

2. **Every time you want to add/update PGS files/traits:**

   ```
   ./prep_pgs.sh weights_manifest.csv
   ```

3. **Per each user:**

```
./user.sh users/USER123/USER123_23andme.txt.gz
# or if you already have the imputed VCF:
./user.sh users/USER123/USER123_imputed_all.vcf.gz --traits "Height,HeartRate"
```
