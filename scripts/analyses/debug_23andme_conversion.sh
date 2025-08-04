#!/usr/bin/env bash
###############################################################################
# Fixed section for TSV creation in impute.sh
###############################################################################

# Replace the TSV creation section in your script with this:

if [[ "$INPUT_BUILD" == 37 ]]; then
  echo "==> prepare 23andMe 4-column TSV for bcftools tsv2vcf"
  
  # Create uncompressed TSV first for debugging
  TSV_TEMP="${OUT_DIR}/${STEM}.23andme.4col.tsv"
  
  if [[ $IN_TXT == *.gz ]]; then
    zcat "$RAW_GZ"
  else
    cat "$RAW_GZ"
  fi | awk 'BEGIN{OFS="\t"} 
    # Skip comment lines
    /^#/ { next }
    # Process data lines with at least 4 fields
    NF >= 4 {
      rsid = $1
      chr = $2
      pos = $3
      genotype = $4
      
      # Clean chromosome name
      gsub(/^chr/, "", chr)
      if (chr == "23") chr = "X"
      else if (chr == "24") chr = "Y" 
      else if (chr == "M" || chr == "MT") chr = "MT"
      
      # Skip if chromosome is not standard
      if (chr !~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/) next
      
      # Skip if position is not numeric
      if (pos !~ /^[0-9]+$/) next
      
      # Skip if genotype is missing or invalid
      if (genotype == "--" || genotype == "") next
      
      print rsid, chr, pos, genotype
    }' > "$TSV_TEMP"
  
  # Check if TSV was created successfully
  if [[ ! -s "$TSV_TEMP" ]]; then
    echo "ERROR: TSV conversion failed - no data in output"
    echo "Checking input file format..."
    
    if [[ $IN_TXT == *.gz ]]; then
      echo "First 10 lines of input:"
      zcat "$RAW_GZ" | head -10
      echo "Non-comment lines count:"
      zcat "$RAW_GZ" | grep -v '^#' | wc -l
      echo "Sample data line:"
      zcat "$RAW_GZ" | grep -v '^#' | head -1 | while read line; do
        echo "Line: $line"
        echo "Fields: $(echo "$line" | awk '{print NF}')"
      done
    else
      echo "First 10 lines of input:"
      head -10 "$RAW_GZ"
      echo "Non-comment lines count:"
      grep -v '^#' "$RAW_GZ" | wc -l
      echo "Sample data line:"
      grep -v '^#' "$RAW_GZ" | head -1 | while read line; do
        echo "Line: $line"
        echo "Fields: $(echo "$line" | awk '{print NF}')"
      done
    fi
    exit 1
  fi
  
  # Compress with bgzip and index
  bgzip -c "$TSV_TEMP" > "$TSV4"
  tabix -f -s 2 -b 3 -e 3 "$TSV4"
  
  # Clean up temporary file
  rm "$TSV_TEMP"

  echo "DEBUG: RAW_GZ = $RAW_GZ"
  echo "DEBUG: TSV4 = $TSV4" 
  echo "DEBUG: TSV4 exists: $(test -f "$TSV4" && echo "YES" || echo "NO")"
  echo "DEBUG: TSV4 size: $(ls -lh "$TSV4" 2>/dev/null || echo "FILE NOT FOUND")"
  echo "DEBUG: TSV4 line count (bgzip): $(bgzip -dc "$TSV4" 2>/dev/null | wc -l || echo "CANNOT READ")"
  echo "DEBUG: TSV4 line count (zcat): $(zcat "$TSV4" 2>/dev/null | wc -l || echo "CANNOT READ")"
  echo "DEBUG: First 3 TSV4 lines:"
  bgzip -dc "$TSV4" 2>/dev/null | head -3 || echo "CANNOT READ TSV4"

  echo "==> TSV → VCF (GRCh37)"
  bcftools convert --tsv2vcf "$TSV4" -f "$FASTA_37" -s "$STEM" -Oz -o "$VCF_GZ" \
    2>&1 | tee "${OUT_DIR}/${STEM}.tsv2vcf.b37.log"
  
  # Check if VCF was created successfully
  if [[ ! -s "$VCF_GZ" ]]; then
    echo "ERROR: VCF conversion failed"
    echo "bcftools log:"
    cat "${OUT_DIR}/${STEM}.tsv2vcf.b37.log"
    exit 1
  fi
  
  tabix -f -p vcf "$VCF_GZ"

  # Rest of the b37 processing...
  echo "==> add chr-prefixes (b37)"
  bcftools annotate --rename-chrs "$ADDCHR_MAP" -Oz -o "$VCF_CHR_GZ" "$VCF_GZ"
  tabix -f -p vcf "$VCF_CHR_GZ"

  echo "==> patch missing ALT alleles"
  python scripts/helpers/alt_fix.py "$VCF_CHR_GZ"
  [[ -f "$RAW_VCF" ]] || { echo "Expected $RAW_VCF after alt_fix.py, not found"; exit 1; }

  echo "==> CrossMap liftover to GRCh38"
  crossmap vcf "$CHAIN" "$RAW_VCF" "$FASTA_38" "$LIFT_VCF"

  echo "==> sort + normalize + fixref + keep primary contigs"
  bcftools sort "$LIFT_VCF" -Oz -o "$SORT_VCF" && tabix -f -p vcf "$SORT_VCF"
  bcftools norm -f "$FASTA_38" -m -both "$SORT_VCF" -Oz -o "$NORM_VCF" && tabix -f -p vcf "$NORM_VCF"
  bcftools +fixref "$NORM_VCF" -- -f "$FASTA_38" -m flip -Oz -o "$FIXREF_VCF" && tabix -f -p vcf "$FIXREF_VCF"
  CONTIGS=$(printf 'chr%s,' {1..22} X Y M); CONTIGS=${CONTIGS%,}
  bcftools view -r "$CONTIGS" -Oz -o "$PRIM_VCF" "$FIXREF_VCF" && tabix -f -p vcf "$PRIM_VCF"
  CHR_VCF="$PRIM_VCF"

elif [[ "$INPUT_BUILD" == 38 ]]; then
  # Similar fix for build 38 - create temporary uncompressed TSV first
  echo "==> prepare 23andMe 4-column TSV for bcftools tsv2vcf"
  
  TSV_TEMP="${OUT_DIR}/${STEM}.23andme.4col.tsv"
  
  if [[ $IN_TXT == *.gz ]]; then
    zcat "$RAW_GZ"
  else
    cat "$RAW_GZ"
  fi | awk 'BEGIN{OFS="\t"} 
    /^#/ { next }
    NF >= 4 {
      rsid = $1; chr = $2; pos = $3; genotype = $4
      gsub(/^chr/, "", chr)
      if (chr == "23") chr = "X"
      else if (chr == "24") chr = "Y" 
      else if (chr == "M" || chr == "MT") chr = "MT"
      if (chr !~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/) next
      if (pos !~ /^[0-9]+$/) next
      if (genotype == "--" || genotype == "") next
      print rsid, chr, pos, genotype
    }' > "$TSV_TEMP"
  
  if [[ ! -s "$TSV_TEMP" ]]; then
    echo "ERROR: TSV conversion failed - no data in output"
    exit 1
  fi
  
  bgzip -c "$TSV_TEMP" > "$TSV4"
  tabix -f -s 2 -b 3 -e 3 "$TSV4"
  rm "$TSV_TEMP"

  # Continue with build 38 processing...
  echo "==> TSV → VCF (GRCh38)"
  VCF38_GZ="${OUT_DIR}/${STEM}.build38.vcf.gz"
  VCF38_CHR_GZ="${OUT_DIR}/${STEM}.build38.chr.vcf.gz"
  RAW38_VCF="${OUT_DIR}/${STEM}.build38.chr.alt.vcf.gz"

  bcftools convert --tsv2vcf "$TSV4" -f "$FASTA_38" -s "$STEM" -Oz -o "$VCF38_GZ" \
    2>&1 | tee "${OUT_DIR}/${STEM}.tsv2vcf.b38.log"
  
  if [[ ! -s "$VCF38_GZ" ]]; then
    echo "ERROR: VCF conversion failed"
    exit 1
  fi
  
  tabix -f -p vcf "$VCF38_GZ"

  echo "==> add chr-prefixes (b38)"
  bcftools annotate --rename-chrs "$ADDCHR_MAP" -Oz -o "$VCF38_CHR_GZ" "$VCF38_GZ"
  tabix -f -p vcf "$VCF38_CHR_GZ"

  echo "==> patch missing ALT alleles"
  python static_files/scripts/alt_fix.py "$VCF38_CHR_GZ"
  [[ -f "$RAW38_VCF" ]] || { echo "Expected $RAW38_VCF after alt_fix.py, not found"; exit 1; }

  echo "==> normalize + fixref + keep primary contigs"
  bcftools norm -f "$FASTA_38" -m -both "$RAW38_VCF" -Oz -o "$NORM_VCF" && tabix -f -p vcf "$NORM_VCF"
  bcftools +fixref "$NORM_VCF" -- -f "$FASTA_38" -m flip -Oz -o "$FIXREF_VCF" && tabix -f -p vcf "$FIXREF_VCF"
  CONTIGS=$(printf 'chr%s,' {1..22} X Y M); CONTIGS=${CONTIGS%,}
  bcftools view -r "$CONTIGS" -Oz -o "$PRIM_VCF" "$FIXREF_VCF" && tabix -f -p vcf "$PRIM_VCF"
  CHR_VCF="$PRIM_VCF"
else
  echo "✗ Unexpected INPUT_BUILD=$INPUT_BUILD"; exit 1
fi