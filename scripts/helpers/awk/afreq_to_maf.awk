# Converts plink2 .afreq into chr pos ref alt maf, inferring columns by header names.
# Robust to presence/absence of MAF and A1_FREQ.
BEGIN {
  OFS = (OFS=="\t" ? "\t" : OFS)
}
NR == 1 {
  for (i = 1; i <= NF; i++) {
    if ($i == "#CHROM" || $i == "CHROM") C = i
    else if ($i == "POS") P = i
    else if ($i == "REF") R = i
    else if ($i == "ALT") A = i
    else if ($i == "MAF") M = i
    else if ($i == "A1_FREQ") F = i
  }
  next
}
{
  chr = $C
  if (chr !~ /^chr/) chr = "chr" chr
  maf = (M ? $M : (F ? (($F > 0.5) ? 1 - $F : $F) : ""))
  print chr, $P, $R, $A, maf
}
