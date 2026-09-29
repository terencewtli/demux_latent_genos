#!/bin/bash
# Leakage test input for the TOPMed Imputation Server: are the 1000G pool donors
# in the TOPMed r3 panel?
#
# Input = the donors' TRUE 1000G 30x genotypes on chr20, unphased, thinned to
# array-like density: every 10th biallelic SNV with AF > 5% (~20k sites). After
# imputation, look at variants that are rare in 1000G and carried by a donor.
# If the donors are in the panel, those come back almost perfectly (r2 ~ 1);
# if not, rare-variant r2 should look like normal out-of-panel imputation.
# Calibrate by imputing the same input locally with the 1000G panel with and
# without the donors (A04c machinery).
#
# Donors = all unique donors in txt/latent_genos_souporcell_tasks.txt pools.
# Output: latent_genos/results/topmed_leakage/input/chr20.vcf.gz
# Run directly (small): bash scripts/latent_genos/qsub/A04d_topmed_leakage_prep.sh

set -euo pipefail
module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
OUT=/u/project/cluo/terencew/claude/project_ideas/latent_genos/results/topmed_leakage
LIST=$PROJDIR/txt/latent_genos_souporcell_tasks.txt
CHR=chr20
THIN=10

mkdir -p "$OUT/input"
if [ -s "$OUT/input/$CHR.vcf.gz" ]; then
    echo "$(date): input exists, skipping"
    exit 0
fi

awk '{print $1"/vcf/"$2".vcf.gz"}' "$LIST" | sort -u | while read -r V; do
    bcftools query -l "$PROJDIR/$V"
done | sort -u > "$OUT/donors.txt"
echo "$(date): $(wc -l < "$OUT/donors.txt") unique donors"

time bcftools view -S "$OUT/donors.txt" -m2 -M2 -v snps -i 'INFO/AF>0.05 && INFO/AF<0.95' \
    "$SRC/1000G.$CHR.vcf.gz" -Ou \
  | bcftools annotate -x INFO,^FORMAT/GT -Ou \
  | bcftools +setGT -Ou -- -t a -n u 2>/dev/null \
  | bcftools view -H -Ov \
  | awk -v t=$THIN 'NR % t == 1' > "$OUT/input/body.tmp"

bcftools view -h -S "$OUT/donors.txt" "$SRC/1000G.$CHR.vcf.gz" \
  | grep -vE '^##(INFO|bcftools|source)' > "$OUT/input/header.tmp"
cat "$OUT/input/header.tmp" "$OUT/input/body.tmp" | bcftools view -Oz -o "$OUT/input/$CHR.vcf.gz"
bcftools index -t "$OUT/input/$CHR.vcf.gz"
rm "$OUT/input/"*.tmp
echo "$(date): $(bcftools index -n "$OUT/input/$CHR.vcf.gz") sites x $(bcftools query -l "$OUT/input/$CHR.vcf.gz" | wc -l) samples"
