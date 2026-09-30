#!/bin/bash
#$ -N A04i_fix_chr1_vcfs
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04i_fix_chr1_vcfs.$JOB_ID.$TASK_ID
#$ -j y

# Rebuild chr1 in the ambisim pool truth VCFs. The source 1000G chr1 (demux_benchmark .../by_chrom/1000G.chr1.vcf.gz)
# is corrupt (gzip crc error), so every pool VCF stops at chr1:190,673,377. This writes <tree>/vcf_fixed/<pool>.vcf.gz
# = clean chr1 (re-downloaded, A04h) + the original chr2-22/X records. The originals are untouched.
# NOTE: ambisim reads were simulated from the truncated genotypes, so the ORIGINAL vcf/ files remain the correct truth
# for the existing simulations; vcf_fixed/ is for re-simulation.
# Recipe (from the pool VCF headers): 1000G 30x phased panel -> view -v snps -m2 -M2 -i 'INFO/AF>=0.05 && INFO/AF<=.95'
#   -> view -S <pool donors, same order> (AC / AN recomputed) -> view -g ^miss.
# Task 0 (qsub -t 1 MODE=prep): filtered chr1 panel once. Tasks (MODE=pool, -t 1-N): one row of txt/pool_vcfs_chr1_fix.txt.
# Usage: j=$(qsub -terse -v MODE=prep -t 1 A04i_fix_chr1_vcfs.sh); qsub -v MODE=pool -hold_jid $j -t 1-266 A04i_fix_chr1_vcfs.sh
source ~/.bashrc
set -euo pipefail
module load bcftools htslib
PD=/u/project/cluo/terencew/claude/project_ideas/pool_design
REF=/u/project/cluo/terencew/claude/project_ideas/latent_genos/reference/1000G_30x
CHR1=$REF/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
CHR1F=$REF/1000G.chr1.common_biallelic.vcf.gz
MODE=${MODE:-pool}

# ID=1
ID=$SGE_TASK_ID
echo "$(date): $MODE task $ID on $(hostname -s)"
if [ "$MODE" = prep ]; then
    if [ ! -s "$CHR1F.tbi" ]; then
        time bcftools view --threads 4 -v snps -m2 -M2 -i 'INFO/AF>=0.05 && INFO/AF<=.95' "$CHR1" -Oz -o "$CHR1F.tmp.vcf.gz"
        mv "$CHR1F.tmp.vcf.gz" "$CHR1F"
        tabix -f -p vcf "$CHR1F"
    fi
    echo "chr1 common biallelic records: $(bcftools index -n "$CHR1F"); last: $(tabix "$CHR1F" chr1:248000000-249000000 | tail -1 | cut -f2)"
    exit 0
fi

IN=$PD/$(sed -n "${ID}p" $PD/txt/pool_vcfs_chr1_fix.txt)
TREE=$(dirname "$(dirname "$IN")")
OUTDIR=$TREE/vcf_fixed
OUT=$OUTDIR/$(basename "$IN")
mkdir -p "$OUTDIR"
if [ -s "$OUT.tbi" ]; then echo "exists: $OUT"; exit 0; fi
TMP=$(mktemp -d /u/project/cluo_scratch/terencew/claude/latent_genos/A04i.XXXXXX)
bcftools query -l "$IN" > "$TMP/samples.txt"
time bcftools view -S "$TMP/samples.txt" "$CHR1F" -Ou | bcftools view -g ^miss -Oz -o "$TMP/chr1.vcf.gz"
bcftools view -t ^chr1 "$IN" -Oz -o "$TMP/rest.vcf.gz"
# header of the rest (keeps the original INFO / contig / provenance lines); chr1 records are a subset of the same fields
bcftools concat "$TMP/chr1.vcf.gz" "$TMP/rest.vcf.gz" -Oz -o "$OUT.tmp.vcf.gz"
mv "$OUT.tmp.vcf.gz" "$OUT"
tabix -f -p vcf "$OUT"
echo "chr1 records: before $(tabix "$IN" chr1 | wc -l), after $(tabix "$OUT" chr1 | wc -l); chr1 last: $(tabix "$OUT" chr1:248000000-249000000 | tail -1 | cut -f2)"
echo "samples identical: $(cmp -s <(bcftools query -l "$IN") <(bcftools query -l "$OUT") && echo yes || echo NO)"
rm -rf "$TMP"
echo "$(date): done $OUT"
