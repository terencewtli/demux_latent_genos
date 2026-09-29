#!/bin/bash
#$ -N A04a_build_panel_msav
#$ -cwd
#$ -l h_data=2G,h_rt=8:00:00
#$ -pe shared 10
#$ -t 1-22:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04a_build_panel_msav.$JOB_ID.$TASK_ID
#$ -j y

# Build the local Minimac4 v4 reference for one chromosome from the NYGC 1000G
# 30x phased panel (3,202 samples). This is the FULL panel. Leave-pool-out
# happens at run time in A04c: minimac4 --sample-ids-file keeps only the listed
# samples, and Eagle gets a per-pool subset. So one msav per chromosome serves
# every pool.
#
# Filters: biallelic SNVs + indels only. The panel's symbolic SVs are dropped.
# Outputs (in $PANEL):
#   1000G_30x.chrN.msav          minimac4 reference
#   1000G_30x.chrN.sites.vcf.gz  sites-only, used by A04b to allele-match targets
#
# Submit with: qsub scripts/latent_genos/qsub/A04a_build_panel_msav.sh

source ~/.bashrc

set -euo pipefail

module load bcftools

SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
PANEL=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
MINIMAC4=/u/project/cluo/terencew/programs/build/Minimac4/build/minimac4
TMP=/u/project/cluo_scratch/terencew/claude/latent_genos/A04a

# ID=22
ID=$SGE_TASK_ID
CHR=chr$ID
IN=$SRC/1000G.$CHR.vcf.gz
MSAV=$PANEL/1000G_30x.$CHR.msav
SITES=$PANEL/1000G_30x.$CHR.sites.vcf.gz

echo "Start: $(date)  host=$(hostname -s)  chr=$CHR slots=${NSLOTS:-1}"
mkdir -p "$PANEL" "$TMP"

if [ -s "$MSAV" ] && [ -s "$SITES.tbi" ]; then
    echo "$(date): $MSAV exists, skipping"
    exit 0
fi

BCF=$TMP/1000G_30x.$CHR.biallelic.bcf
time bcftools view --threads "${NSLOTS:-1}" -m2 -M2 -v snps,indels "$IN" -Ob -o "$BCF"
bcftools index -f "$BCF"
echo "$(date): kept $(bcftools index -n "$BCF") biallelic SNV/indel records; samples $(bcftools query -l "$BCF" | wc -l)"

time bcftools view -G "$BCF" -Oz -o "$SITES.tmp.vcf.gz"
mv "$SITES.tmp.vcf.gz" "$SITES"
bcftools index -t -f "$SITES"

time $MINIMAC4 --compress-reference "$BCF" -t "${NSLOTS:-1}" -o "$MSAV.tmp"
mv "$MSAV.tmp" "$MSAV"

rm -f "$BCF" "$BCF.csi"
ls -la "$MSAV" "$SITES"
echo "End: $(date)  chr=$CHR"
