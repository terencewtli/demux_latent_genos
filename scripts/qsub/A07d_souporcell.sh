#!/bin/bash
#$ -N A07d_souporcell
#$ -cwd
#$ -l h_data=8G,h_rt=16:00:00
#$ -pe shared 4
#$ -t 1-5:1
#$ -tc 5
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07d_souporcell.$JOB_ID.$TASK_ID
#$ -j y

# OneK1K pilot: genotype-free souporcell on one pool's STARsolo BAM, k = number of pooled donors (GEO).
# Same settings as the simulation GEX arm (pool_design/scripts/latent_genos/qsub/A03a_souporcell.sh):
# --skip_remap with the 1000G AF > 5% common-variant list, 200 restarts. Fasta = the genome the
# STAR index was built from (hg38_igvf no-alt analysis set).
# Output: $SCR/A07/souporcell/pool<P>/{clusters.tsv,cluster_genotypes.vcf,ambient_rna.txt,...}
# Requires A07b. Submit with (from latent_genos/): qsub -hold_jid A07b_starsolo scripts/qsub/A07d_souporcell.sh

source ~/.bashrc
conda activate souporcell
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos/A07
SOUPORCELL=/u/project/cluo/terencew/programs/souporcell/souporcell_pipeline.py
FASTA=/u/project/cluo/terencew/reference/hg38_igvf/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
COMMON=/u/project/cluo/terencew/reference/cellsnp_db/hg38.AF5e2.possorted.chr_prefix.vcf
# ID=1
ID=$SGE_TASK_ID
read -r POOL GSM K <<< "$(sed -n "${ID}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
BAM=$SCR/starsolo/pool$POOL/Aligned.sortedByCoord.out.bam
BC=$SCR/starsolo/pool$POOL/Solo.out/Gene/filtered/barcodes.tsv
OUT=$SCR/souporcell/pool$POOL

echo "Start: $(date)  host=$(hostname -s)  pool=$POOL k=$K slots=${NSLOTS:-1}"
[ -s "$BAM.bai" ] && [ -s "$BC" ] || { echo "ERROR: STARsolo output missing for pool $POOL" >&2; exit 1; }
if [ -s "$OUT/clusters.tsv" ] && [ -s "$OUT/cluster_genotypes.vcf" ]; then echo "$(date): done, skipping"; exit 0; fi
mkdir -p "$OUT"
echo "cells: $(wc -l < "$BC")"

time $SOUPORCELL -i "$BAM" -b "$BC" -f "$FASTA" -t "${NSLOTS:-1}" -o "$OUT" -k "$K" \
    --common_variants "$COMMON" --restarts 200 --skip_remap True

awk 'NR>1 {n[$2]++} END {for (s in n) print s, n[s]}' "$OUT/clusters.tsv"
echo "End: $(date)  pool=$POOL"
