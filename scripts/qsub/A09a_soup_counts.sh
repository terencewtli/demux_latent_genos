#!/bin/bash
#$ -N A09a_soup_counts
#$ -cwd
#$ -l h_data=8G,h_rt=6:00:00
#$ -pe shared 4
#$ -t 1-20
#$ -tc 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09a_soup_counts.$JOB_ID.$TASK_ID
#$ -j y
# IGVF fibroblast->iPSC multiome (10 pools x 4 PGP donors): souporcell common-variants mode, CPU stages
# (covered common variants, allele counts with sc_bam) via `souporcell_gpu pipeline --steps variants,counts`.
# Settings as A03a / A07d: 1000G AF > 5% common variants, min_ref / min_alt 10, --no_umi for ATAC.
# Rows: txt/A09_runs.txt (<pool> <gex|atac>). Then A09b (GPU: cluster + troublet), A09c (consensus).
# Output: $SCR/<pool>_<mod>/ (souporcell layout and .done files). souporcell_v2: the first run (souporcell/) wrote
# empty outputs for 18 runs (stale BAM indexes); that directory has since been emptied.
source ~/.bashrc
conda activate souporcell          # samtools 1.9 (as souporcell), bedtools / vartrix in ~/bin
set -euo pipefail
PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos/A09/souporcell_v2
CR=/u/project/cluo/terencew/igvf/2023_YR2/multiome/mapping/cr_arc/igvf_ref/default
FASTA=/u/project/cluo/terencew/reference/hg38_igvf/IGVF_hg38/fasta/genome.fa
COMMON=/u/project/cluo/terencew/reference/cellsnp_db/hg38.AF5e2.possorted.chr_prefix.vcf
PY=/u/home/t/terencew/project-cluo/miniconda3/envs/soupgpu/bin/python
export PYTHONPATH=/u/project/cluo/terencew/claude/exploration/engineer/code/souporcell_gpu/src
# ID=1
ID=$SGE_TASK_ID
read -r POOL MOD <<< "$(sed -n "${ID}p" "$PROJ/txt/A09_runs.txt")"
SRC=$CR/$POOL/outs/${MOD}_possorted_bam.bam
# 9 of the 10 pools' BAMs were rewritten after their .bai (index older than the BAM -> samtools region reads fail),
# so read each BAM through a scratch symlink with a fresh index next to it; the originals are not touched.
BAMDIR=/u/project/cluo_scratch/terencew/claude/latent_genos/A09/bam
mkdir -p "$BAMDIR"
BAM=$BAMDIR/${POOL}_${MOD}.bam
ln -sfn "$SRC" "$BAM"
if [ ! -s "$BAM.bai" ]; then
    if [ "$SRC" -nt "$SRC.bai" ]; then
        time samtools index -@ "${NSLOTS:-1}" "$SRC" "$BAM.bai.tmp" && mv "$BAM.bai.tmp" "$BAM.bai"
    else
        ln -sfn "$SRC.bai" "$BAM.bai"
    fi
fi
echo "reads in chr22:20.0-20.1 Mb: $(samtools view -c "$BAM" chr22:20000000-20100000)"      # fails here if the index is still bad
BC=$CR/$POOL/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
EXTRA=()
[ "$MOD" = "atac" ] && EXTRA=(--no_umi)
echo "Start: $(date)  host=$(hostname -s)  pool=$POOL mod=$MOD slots=${NSLOTS:-1}"
time $PY -m souporcell_gpu pipeline -i "$BAM" -b "$BC" -f "$FASTA" -k 4 -o "$SCR/${POOL}_${MOD}" \
    --common_variants "$COMMON" -t "${NSLOTS:-1}" --restarts 200 --steps variants,counts \
    --samtools "$(which samtools)" --bedtools "$(which bedtools)" ${EXTRA[@]+"${EXTRA[@]}"}
echo "End: $(date)"
