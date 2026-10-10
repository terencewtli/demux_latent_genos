#!/bin/bash
#$ -N A09c_soup_consensus
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-20
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09c_soup_consensus.$JOB_ID.$TASK_ID
#$ -j y
# souporcell consensus.py (cluster genotypes + ambient RNA) per A09 run, via `souporcell_gpu pipeline --steps consensus`.
# Submit: qsub -hold_jid A09b_soup_cluster_gpu scripts/qsub/A09c_soup_consensus.sh
set -euo pipefail
PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos/A09/souporcell_v2
FASTA=/u/project/cluo/terencew/reference/hg38_igvf/IGVF_hg38/fasta/genome.fa
COMMON=/u/project/cluo/terencew/reference/cellsnp_db/hg38.AF5e2.possorted.chr_prefix.vcf
PY=/u/home/t/terencew/project-cluo/miniconda3/envs/soupgpu/bin/python
export PYTHONPATH=/u/project/cluo/terencew/claude/exploration/engineer/code/souporcell_gpu/src
# ID=1
ID=$SGE_TASK_ID
read -r POOL MOD <<< "$(sed -n "${ID}p" "$PROJ/txt/A09_runs.txt")"
O=$SCR/${POOL}_${MOD}
echo "Start: $(date)  host=$(hostname -s)  pool=$POOL mod=$MOD"
time $PY -m souporcell_gpu pipeline -i NA -b "$O/barcodes.tsv" -f "$FASTA" -k 4 -o "$O" --common_variants "$COMMON" \
    --steps consensus --device cpu
cat "$O/ambient_rna.txt"
echo "End: $(date)"
