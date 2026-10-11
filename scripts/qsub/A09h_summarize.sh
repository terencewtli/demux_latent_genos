#!/bin/bash
#$ -N A09h_summarize
#$ -cwd
#$ -l h_data=8G,h_rt=1:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09h_summarize.$JOB_ID
#$ -j y

# A09h summary over all chromosomes and modalities; see scripts/A09h_summarize.py. Writes results/A09_multiome/consistency/*_summary.tsv.

source ~/.bashrc
conda activate allcools
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
echo "Start: $(date)"
time python "$PROJ/scripts/A09h_summarize.py" --mods multi atac gex
echo "End: $(date)"
