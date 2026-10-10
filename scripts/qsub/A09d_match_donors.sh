#!/bin/bash
#$ -N A09d_match_donors
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09d_match_donors.$JOB_ID
#$ -j y

# A09d: souporcell cluster->donor map (WGS concordance) + agreement with ambimux and GEX vs ATAC; see scripts/A09d_match_donors.py.
# Requires A09c. Writes results/A09_multiome/match/.
# Submit with (from latent_genos/): qsub scripts/qsub/A09d_match_donors.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
echo "Start: $(date)  host=$(hostname -s)"
time python "$PROJ/scripts/A09d_match_donors.py"
echo "End: $(date)"
