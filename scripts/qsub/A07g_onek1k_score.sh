#!/bin/bash
#$ -N A07g_onek1k_score
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-5:1
#$ -tc 5
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07g_onek1k_score.$JOB_ID.$TASK_ID
#$ -j y

# OneK1K pilot scoring for one pool (row of txt/onek1k_pilot_pools.txt); see scripts/A07g_onek1k_score.py.
# Requires A07c (truth), A07d (souporcell), A07f (imputation, all chromosomes of the pool).
# Submit with (from latent_genos/): qsub -hold_jid A07f_glimpse_impute scripts/qsub/A07g_onek1k_score.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
# ID=1
ID=$SGE_TASK_ID
read -r POOL GSM K <<< "$(sed -n "${ID}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
echo "Start: $(date)  host=$(hostname -s)  pool=$POOL"
time python "$PROJ/scripts/A07g_onek1k_score.py" --pool "$POOL"
echo "End: $(date)  pool=$POOL"
