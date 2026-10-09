#!/bin/bash
#$ -N A07i_raw_arms
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-5:1
#$ -tc 5
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07i_raw_arms.$JOB_ID.$TASK_ID
#$ -j y

# No-imputation baseline arms (raw_gt, raw_gl) for one pool (row of txt/onek1k_pilot_pools.txt); see
# scripts/A07i_raw_arms.py. Requires A07d (souporcell) and A07f (target.bcf) for the pool.
# SFX=_v2 (qsub -v SFX=_v2): the A07d rerun pools.
# Submit with (from latent_genos/): qsub -t 2-3 scripts/qsub/A07i_raw_arms.sh            (pools 1, 55; v1)
#                                   qsub -v SFX=_v2 -t 1 (then -t 4-5) -hold_jid A07f_glimpse_impute scripts/qsub/A07i_raw_arms.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
# ID=2
ID=$SGE_TASK_ID
read -r POOL GSM K <<< "$(sed -n "${ID}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
echo "Start: $(date)  host=$(hostname -s)  pool=$POOL${SFX:-}"
time python "$PROJ/scripts/A07i_raw_arms.py" --pool "$POOL" --suffix "${SFX:-}"
echo "End: $(date)  pool=$POOL"
