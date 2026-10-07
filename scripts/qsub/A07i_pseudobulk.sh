#!/bin/bash
#$ -N A07i_pseudobulk
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07i_pseudobulk.$JOB_ID
#$ -j y

# A07i step 1: donor x cell-type pseudobulk, oracle (GEO) and soup (souporcell) cell labels.
# Submit with: qsub -hold_jid A07g_onek1k_score scripts/qsub/A07i_pseudobulk.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
echo "Start: $(date)  host=$(hostname -s)"
for L in oracle soup; do
    time python $PROJ/scripts/A07i_pseudobulk.py --labels $L
done
echo "End: $(date)"
