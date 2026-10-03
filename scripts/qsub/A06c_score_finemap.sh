#!/bin/bash
#$ -N A06c_score_finemap
#$ -cwd
#$ -l h_data=8G,h_rt=1:00:00
#$ -pe shared 1
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A06c_score_finemap.$JOB_ID
#$ -j y

# A06c: score A06b output (truth vs arms) and apply the go/no-go rule.
# Usage: qsub [-hold_jid <A06b job>] scripts/qsub/A06c_score_finemap.sh <a06b_outdir> <results_dir>
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/latent_genos
echo "$(date): A06c $1 -> $2 on $(hostname -s)"
time /u/home/t/terencew/project-cluo/miniconda3/envs/scFates/bin/python $PROJDIR/scripts/A06c_score_finemap.py $1 $2
echo "$(date): done"
