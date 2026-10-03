#!/bin/bash
#$ -N A06d_emulate_targets
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 1
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A06d_emulate_targets.$JOB_ID
#$ -j y

# A06d: Level B GL targets for every region of <a06a_outdir> (skips regions already done).
# Usage: qsub scripts/qsub/A06d_emulate_targets.sh <a06a_outdir>
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/latent_genos
echo "$(date): A06d $1 on $(hostname -s)"
time /u/home/t/terencew/project-cluo/miniconda3/envs/scFates/bin/python $PROJDIR/scripts/A06d_emulate_targets.py $1
echo "$(date): done"
