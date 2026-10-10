#!/bin/bash
#$ -N A09g_score
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-31:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09g_score.$JOB_ID.$TASK_ID
#$ -j y

# A09g: score one A09e target (row of txt/A09_targets.txt) against WGS (FT=PASS, GQ>=20): aggregate r2 / concordance
# by typed vs untyped and panel MAF bin, plus the no-imputation (PL argmax) baseline at typed sites; see scripts/A09g_score.py.
# Writes results/A09_multiome/score/<target>/{sums,metrics}.tsv (overwrites). Requires A09e for that target.
# Submit with (from latent_genos/): qsub -hold_jid A09e_glimpse_impute scripts/qsub/A09g_score.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
# ID=1
ID=$SGE_TASK_ID
TARGET=$(sed -n "${ID}p" "$PROJ/txt/A09_targets.txt")
echo "Start: $(date)  host=$(hostname -s)  target=$TARGET"
time python "$PROJ/scripts/A09g_score.py" --target "$TARGET"
echo "End: $(date)  target=$TARGET"
