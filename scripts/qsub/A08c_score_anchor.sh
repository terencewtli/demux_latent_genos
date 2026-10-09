#!/bin/bash
#$ -N A08c_score_anchor
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A08c_score_anchor.$JOB_ID
#$ -j y

# A08c: eQTL-anchored accuracy + LD fidelity over all A08b chromosomes; see scripts/A08c_score_anchor.py.
# Requires A08a (leads.tsv) and A08b. Writes results/A08_eqtl_anchor/score/ (overwrites; ~1 min per chromosome).
# Submit with (from latent_genos/): qsub -hold_jid A08b_extract_dosages scripts/qsub/A08c_score_anchor.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
echo "Start: $(date)  host=$(hostname -s)"
time python "$PROJ/scripts/A08c_score_anchor.py"
echo "End: $(date)"
