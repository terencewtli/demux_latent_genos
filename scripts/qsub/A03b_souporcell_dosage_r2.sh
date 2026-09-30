#!/bin/bash
#$ -N A03b_souporcell_dosage_r2
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A03b_souporcell_dosage_r2.$JOB_ID
#$ -j y

# A03b on every finished souporcell row (rows with assign.tsv are skipped); writes the cluster -> donor assign.tsv
# that A04g scoring needs, plus naive dosage r2 summaries.
# Usage: qsub scripts/qsub/A03b_souporcell_dosage_r2.sh
source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools
echo "$(date): A03b on $(hostname -s), $(bcftools --version | head -1)"
time python /u/project/cluo/terencew/claude/project_ideas/latent_genos/scripts/A03b_souporcell_dosage_r2.py
echo "$(date): done"
