#!/bin/bash
#$ -N A08b_extract_dosages
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-22:1
#$ -tc 22
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A08b_extract_dosages.$JOB_ID.$TASK_ID
#$ -j y

# A08b: per chromosome, latent (A07f GLIMPSE2) and OneK1K dosages on the same pilot donors + souporcell site context;
# see scripts/A08b_extract_dosages.py. Pools/suffixes are fixed in POOLS there (70/11/19 = _v2); run after every
# A07f pool is imputed. Skips a chromosome whose output exists (remove it to rerun).
# Submit with (from latent_genos/): qsub -hold_jid <A07f v2 job ids> scripts/qsub/A08b_extract_dosages.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
# ID=22
ID=$SGE_TASK_ID
CHR=chr$ID
echo "Start: $(date)  host=$(hostname -s)  $CHR"
if [ -s "$PROJ/results/A08_eqtl_anchor/dosage/$CHR.donors.txt" ]; then echo "$CHR exists, skipping"; exit 0; fi
time python "$PROJ/scripts/A08b_extract_dosages.py" --chrom "$CHR"
echo "End: $(date)  $CHR"
