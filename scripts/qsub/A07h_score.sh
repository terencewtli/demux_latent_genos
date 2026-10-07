#!/bin/bash
#$ -N A07h_score
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07h_score.$JOB_ID
#$ -j y

# A07h scoring of every imputation arm present (A07f, GLIMPSE2-on-BAM soup / oracle, QUILT2 soup / oracle) against
# the OneK1K array on chr6 + chr22; arms without output are skipped, so it can be rerun as arms finish.
# Submit with: qsub -hold_jid A07h_glimpse_bam,A07h_quilt2 scripts/qsub/A07h_score.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools
PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
echo "Start: $(date)  host=$(hostname -s)"
time python $PROJ/scripts/A07h_score.py --chroms chr6,chr22
echo "End: $(date)"
