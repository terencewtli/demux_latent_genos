#!/bin/bash
#$ -N A07h_split_bams
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-10:1
#$ -tc 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07h_split_bams.$JOB_ID.$TASK_ID
#$ -j y

# A07h step 1: per-group BAMs (chr6, chr22) for one pool x labelling (soup = souporcell singlets, oracle = GEO
# donor labels). ~2.5 min per chromosome per pool (chr22 test). Task -> row = (ID-1)/2 + 1, mode = soup|oracle.
# Requires A07d (souporcell) and A07g (assign.tsv: oracle naming convention). Submit with:
#   qsub -hold_jid A07g_onek1k_score scripts/qsub/A07h_split_bams.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools
PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
# ID=1
ID=$SGE_TASK_ID
ROW=$(( (ID - 1) / 2 + 1 )); MODES=(soup oracle); MODE=${MODES[$(( (ID - 1) % 2 ))]}
read -r POOL GSM K <<< "$(sed -n "${ROW}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
echo "Start: $(date)  host=$(hostname -s)  pool=$POOL mode=$MODE"
time python $PROJ/scripts/A07h_split_bams.py --pool "$POOL" --mode "$MODE" --chroms chr6,chr22
echo "End: $(date)"
