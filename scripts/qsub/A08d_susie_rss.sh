#!/bin/bash
#$ -N A08d_susie_rss
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-5:1
#$ -tc 5
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A08d_susie_rss.$JOB_ID.$TASK_ID
#$ -j y

# A08d: SuSiE-RSS per locus x four LD arms, one cell type per task; see scripts/A08d_susie_rss.R (~30 s per locus).
# Skips loci with susie_summary.tsv. After all tasks: python scripts/A08d_summarize.py (seconds).
# Submit with (from latent_genos/): qsub -hold_jid A08d_prep_loci [-v LAMBDA=0.1] scripts/qsub/A08d_susie_rss.sh

source /u/local/Modules/default/init/modules.sh
module load intel/2020.4 gcc/10.2.0 R/4.1.0
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
CTS=(cd4nc cd8et nk bin monoc)
# ID=1
ID=$SGE_TASK_ID
CT=${CTS[$((ID - 1))]}
echo "Start: $(date)  host=$(hostname -s)  ct=$CT"
time Rscript "$PROJ/scripts/A08d_susie_rss.R" --ct "$CT" --lambda "${LAMBDA:-0.1}"
echo "End: $(date)  ct=$CT"
