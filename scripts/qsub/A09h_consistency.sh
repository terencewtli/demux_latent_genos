#!/bin/bash
#$ -N A09h_consistency
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-66:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09h_consistency.$JOB_ID.$TASK_ID
#$ -j y

# A09h: cross-pool (= time point) consistency of each donor's imputed genotypes, one modality x chromosome per task
# (row of txt/A09h_tasks.txt: mod<TAB>chrom); see scripts/A09h_consistency.py.
# Writes results/A09_multiome/consistency/<mod>/<chrom>.{pairs,kwrong,consensus}.tsv (overwrites). Needs A09e.
# Then: qsub -hold_jid A09h_consistency scripts/qsub/A09h_summarize.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
# ID=1
ID=$SGE_TASK_ID
read MOD CHROM < <(sed -n "${ID}p" "$PROJ/txt/A09h_tasks.txt")
echo "Start: $(date)  host=$(hostname -s)  mod=$MOD chrom=$CHROM"
time python "$PROJ/scripts/A09h_consistency.py" --mod "$MOD" --chrom "$CHROM"
echo "End: $(date)"
