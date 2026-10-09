#!/bin/bash
#$ -N A08d_prep_loci
#$ -cwd
#$ -l h_data=8G,h_rt=12:00:00
#$ -pe shared 4
#$ -t 1-5:1
#$ -tc 5
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A08d_prep_loci.$JOB_ID.$TASK_ID
#$ -j y

# A08d prep, one cell type per task: z-scores + four LD-source genotype matrices for up to MAX_LOCI random loci
# (A08a loci.tsv); see scripts/A08d_prep_loci.py. ~3 min table read + ~1 min per locus (OneK1K all-sample query).
# Requires A08a, A08b (all chromosomes). Submit with (from latent_genos/):
#   qsub -hold_jid A08b_extract_dosages [-v MAX_LOCI=200] scripts/qsub/A08d_prep_loci.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
CTS=(cd4nc cd8et nk bin monoc)
# ID=1
ID=$SGE_TASK_ID
CT=${CTS[$((ID - 1))]}
echo "Start: $(date)  host=$(hostname -s)  ct=$CT"
time python "$PROJ/scripts/A08d_prep_loci.py" --ct "$CT" --max-loci "${MAX_LOCI:-200}"
echo "End: $(date)  ct=$CT"
