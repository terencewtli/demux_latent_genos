#!/bin/bash
#$ -N A07i_eqtl
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07i_eqtl.$JOB_ID
#$ -j y

# A07i step 2: tensorQTL cis-eQTL per cell type with the array genotypes and every imputation arm present, on
# chr6 + chr22. Phenotypes from oracle labels (genotype arms compared on identical phenotypes), then the fully
# genotype-free version (soup labels + soup arms).
# Submit with: qsub -hold_jid A07i_pseudobulk,A07h_glimpse_bam,A07h_quilt2 scripts/qsub/A07i_eqtl.sh

source ~/.bashrc
conda activate tensorqtl
set -euo pipefail
module load bcftools
PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
echo "Start: $(date)  host=$(hostname -s)"
time python $PROJ/scripts/A07i_eqtl.py --chroms chr6,chr22 --labels oracle
time python $PROJ/scripts/A07i_eqtl.py --chroms chr6,chr22 --labels soup --arms soupsites_glimpse,bam_glimpse_soup,bam_quilt_soup
echo "End: $(date)"
