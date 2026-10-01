#!/bin/bash
#$ -N A04j_leakage_score
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-3:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04j_leakage_score.$JOB_ID.$TASK_ID
#$ -j y

# TOPMed leakage test scoring with the FIXED lib/score_imputation.py (2026-09-30: bcftools 1.11 query has no -v/-m/-M,
# so the A04e scores silently covered 0 sites). Same chr20 input, three imputations:
#   task 1 topmed       TOPMed r3 server output (job-20260929-163357-077)
#   task 2 with_donors  local 1000G 30x panel WITH the 69 donors (what leakage looks like)
#   task 3 loo          local panel without donors + relatives (no leakage)
# Output: latent_genos/results/topmed_leakage/scores_v2/<arm>.chr20.{sites.tsv.gz,bins.tsv}
# Submit with: qsub scripts/latent_genos/qsub/A04j_leakage_score.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
LIB=$PROJDIR/scripts/latent_genos/lib
SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
LEAK=/u/project/cluo/terencew/claude/project_ideas/latent_genos/results/topmed_leakage
CHR=chr20
# ID=1
ID=$SGE_TASK_ID
ARMS=(topmed with_donors loo)
ARM=${ARMS[$((ID - 1))]}
if [ "$ARM" = topmed ]; then DOSE=$(ls $LEAK/topmed/job-*/output/$CHR.dose.vcf.gz); else DOSE=$LEAK/local/$ARM/$CHR.dose.vcf.gz; fi
OUT=$LEAK/scores_v2
mkdir -p $OUT
echo "Start: $(date)  host=$(hostname -s)  arm=$ARM dose=$DOSE"
[ -s "$DOSE.tbi" ] || bcftools index -t "$DOSE"
time python "$LIB/score_imputation.py" --dose "$DOSE" --typed "$LEAK/input/$CHR.vcf.gz" \
    --truth "$SRC/1000G.$CHR.vcf.gz" --samples "$LEAK/donors.txt" --region "$CHR" --out "$OUT/$ARM.$CHR"
echo "End: $(date)  arm=$ARM"
