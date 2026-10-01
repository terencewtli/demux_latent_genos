#!/bin/bash
#$ -N A04k_score_minimac4
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-264:1
#$ -tc 30
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04k_score_minimac4.$JOB_ID.$TASK_ID
#$ -j y

# Score the Minimac4 arm (A04c) with the SAME scorer / truth as the GLIMPSE2 arm (A04g step 4,
# lib/score_latent_imputation.py): typed-site r2 vs souporcell's naive calls, untyped r2 by 1000G MAF bin.
# Task = row of txt/latent_genos_souporcell_tasks.txt x chr1-22 (same scheme as A04c / A04g); both QC sets.
# Output: latent_genos/results/A04k_minimac4_score/<qc>/<tree>__<pool>__<mod>/chrN.score.*
# Submit with: qsub -hold_jid <A04c job> scripts/latent_genos/qsub/A04k_score_minimac4.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
LIB=$PROJDIR/scripts/latent_genos/lib
PANEL=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
RES=/u/project/cluo/terencew/claude/project_ideas/latent_genos/results
LIST=$PROJDIR/txt/latent_genos_souporcell_tasks.txt
# ID=240
ID=$SGE_TASK_ID
ROW=$(( (ID - 1) / 22 + 1 ))
C=$(( (ID - 1) % 22 + 1 ))
CHR=chr$C
read -r TREE POOL K MOD <<< "$(sed -n "${ROW}p" "$LIST")"
TAG=${TREE}__${POOL}__${MOD}
IMP=$PROJDIR/$TREE/$POOL/impute/$MOD
SOUP=$PROJDIR/$TREE/$POOL/demux/souporcell/$MOD/$SAMPLE/cluster_genotypes.vcf
ASSIGN=$RES/A03b_souporcell_dosage_r2/$TAG/assign.tsv
echo "Start: $(date)  host=$(hostname -s)  tag=$TAG chr=$CHR"
[ -s "$ASSIGN" ] || { echo "no A03b assign.tsv for $TAG"; exit 1; }
for QC in server all; do
    DOSE=$IMP/$QC/$CHR.dose.vcf.gz
    OUT=$RES/A04k_minimac4_score/$QC/$TAG
    mkdir -p "$OUT"
    if [ ! -s "$DOSE.tbi" ]; then echo "$(date): [$QC] no dose for $CHR, skipping"; continue; fi
    if [ -s "$OUT/$CHR.score.donor.tsv" ]; then echo "$(date): [$QC] scores exist, skipping"; continue; fi
    time python "$LIB/score_latent_imputation.py" --dose "$DOSE" --typed "$IMP/target.$QC.vcf.gz" \
        --soup "$SOUP" --assign "$ASSIGN" --truth "$PROJDIR/$TREE/vcf/$POOL.vcf.gz" \
        --panel-sites "$PANEL/1000G_30x.$CHR.sites.vcf.gz" --region "$CHR" --out "$OUT/$CHR.score"
done
echo "End: $(date)  tag=$TAG chr=$CHR"
