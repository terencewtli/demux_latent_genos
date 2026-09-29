#!/bin/bash
#$ -N A04e_leakage_local
#$ -cwd
#$ -l h_data=2G,h_rt=6:00:00
#$ -pe shared 10
#$ -t 1-2:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04e_leakage_local.$JOB_ID.$TASK_ID
#$ -j y

# Local calibration for the TOPMed leakage test. Impute the same thinned chr20
# truth input (A04d) with the local 1000G 30x panel, two ways:
#   task 1  with_donors: the 69 donors ARE in the panel (what leakage looks like)
#   task 2  loo:         donors + 1000G relatives left out (what no leakage looks like)
# Then score both, the same way the TOPMed server output will be scored
# (lib/score_imputation.py): r2 by 1000G AF bin, plus variants private to the
# donors within 1000G.
#
# Requires A04a chr20 msav and A04d input.
# Submit with: qsub -hold_jid A04a_build_panel_msav scripts/latent_genos/qsub/A04e_leakage_local.sh

source ~/.bashrc

conda activate allcools

set -euo pipefail

module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
LIB=$PROJDIR/scripts/latent_genos/lib
PANEL=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
EAGLE=/u/project/cluo/terencew/programs/Eagle_v2.4.1/eagle
MAP=/u/project/cluo/terencew/programs/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
MINIMAC4=/u/project/cluo/terencew/programs/build/Minimac4/build/minimac4
LEAK=/u/project/cluo/terencew/claude/project_ideas/latent_genos/results/topmed_leakage
TMPROOT=/u/project/cluo_scratch/terencew/claude/latent_genos/A04e
CHR=chr20
N=${NSLOTS:-1}

# ID=2
ID=$SGE_TASK_ID
ARMS=(with_donors loo)
ARM=${ARMS[$((ID - 1))]}

TGT=$LEAK/input/$CHR.vcf.gz
MSAV=$PANEL/1000G_30x.$CHR.msav
OUT=$LEAK/local/$ARM
TMP=$TMPROOT/$ARM

echo "Start: $(date)  host=$(hostname -s)  arm=$ARM slots=$N"
[ -s "$TGT" ] || { echo "ERROR: $TGT missing (run A04d)" >&2; exit 1; }
[ -s "$MSAV" ] || { echo "ERROR: $MSAV missing (run A04a)" >&2; exit 1; }

mkdir -p "$OUT/loo" "$TMP"
DOSE=$OUT/$CHR.dose.vcf.gz

if [ -s "$DOSE.tbi" ]; then
    echo "$(date): $DOSE exists, skipping imputation"
else
    bcftools query -l "$SRC/1000G.$CHR.vcf.gz" > "$OUT/loo/panel_samples.txt"
    if [ "$ARM" = "loo" ]; then
        python "$LIB/leave_pool_out.py" --donors "$LEAK/donors.txt" \
            --panel-samples "$OUT/loo/panel_samples.txt" \
            --ped "$PANEL/20130606_g1k_3202_samples_ped_population.txt" --outdir "$OUT/loo"
    else
        : > "$OUT/loo/exclude.txt"
        cp "$OUT/loo/panel_samples.txt" "$OUT/loo/keep.txt"
    fi
    echo "panel samples kept: $(wc -l < "$OUT/loo/keep.txt")"

    REF=$TMP/ref.bcf
    time bcftools view --threads "$N" -S "$OUT/loo/keep.txt" -T "$TGT" -m2 -M2 -v snps \
        "$SRC/1000G.$CHR.vcf.gz" -Ob -o "$REF"
    bcftools index -f "$REF"

    time $EAGLE --vcfRef "$REF" --vcfTarget "$TGT" --geneticMapFile "$MAP" \
        --outPrefix "$TMP/phased" --vcfOutFormat z --numThreads "$N"
    bcftools index -t -f "$TMP/phased.vcf.gz"

    time $MINIMAC4 "$MSAV" "$TMP/phased.vcf.gz" \
        --sample-ids-file "$OUT/loo/keep.txt" \
        -f GT,DS,HDS,GP,SD -a \
        -e "$OUT/$CHR.empirical.vcf.gz" \
        -o "$TMP/dose.vcf.gz" -O vcf.gz \
        -t "$N" --temp-prefix "$TMP/m4_"
    mv "$TMP/phased.vcf.gz" "$OUT/$CHR.phased.vcf.gz"
    mv "$TMP/dose.vcf.gz" "$DOSE"
    bcftools index -t -f "$DOSE"
    rm -f "$REF" "$REF.csi" "$TMP/phased.vcf.gz.tbi"
fi

time python "$LIB/score_imputation.py" --dose "$DOSE" --typed "$TGT" \
    --truth "$SRC/1000G.$CHR.vcf.gz" --samples "$LEAK/donors.txt" \
    --region "$CHR" --out "$OUT/$CHR.score"

echo "End: $(date)  arm=$ARM"
