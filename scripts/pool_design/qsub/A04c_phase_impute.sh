#!/bin/bash
#$ -N A04c_phase_impute
#$ -cwd
#$ -l h_data=2G,h_rt=4:00:00
#$ -pe shared 10
#$ -t 1-264:1
#$ -tc 50
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04c_phase_impute.$JOB_ID.$TASK_ID
#$ -j y

# Eagle2 phasing + Minimac4 imputation of one (souporcell run x chromosome),
# using the local 1000G 30x panel with the pool's donors and relatives left out.
# This mirrors the TOPMed server pipeline (Eagle v2.4 reference-based
# pre-phasing, then Minimac4 with 20 Mb chunks). The panel differs: 1000G 30x
# leave-pool-out instead of TOPMed r3.
#
# Task ID -> row of txt/latent_genos_souporcell_tasks.txt (12 rows) x chr1-22:
#   row = (ID-1)/22 + 1, chr = (ID-1)%22 + 1
# Each task runs both QC sets from A04b (server, all).
#
# Leave-pool-out:
#   Eagle:    reference = panel subset to loo/keep.txt at the target's sites
#   Minimac4: full-panel msav + --sample-ids-file loo/keep.txt
#
# Output: <tree>/<pool>/impute/<mod>/<qc>/chrN.{phased,dose,empirical}.vcf.gz
#   dose: GT,DS,HDS,GP,SD for every panel site, plus typed-only sites (-a);
#         INFO R2 = minimac Rsq
#   empirical: leave-one-out dosages at typed sites (the server's empiricalDose)
#
# Requires A04a and A04b to be finished.
# Submit with: qsub scripts/latent_genos/qsub/A04c_phase_impute.sh

source ~/.bashrc

set -euo pipefail

module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
PANEL=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
# chr1: the by_chrom copy is corrupt (gzip crc; VCF parse error in Eagle, 2026-10-01), use the clean EBI re-download (A04h)
CHR1=/u/project/cluo/terencew/claude/project_ideas/latent_genos/reference/1000G_30x/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
EAGLE=/u/project/cluo/terencew/programs/Eagle_v2.4.1/eagle
MAP=/u/project/cluo/terencew/programs/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
MINIMAC4=/u/project/cluo/terencew/programs/build/Minimac4/build/minimac4
TMPROOT=/u/project/cluo_scratch/terencew/claude/latent_genos/A04c

LIST=$PROJDIR/txt/latent_genos_souporcell_tasks.txt
# ID=22
ID=$SGE_TASK_ID
ROW=$(( (ID - 1) / 22 + 1 ))
C=$(( (ID - 1) % 22 + 1 ))
CHR=chr$C
read -r TREE POOL K MOD <<< "$(sed -n "${ROW}p" "$LIST")"

IMP=$PROJDIR/$TREE/$POOL/impute/$MOD
KEEP=$IMP/loo/keep.txt
EXCL=$IMP/loo/exclude.txt
MSAV=$PANEL/1000G_30x.$CHR.msav
IN=$SRC/1000G.$CHR.vcf.gz
[ "$CHR" = chr1 ] && IN=$CHR1
TMP=$TMPROOT/$TREE.$POOL.$MOD.$CHR
N=${NSLOTS:-1}

echo "Start: $(date)  host=$(hostname -s)  tree=$TREE pool=$POOL mod=$MOD chr=$CHR slots=$N"

if [ ! -s "$IMP/target.all.vcf.gz.tbi" ]; then
    echo "$(date): A04b targets missing ($IMP), skipping"
    exit 0
fi
[ -s "$MSAV" ] || { echo "ERROR: $MSAV missing (run A04a)" >&2; exit 1; }
echo "leave-pool-out: excluding $(wc -l < "$EXCL") samples, keeping $(wc -l < "$KEEP")"

mkdir -p "$TMP"

for QC in server all; do
    OUT=$IMP/$QC
    mkdir -p "$OUT"
    DOSE=$OUT/$CHR.dose.vcf.gz
    if [ -s "$DOSE.tbi" ]; then
        echo "$(date): [$QC] $DOSE exists, skipping"
        continue
    fi

    TGT=$TMP/$QC.target.vcf.gz
    bcftools view -r "$CHR" "$IMP/target.$QC.vcf.gz" -Oz -o "$TGT"
    bcftools index -t -f "$TGT"
    NT=$(bcftools index -n "$TGT")
    echo "$(date): [$QC] $NT target sites on $CHR"
    if [ "$NT" -lt 3 ]; then
        echo "$(date): [$QC] fewer than 3 target sites, skipping $CHR"
        continue
    fi

    # per-pool leave-out Eagle reference, restricted to the target's sites
    REF=$TMP/$QC.ref.bcf
    time bcftools view --threads "$N" -S "^$EXCL" --force-samples -T "$TGT" \
        -m2 -M2 -v snps "$IN" -Ob -o "$REF"
    bcftools index -f "$REF"

    time $EAGLE --vcfRef "$REF" --vcfTarget "$TGT" --geneticMapFile "$MAP" \
        --outPrefix "$TMP/$QC.phased" --vcfOutFormat z --numThreads "$N"
    mv "$TMP/$QC.phased.vcf.gz" "$OUT/$CHR.phased.vcf.gz"
    bcftools index -t -f "$OUT/$CHR.phased.vcf.gz"

    time $MINIMAC4 "$MSAV" "$OUT/$CHR.phased.vcf.gz" \
        --sample-ids-file "$KEEP" \
        -f GT,DS,HDS,GP,SD -a \
        -e "$OUT/$CHR.empirical.vcf.gz" \
        -o "$TMP/$QC.dose.vcf.gz" -O vcf.gz \
        -t "$N" --temp-prefix "$TMP/m4_" --min-ratio-behavior skip
    mv "$TMP/$QC.dose.vcf.gz" "$DOSE"
    bcftools index -t -f "$DOSE"
    echo "$(date): [$QC] imputed $(bcftools index -n "$DOSE") sites"
    rm -f "$TGT" "$TGT.tbi" "$REF" "$REF.csi"
done

rmdir "$TMP" 2>/dev/null || true
echo "End: $(date)  tree=$TREE pool=$POOL mod=$MOD chr=$CHR"
