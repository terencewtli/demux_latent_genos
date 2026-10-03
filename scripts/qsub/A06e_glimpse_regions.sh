#!/bin/bash
#$ -N A06e_glimpse_regions
#$ -cwd
#$ -l h_data=4G,h_rt=8:00:00
#$ -pe shared 4
# GLIMPSE2 static binaries need AVX2
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-100
#$ -tc 40
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A06e_glimpse_regions.$JOB_ID.$TASK_ID
#$ -j y

# A06e: Level B imputation for one A06 region (docs/A06_FINEMAP_PLAN.md).
#   1. leave-out reference: 1000G 30x minus the 300 test donors and their relatives (<a06a>/loo/exclude.txt; 2,814
#      kept), biallelic SNVs + indels, AC >= 1, region TSS +- (half + buffer)
#   2. GLIMPSE2_split_reference: input region = TSS +- 1 Mb, output region = TSS +- 500 kb (one chunk per region)
#   3. per arm (gex, atac, multiome): A06d target allele-matched to the panel sites (as A04g), GLIMPSE2_phase
#      --burnin 2 --main 5 (the A04g b2m5 setting), then A06e_extract_ds.py -> <region>/imp_<arm>.tsv.gz
# Task = row of <a06a>/regions.tsv. Usage: qsub [-t ...] scripts/qsub/A06e_glimpse_regions.sh <a06a_outdir>
source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
SITES=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
GBIN=$SCR/bin
A06A=$1
N=${NSLOTS:-4}
HALF=500000
BUF=500000
# ID=1
ID=$SGE_TASK_ID
read -r CHR TSS GENE START END REGION <<< "$(awk -F'\t' -v i=$((ID + 1)) 'NR==i{print $1, $2, $3, $4, $5, $6}' $A06A/regions.tsv)"
RD=$A06A/$REGION
TMP=$SCR/A06e/$REGION
mkdir -p $TMP
echo "$(date): A06e $REGION $GENE $CHR:$TSS on $(hostname -s), $N slots"
if [ -s $RD/imp_multiome.tsv.gz ] && [ -s $RD/imp_atac.tsv.gz ] && [ -s $RD/imp_gex.tsv.gz ]; then echo "done already"; exit 0; fi

IS=$(( TSS - HALF - BUF )); [ $IS -lt 1 ] && IS=1
IE=$(( TSS + HALF + BUF ))
OS=$(( TSS - HALF )); [ $OS -lt 1 ] && OS=1
OE=$(( TSS + HALF ))
BIN=$(ls $TMP/ref_${CHR}_*.bin 2>/dev/null | head -1 || true)   # GLIMPSE2 names it after the INPUT region
if [ -z "$BIN" ] || [ ! -s "$BIN" ]; then
    time bcftools view --threads $N -r $CHR:$IS-$IE -S ^$A06A/loo/exclude.txt --force-samples \
        -m2 -M2 -v snps,indels -c1 $SRC/1000G.$CHR.vcf.gz -Ob -o $TMP/ref.bcf
    bcftools index -f $TMP/ref.bcf
    time $GBIN/GLIMPSE2_split_reference_static --reference $TMP/ref.bcf --map $SCR/glimpse_maps/b38/$CHR.b38.gmap.gz \
        --input-region $CHR:$IS-$IE --output-region $CHR:$OS-$OE --output $TMP/ref --threads $N > $TMP/split.log
    rm -f $TMP/ref.bcf $TMP/ref.bcf.csi
    BIN=$(ls $TMP/ref_${CHR}_*.bin | head -1)
fi
[ -s "$BIN" ] || { echo "ERROR: missing $BIN" >&2; ls $TMP >&2; exit 1; }

for ARM in gex atac multiome; do
    OUT=$RD/imp_$ARM.tsv.gz
    [ -s $OUT ] && { echo "$ARM done already"; continue; }
    if [ -e $RD/target_$ARM.none ]; then
        echo "$ARM: no target sites, writing no-information dosages"
        python $PROJDIR/scripts/A06e_extract_ds.py NONE $A06A $REGION $OUT
        continue
    fi
    bcftools index -f $RD/target_$ARM.vcf.gz
    bcftools isec -c none -n =2 -w 1 $RD/target_$ARM.vcf.gz $SITES/1000G_30x.$CHR.sites.vcf.gz -Ob -o $TMP/tgt_$ARM.bcf
    bcftools index -f $TMP/tgt_$ARM.bcf
    echo "$(date): $ARM target sites $(bcftools index -n $TMP/tgt_$ARM.bcf)"
    time $GBIN/GLIMPSE2_phase_static --input-gl $TMP/tgt_$ARM.bcf --reference $BIN --impute-reference-only-variants \
        --burnin 2 --main 5 --threads $N --output $TMP/imp_$ARM.bcf --log $TMP/imp_$ARM.log > /dev/null
    bcftools index -f $TMP/imp_$ARM.bcf
    python $PROJDIR/scripts/A06e_extract_ds.py $TMP/imp_$ARM.bcf $A06A $REGION $OUT
done
echo "$(date): done"
