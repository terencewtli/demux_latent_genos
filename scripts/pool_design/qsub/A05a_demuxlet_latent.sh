#!/bin/bash
#$ -N A05a_demuxlet_latent
#$ -cwd
#$ -l h_data=8G,h_rt=12:00:00
#$ -pe shared 4
#$ -t 1-12:1
#$ -tc 12
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A05a_demuxlet_latent.$JOB_ID.$TASK_ID
#$ -j y

# Step 3 (close the loop): re-demultiplex each souporcell run with demuxlet using genotypes recovered from the
# pooled reads instead of the true donor genotypes.
#
# Task ID -> row of txt/latent_genos_souporcell_tasks.txt (12 rows: tree pool k mod).
# Arms (ARMS, default all three; samples are souporcell clusters c0..c{k-1}, mapped to donors at scoring
# time with A03b assign.tsv):
#   latent_GT   souporcell cluster_genotypes.vcf hard calls (no imputation)
#   imputed_GT  GLIMPSE2 ambient_b2m5 (A04g) hard calls
#   imputed_GP  GLIMPSE2 ambient_b2m5 genotype posteriors (demuxlet --field GP)
# Baseline = the existing pool_design demuxlet .best on TRUE genotypes (same pileup, same sites).
#
# Sites: every arm is restricted (exact CHROM/POS/REF/ALT) to the pool VCF sites, i.e. the sites of the existing
# popscle pileup = 1000G global AF >= 5% biallelic SNPs (B01). That list does not depend on who is in the pool
# (14% of chr22 sites are monomorphic in the n16 random rep1 pool), so no arm gets truth-selected sites. The
# latent arm only has genotypes where souporcell's clusters had reads; demuxlet uses the overlap.
#
# Reuses <tree>/<pool>/demux/demuxlet/<mod>/<SAMPLE>.pileup.* (no new pileup). ambisim_final greedy_maxmin rep1
# has no pileup (never demuxlet'd in pool_design); that task exits with a message.
#
# Output: latent_genos/results/A05a_demuxlet_latent/<tree>__<pool>__<mod>/{<arm>.vcf.gz,<arm>.best,...}
# Submit: qsub scripts/latent_genos/qsub/A05a_demuxlet_latent.sh   (one arm: qsub -v ARMS=imputed_GP ...)

source ~/.bashrc
set -euo pipefail
module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
LGDIR=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SAMPLE=20220928-IGVF-D0
POPSCLE=/u/project/cluo/terencew/programs/popscle/bin/popscle
ARMS=${ARMS:-latent_GT imputed_GT imputed_GP}
ARMS=${ARMS//,/ }

ID=${SGE_TASK_ID:-1}
read TREE POOL K MOD < <(sed -n "${ID}p" $PROJDIR/txt/latent_genos_souporcell_tasks.txt)
TAG=${TREE}__${POOL}__${MOD}
OUT=$LGDIR/results/A05a_demuxlet_latent/$TAG; mkdir -p $OUT
PLP=$PROJDIR/$TREE/$POOL/demux/demuxlet/$MOD/$SAMPLE.pileup
BARCODES=$PROJDIR/$TREE/$POOL/cr_arc/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
POOLVCF=$PROJDIR/$TREE/vcf/$POOL.vcf.gz
SOUP=$PROJDIR/$TREE/$POOL/demux/souporcell/$MOD/$SAMPLE/cluster_genotypes.vcf
IMP=$LGDIR/results/A04g_glimpse_impute/ambient_b2m5/$TAG
echo "$(date): $TAG (k=$K) arms: $ARMS on $(hostname -s)"

if [ ! -s $PLP.var.gz ]; then echo "no pileup for $TAG ($PLP.var.gz); skipping"; exit 0; fi
TMP=${TMPDIR:-/tmp}/A05a.$TAG.$$; mkdir -p $TMP; trap 'rm -rf $TMP' EXIT

# pool sites (sites-only), shared by both builds
if [ ! -s $OUT/pool_sites.vcf.gz ]; then
    time bcftools view -G -Oz -o $TMP/s.vcf.gz $POOLVCF && bcftools index -t $TMP/s.vcf.gz
    mv $TMP/s.vcf.gz $OUT/pool_sites.vcf.gz; mv $TMP/s.vcf.gz.tbi $OUT/pool_sites.vcf.gz.tbi
fi

# exact-allele intersection with the pool sites; keep FORMAT from the first file
restrict() {  # in.vcf.gz out.vcf.gz
    bcftools isec -c none -n =2 -w 1 -Oz -o $TMP/r.vcf.gz $1 $OUT/pool_sites.vcf.gz
    mv $TMP/r.vcf.gz $2; bcftools index -t -f $2
}

if [[ " $ARMS " == *" latent_GT "* ]] && [ ! -s $OUT/latent_GT.vcf.gz ]; then
    # souporcell names its samples 0..k-1 (renamed c0.. to match A03b / A04g) and leaves its FORMAT tags undeclared,
    # which bcftools 1.11 cannot subset; keep GT only, in awk
    time awk -v OFS='\t' '/^##FORMAT/{next} /^##/{print;next}
        /^#CHROM/{print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
                  print "##FILTER=<ID=BACKGROUND,Description=\"souporcell background\">"
                  for(i=10;i<=NF;i++)$i="c"$i; print; next}
        {$9="GT"; for(i=10;i<=NF;i++){split($i,a,":"); $i=a[1]} print}' $SOUP |
        bcftools view -m2 -M2 -v snps -Oz -o $TMP/l.vcf.gz -
    bcftools index -t $TMP/l.vcf.gz
    restrict $TMP/l.vcf.gz $OUT/latent_GT.vcf.gz
fi
if [[ " $ARMS " == *" imputed_"* ]] && [ ! -s $OUT/imputed.vcf.gz ]; then
    ls $IMP/chr{1..22}.imputed.bcf > $TMP/chroms.txt
    time bcftools concat -f $TMP/chroms.txt --naive -Ob -o $TMP/all.bcf
    bcftools view -m2 -M2 -v snps -Oz -o $TMP/i.vcf.gz $TMP/all.bcf && bcftools index -t $TMP/i.vcf.gz
    restrict $TMP/i.vcf.gz $OUT/imputed.vcf.gz
fi
for f in $OUT/*.vcf.gz; do echo "$(basename $f): $(bcftools index -n $f) sites"; done

for ARM in $ARMS; do
    if [ -s $OUT/$ARM.best ]; then echo "$ARM done already"; continue; fi
    case $ARM in
        latent_GT)  VCF=$OUT/latent_GT.vcf.gz; FIELD=GT ;;
        imputed_GT) VCF=$OUT/imputed.vcf.gz;   FIELD=GT ;;
        imputed_GP) VCF=$OUT/imputed.vcf.gz;   FIELD=GP ;;
        *) echo "unknown arm $ARM" >&2; exit 1 ;;
    esac
    echo "$(date): demuxlet $ARM (--field $FIELD)"
    time $POPSCLE demuxlet --plp $PLP --group-list $BARCODES --field $FIELD --vcf $VCF --out $OUT/$ARM.tmp
    mv $OUT/$ARM.tmp.best $OUT/$ARM.best
done
echo "$(date): done"
