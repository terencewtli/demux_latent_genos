#!/bin/bash
#$ -N A04b_prep_targets
#$ -cwd
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 2
#$ -t 1-12:1
#$ -tc 12
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04b_prep_targets.$JOB_ID.$TASK_ID
#$ -j y

# Turn one souporcell run (a row of txt/latent_genos_souporcell_tasks.txt) into
# imputation targets, plus the leave-pool-out sample lists A04c uses.
#   1. loo/{exclude,keep}.txt: pool donors + their 1000G relatives vs the panel
#   2. target.raw.vcf.gz: lib/souporcell_to_target.py (GT + AD + GL, BACKGROUND
#      loci dropped)
#   3. target.matched.vcf.gz: sites whose CHROM/POS/REF/ALT exactly match the
#      local 1000G 30x panel (A04a sites files); mismatches are dropped, not flipped
#   4. two QC sets, applied to the matched sites:
#        server: CR >= 0.9 and not monomorphic, mimicking the TOPMed server's
#                SNP call-rate < 90% and monomorphic exclusions
#        all:    not monomorphic, with no call-rate filter
# Counts at each step: impute/<mod>/target_counts.tsv
#
# Requires A03a (souporcell) and A04a (panel sites) to be finished.
# Submit with: qsub scripts/latent_genos/qsub/A04b_prep_targets.sh

source ~/.bashrc

conda activate allcools

set -euo pipefail

module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
LIB=$PROJDIR/scripts/latent_genos/lib
PANEL=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
PED=$PANEL/20130606_g1k_3202_samples_ped_population.txt
PED_URL=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

LIST=$PROJDIR/txt/latent_genos_souporcell_tasks.txt
# ID=1
ID=$SGE_TASK_ID
read -r TREE POOL K MOD <<< "$(sed -n "${ID}p" "$LIST")"

SOUP=$PROJDIR/$TREE/$POOL/demux/souporcell/$MOD/$SAMPLE/cluster_genotypes.vcf
POOL_VCF=$PROJDIR/$TREE/vcf/$POOL.vcf.gz
OUT=$PROJDIR/$TREE/$POOL/impute/$MOD

echo "Start: $(date)  tree=$TREE pool=$POOL k=$K mod=$MOD"

if [ ! -s "$SOUP" ]; then
    echo "$(date): souporcell output missing ($SOUP), skipping (run A03a first)"
    exit 0
fi
if [ -s "$OUT/target.all.vcf.gz.tbi" ] && [ -s "$OUT/target.server.vcf.gz.tbi" ]; then
    echo "$(date): targets exist, skipping ($OUT)"
    exit 0
fi

mkdir -p "$OUT/loo" "$PANEL"

if [ ! -s "$PED" ]; then
    curl -fsL "$PED_URL" -o "$PED.tmp.$ID" && mv "$PED.tmp.$ID" "$PED"
fi

### 1. leave-pool-out lists
bcftools query -l "$POOL_VCF" > "$OUT/loo/donors.txt"
bcftools query -l "$SRC/1000G.chr22.vcf.gz" > "$OUT/loo/panel_samples.txt"
python "$LIB/leave_pool_out.py" --donors "$OUT/loo/donors.txt" \
    --panel-samples "$OUT/loo/panel_samples.txt" --ped "$PED" --outdir "$OUT/loo"

### 2. convert souporcell calls
time python "$LIB/souporcell_to_target.py" "$SOUP" 2> "$OUT/convert.log" \
    | bcftools sort -Oz -o "$OUT/target.raw.vcf.gz"
bcftools index -t -f "$OUT/target.raw.vcf.gz"
cat "$OUT/convert.log"

### 3. exact allele match to the panel, per chromosome
PARTS=()
for C in $(seq 1 22); do
    S=$PANEL/1000G_30x.chr$C.sites.vcf.gz
    [ -s "$S.tbi" ] || { echo "ERROR: panel sites missing: $S (run A04a)" >&2; exit 1; }
    P=$OUT/target.matched.chr$C.vcf.gz
    bcftools isec -c none -n =2 -w 1 -r chr$C "$OUT/target.raw.vcf.gz" "$S" -Oz -o "$P"
    PARTS+=("$P")
done
bcftools concat "${PARTS[@]}" -Oz -o "$OUT/target.matched.vcf.gz"
bcftools index -t -f "$OUT/target.matched.vcf.gz"
rm -f "${PARTS[@]}"

### 4. QC sets
bcftools view -i 'INFO/CR>=0.9 && INFO/MONO=0' "$OUT/target.matched.vcf.gz" -Oz -o "$OUT/target.server.vcf.gz"
bcftools view -i 'INFO/MONO=0' "$OUT/target.matched.vcf.gz" -Oz -o "$OUT/target.all.vcf.gz"
bcftools index -t -f "$OUT/target.server.vcf.gz"
bcftools index -t -f "$OUT/target.all.vcf.gz"

{
    printf "step\tn_sites\n"
    for S in raw matched server all; do
        printf "%s\t%s\n" "$S" "$(bcftools index -n "$OUT/target.$S.vcf.gz")"
    done
} > "$OUT/target_counts.tsv"
cat "$OUT/target_counts.tsv"

echo "End: $(date)  tree=$TREE pool=$POOL mod=$MOD"
