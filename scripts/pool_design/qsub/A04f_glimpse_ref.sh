#!/bin/bash
#$ -N A04f_glimpse_ref
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
# GLIMPSE2 static binaries need AVX2 (Illegal instruction on older nodes, 2026-10-01)
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-132:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04f_glimpse_ref.$JOB_ID.$TASK_ID
#$ -j y

# Per-pool leave-pool-out GLIMPSE2 reference panel for one chromosome.
# GLIMPSE2 cannot drop reference samples at run time (unlike minimac4
# --sample-ids-file), so each pool gets its own binary panel:
#   1. loo/{exclude,keep}.txt: pool donors + 1000G relatives (lib/leave_pool_out.py,
#      same rule as A04b)
#   2. ref.bcf: 1000G 30x phased panel minus exclude.txt; biallelic SNVs + indels
#      (same filter as A04a), AC >= 1 after the sample drop
#   3. GLIMPSE2_chunk --sequential (default windows) -> chunks.txt
#   4. GLIMPSE2_split_reference per chunk -> bin/ref_<chr>_<start>_<end>.bin
#   ref.bcf is deleted after the split; only the binaries are kept.
#
# Task ID -> unique (tree, pool) of txt/latent_genos_souporcell_tasks.txt (6 pools;
# GEX and ATAC share a panel) x chr1-22:
#   pool = (ID-1)/22 + 1, chr = (ID-1)%22 + 1
#   e.g. ambisim_final EUR_only__greedy_maxmin__rep1 chr20 = 130
# chr1: the by_chrom copy is corrupt (gzip crc error), so chr1 comes from the clean
# EBI re-download (A04h_download_chr1.sh, latent_genos/reference/1000G_30x/); wait
# for A04h before submitting chr1 tasks (1, 23, 45, 67, 89, 111).
#
# Output (scratch): $SCR/glimpse_ref/<tree>__<pool>/chrN/{loo/,chunks.txt,bin/,split.done}
# Binaries: executable copies of programs/GLIMPSE2/*_static (the originals lack +x)
# in $SCR/bin. Maps: GLIMPSE b38 maps (github odelaneau/GLIMPSE maps/genetic_maps.b38).
#
# Submit with: qsub scripts/latent_genos/qsub/A04f_glimpse_ref.sh

source ~/.bashrc

conda activate allcools

set -euo pipefail

module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
LIB=$PROJDIR/scripts/latent_genos/lib
SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
CHR1=/u/project/cluo/terencew/claude/project_ideas/latent_genos/reference/1000G_30x/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
PED=/u/project/cluo/terencew/reference/topmed/local_1000G_30x/20130606_g1k_3202_samples_ped_population.txt
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
GBIN=$SCR/bin
MAPS=$SCR/glimpse_maps/b38

LIST=$PROJDIR/txt/latent_genos_souporcell_tasks.txt
# ID=130
ID=$SGE_TASK_ID
PIDX=$(( (ID - 1) / 22 + 1 ))
C=$(( (ID - 1) % 22 + 1 ))
CHR=chr$C
read -r TREE POOL <<< "$(awk '!seen[$1" "$2]++ {print $1, $2}' "$LIST" | sed -n "${PIDX}p")"

POOL_VCF=$PROJDIR/$TREE/vcf/$POOL.vcf.gz
MAP=$MAPS/$CHR.b38.gmap.gz
IN=$SRC/1000G.$CHR.vcf.gz
[ "$CHR" = chr1 ] && IN=$CHR1
OUT=$SCR/glimpse_ref/${TREE}__${POOL}/$CHR
N=${NSLOTS:-1}

echo "Start: $(date)  host=$(hostname -s)  tree=$TREE pool=$POOL chr=$CHR slots=$N"

if [ -e "$OUT/split.done" ]; then
    echo "$(date): $OUT/split.done exists, skipping"
    exit 0
fi
[ -s "$MAP" ] || { echo "ERROR: genetic map missing: $MAP" >&2; exit 1; }
[ -s "$IN.tbi" ] || { echo "ERROR: source VCF/index missing: $IN" >&2; exit 1; }
mkdir -p "$OUT/loo" "$OUT/bin"

### 1. leave-pool-out lists
bcftools query -l "$POOL_VCF" > "$OUT/loo/donors.txt"
bcftools query -l "$SRC/1000G.chr22.vcf.gz" > "$OUT/loo/panel_samples.txt"
python "$LIB/leave_pool_out.py" --donors "$OUT/loo/donors.txt" \
    --panel-samples "$OUT/loo/panel_samples.txt" --ped "$PED" --outdir "$OUT/loo"

### 2. leave-out reference BCF
BCF=$OUT/ref.bcf
time bcftools view --threads "$N" -S "^$OUT/loo/exclude.txt" --force-samples \
    -m2 -M2 -v snps,indels -c1 "$IN" -Ob -o "$BCF"
bcftools index -f "$BCF"
echo "$(date): ref $(bcftools index -n "$BCF") records, $(bcftools query -l "$BCF" | wc -l) samples"

### 3. chunks
time "$GBIN/GLIMPSE2_chunk_static" --input "$BCF" --region "$CHR" --map "$MAP" \
    --sequential --threads "$N" --output "$OUT/chunks.txt" --log "$OUT/chunk.log"
echo "$(date): $(wc -l < "$OUT/chunks.txt") chunks"

### 4. binary reference per chunk
while read -r CIDX CCHR IRG ORG REST; do
    time "$GBIN/GLIMPSE2_split_reference_static" --reference "$BCF" --map "$MAP" \
        --input-region "$IRG" --output-region "$ORG" --threads "$N" \
        --output "$OUT/bin/ref" --log "$OUT/bin/split_$CIDX.log"
done < "$OUT/chunks.txt"

rm -f "$BCF" "$BCF.csi"
touch "$OUT/split.done"
du -sh "$OUT/bin"
echo "End: $(date)  tree=$TREE pool=$POOL chr=$CHR"
