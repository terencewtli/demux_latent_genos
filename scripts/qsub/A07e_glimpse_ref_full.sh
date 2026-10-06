#!/bin/bash
#$ -N A07e_glimpse_ref_full
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
# GLIMPSE2 static binaries need AVX2 (Illegal instruction on older nodes, 2026-10-01)
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-22:1
#$ -tc 22
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07e_glimpse_ref_full.$JOB_ID.$TASK_ID
#$ -j y

# Full 1000G 30x GLIMPSE2 reference (all 3,202 samples) for one chromosome. OneK1K donors are not in
# 1000G, so no leave-pool-out is needed (unlike A04f). Otherwise identical to A04f: biallelic SNVs +
# indels, AC >= 1, GLIMPSE2_chunk --sequential, split_reference per chunk; ref.bcf is removed after
# the split (same as A04f, it is this job's own intermediate).
# Output: $SCR/glimpse_ref/onek1k_full_1000G/chrN/{chunks.txt,bin/,split.done}
# Submit with (from latent_genos/): qsub scripts/qsub/A07e_glimpse_ref_full.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
CHR1=/u/project/cluo/terencew/claude/project_ideas/latent_genos/reference/1000G_30x/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
GBIN=$SCR/bin
# ID=22
ID=$SGE_TASK_ID
CHR=chr$ID
MAP=$SCR/glimpse_maps/b38/$CHR.b38.gmap.gz
IN=$SRC/1000G.$CHR.vcf.gz
[ "$CHR" = chr1 ] && IN=$CHR1   # by_chrom chr1 is corrupt (stops at 190.67 Mb)
OUT=$SCR/glimpse_ref/onek1k_full_1000G/$CHR
N=${NSLOTS:-1}

echo "Start: $(date)  host=$(hostname -s)  chr=$CHR slots=$N"
if [ -e "$OUT/split.done" ]; then echo "$(date): exists, skipping"; exit 0; fi
[ -s "$MAP" ] || { echo "ERROR: genetic map missing: $MAP" >&2; exit 1; }
mkdir -p "$OUT/bin"

BCF=$OUT/ref.bcf
time bcftools view --threads "$N" -m2 -M2 -v snps,indels -c1 "$IN" -Ob -o "$BCF"
bcftools index -f "$BCF"
echo "$(date): ref $(bcftools index -n "$BCF") records, $(bcftools query -l "$BCF" | wc -l) samples"
time "$GBIN/GLIMPSE2_chunk_static" --input "$BCF" --region "$CHR" --map "$MAP" \
    --sequential --threads "$N" --output "$OUT/chunks.txt" --log "$OUT/chunk.log"
while read -r CIDX CCHR IRG ORG REST; do
    time "$GBIN/GLIMPSE2_split_reference_static" --reference "$BCF" --map "$MAP" \
        --input-region "$IRG" --output-region "$ORG" --threads "$N" \
        --output "$OUT/bin/ref" --log "$OUT/bin/split_$CIDX.log"
done < "$OUT/chunks.txt"
rm -f "$BCF" "$BCF.csi"
touch "$OUT/split.done"
du -sh "$OUT/bin"
echo "End: $(date)  chr=$CHR"
