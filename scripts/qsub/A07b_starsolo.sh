#!/bin/bash
#$ -N A07b_starsolo
#$ -cwd
#$ -l h_data=8G,h_rt=24:00:00
#$ -pe shared 8
#$ -t 1-5:1
#$ -tc 5
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07b_starsolo.$JOB_ID.$TASK_ID
#$ -j y

# OneK1K pilot: STARsolo on one pool (row of txt/onek1k_pilot_pools.txt: <pool> <GSM> <n_donors>),
# all 20 lanes together. 10x 3' v2 (737K-august-2016 whitelist, CB 16 + UMI 10), CellRanger-like
# settings (1MM_multi CB matching, 1MM_CR UMI dedup, EmptyDrops_CR cell calling, CellRanger4 clipping).
# Reference: hg38_igvf STAR index (GRCh38 no-alt analysis set + GENCODE v43, built with STAR 2.7.4a);
# STAR 2.7.10a from ~/bin (the "mapping" env has kb but no STAR).
# BAM: coordinate-sorted, CB/UB tags for souporcell (A07d).
# Output: $SCR/A07/starsolo/pool<P>/{Aligned.sortedByCoord.out.bam,Solo.out/...}
# Requires A07a for this pool (all 20 runs) and A07c (whitelist).
# Submit with (from latent_genos/): qsub -hold_jid A07a_onek1k_download,A07c_onek1k_truth scripts/qsub/A07b_starsolo.sh

source ~/.bashrc
set -euo pipefail
module load samtools 2>/dev/null || true

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos/A07
STAR=/u/home/t/terencew/bin/STAR
INDEX=/u/project/cluo/terencew/reference/hg38_igvf/STAR
WL=$PROJ/reference/10x/737K-august-2016.txt
# ID=1
ID=$SGE_TASK_ID
read -r POOL GSM NDON <<< "$(sed -n "${ID}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
FQ=$SCR/fastq/pool$POOL
OUT=$SCR/starsolo/pool$POOL
N=${NSLOTS:-1}

echo "Start: $(date)  host=$(hostname -s)  pool=$POOL donors=$NDON slots=$N"
if [ -s "$OUT/Aligned.sortedByCoord.out.bam.bai" ]; then echo "$(date): done, skipping"; exit 0; fi
RUNS=$(awk -v p="$POOL" '$1==p {print $2}' "$PROJ/txt/onek1k_pilot_runs.txt")
for r in $RUNS; do [ -s "$FQ/${r}_R2.fastq.gz" ] || { echo "ERROR: missing $FQ/${r}_R2.fastq.gz" >&2; exit 1; }; done
[ -s "$WL" ] || { echo "ERROR: whitelist missing $WL (A07c)" >&2; exit 1; }
R1=$(for r in $RUNS; do printf "%s," "$FQ/${r}_R1.fastq.gz"; done | sed 's/,$//')
R2=$(for r in $RUNS; do printf "%s," "$FQ/${r}_R2.fastq.gz"; done | sed 's/,$//')
mkdir -p "$OUT"

time $STAR --runThreadN "$N" --genomeDir "$INDEX" \
    --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
    --soloType CB_UMI_Simple --soloCBwhitelist "$WL" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloBarcodeReadLength 0 \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR --soloFeatures Gene \
    --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
    --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --limitBAMsortRAM 30000000000 --outBAMsortingThreadN "$N" \
    --outTmpDir "${TMPDIR:-/tmp}/A07b.pool$POOL" --outFileNamePrefix "$OUT/"
time samtools index -@ "$N" "$OUT/Aligned.sortedByCoord.out.bam"
echo "$(date): cells $(wc -l < "$OUT/Solo.out/Gene/filtered/barcodes.tsv")"
grep -E "Number of input reads|Uniquely mapped reads %|Reads With Valid Barcodes|Estimated Number of Cells|Median UMI per Cell" \
    "$OUT/Log.final.out" "$OUT/Solo.out/Gene/Summary.csv" || true
echo "End: $(date)  pool=$POOL"
