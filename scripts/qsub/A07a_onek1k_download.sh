#!/bin/bash
#$ -N A07a_onek1k_download
#$ -cwd
#$ -l h_data=4G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-100:1
#$ -tc 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07a_onek1k_download.$JOB_ID.$TASK_ID
#$ -j y

# OneK1K pilot (A07): download one SRA run (row of txt/onek1k_pilot_runs.txt: <pool> <SRR>)
# and write gzipped 10x 3' v2 R1 (26 bp: 16 CB + 10 UMI) and R2 (98 bp cDNA) to scratch.
# fasterq-dump --include-technical --split-files gives _1 = I1 (8 bp), _2 = R1, _3 = R2
# (checked on SRR18028385); read lengths are re-checked here. The .sra and raw fastqs
# live in the node-local $TMPDIR (removed by SGE at job end); only R1/R2 .gz are kept
# (~3.2 GB per run, ~64 GB per pool).
# Output: $SCR/A07/fastq/pool<P>/<SRR>_R{1,2}.fastq.gz
# Submit with (from latent_genos/): qsub scripts/qsub/A07a_onek1k_download.sh

source ~/.bashrc
conda activate sratools
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos/A07
LIST=$PROJ/txt/onek1k_pilot_runs.txt
# ID=1
ID=$SGE_TASK_ID
read -r POOL RUN <<< "$(sed -n "${ID}p" "$LIST")"
OUT=$SCR/fastq/pool$POOL
TMP=${TMPDIR:-/tmp}/A07a.$RUN
N=${NSLOTS:-1}
mkdir -p "$OUT" "$TMP"

echo "Start: $(date)  host=$(hostname -s)  pool=$POOL run=$RUN  tmp=$TMP ($(df -h "$TMP" | awk 'NR==2{print $4}') free)"
if [ -s "$OUT/${RUN}_R1.fastq.gz" ] && [ -s "$OUT/${RUN}_R2.fastq.gz" ]; then
    echo "$(date): exists, skipping"; exit 0
fi

cd "$TMP"
time prefetch --max-size 50G -O "$TMP" "$RUN"
time fasterq-dump --include-technical --split-files -e "$N" -t "$TMP" -O "$TMP" "$TMP/$RUN"
for i in 1 2 3; do echo "${RUN}_$i: $(sed -n 2p "${RUN}_$i.fastq" | tr -d '\n' | wc -c) bp"; done
[ "$(sed -n 2p "${RUN}_2.fastq" | tr -d '\n' | wc -c)" -eq 26 ] || { echo "ERROR: _2 is not the 26 bp R1" >&2; exit 1; }
time pigz -p "$N" -c "${RUN}_2.fastq" > "$OUT/${RUN}_R1.fastq.gz.part"
time pigz -p "$N" -c "${RUN}_3.fastq" > "$OUT/${RUN}_R2.fastq.gz.part"
mv "$OUT/${RUN}_R1.fastq.gz.part" "$OUT/${RUN}_R1.fastq.gz"
mv "$OUT/${RUN}_R2.fastq.gz.part" "$OUT/${RUN}_R2.fastq.gz"
echo "$(date): reads $(( $(wc -l < "${RUN}_2.fastq") / 4 ))  $(du -ch "$OUT/${RUN}"_R*.fastq.gz | tail -1)"
echo "End: $(date)  pool=$POOL run=$RUN"
