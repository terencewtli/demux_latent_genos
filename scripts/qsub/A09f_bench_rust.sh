#!/bin/bash
#$ -N A09f_bench_rust
#$ -cwd
#$ -l h_data=4G,h_rt=24:00:00
#$ -pe shared 8
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09f_bench_rust.$JOB_ID
#$ -j y
# A09f equal-restart benchmark, arm 1 of 3: original Rust souporcell 2.4 clustering + troublet on ys3j ATAC
# (largest A09 run: 22,079 cells), same settings as A09b (k=4, 200 restarts, min_ref/min_alt 10, seed 4), 8 threads.
# Reuses the A09a matrices. Writes $SCR/bench/rust only. Compare with scripts/A09f_bench_compare.py.
# Submit with (from latent_genos/): qsub scripts/qsub/A09f_bench_rust.sh
set -euo pipefail
O=/u/project/cluo_scratch/terencew/claude/latent_genos/A09/souporcell_v2/ys3j_atac
T=/u/project/cluo_scratch/terencew/claude/latent_genos/A09/bench/rust
SP=/u/project/cluo/terencew/programs/souporcell
N=${NSLOTS:-1}
mkdir -p "$T"
echo "Start: $(date)  host=$(hostname -s)  threads=$N  cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
echo "### cluster"
time $SP/souporcell/target/release/souporcell -r "$O/ref.mtx" -a "$O/alt.mtx" -b "$O/barcodes.tsv" -k 4 \
    --min_ref 10 --min_alt 10 --restarts 200 --seed 4 -t "$N" > "$T/clusters_tmp.tsv" 2> "$T/cluster.err"
tail -3 "$T/cluster.err"
echo "### troublet"
time $SP/troublet/target/release/troublet -r "$O/ref.mtx" -a "$O/alt.mtx" -c "$T/clusters_tmp.tsv" \
    > "$T/clusters.tsv" 2> "$T/troublet.err"
awk 'NR>1 {n[$2]++} END {for (s in n) print "  " s, n[s]}' "$T/clusters.tsv"
echo "End: $(date)"
