#!/bin/bash
#$ -N A09f_bench_cpu
#$ -cwd
#$ -l h_data=4G,h_rt=24:00:00
#$ -pe shared 8
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09f_bench_cpu.$JOB_ID
#$ -j y
# A09f equal-restart benchmark, arm 2 of 3: souporcell_gpu cluster + troublet with --device cpu (8 threads) on ys3j ATAC,
# same settings as the Rust arm (A09f_bench_rust.sh: k=4, 200 restarts, min_ref/min_alt 10, seed 4).
# Writes $SCR/bench/cpu only. Submit with (from latent_genos/): qsub scripts/qsub/A09f_bench_cpu.sh
set -euo pipefail
O=/u/project/cluo_scratch/terencew/claude/latent_genos/A09/souporcell_v2/ys3j_atac
T=/u/project/cluo_scratch/terencew/claude/latent_genos/A09/bench/cpu
PY=/u/home/t/terencew/project-cluo/miniconda3/envs/soupgpu/bin/python
export PYTHONPATH=/u/project/cluo/terencew/claude/exploration/engineer/code/souporcell_gpu/src
export OMP_NUM_THREADS=${NSLOTS:-1}
mkdir -p "$T"
echo "Start: $(date)  host=$(hostname -s)  threads=$OMP_NUM_THREADS  cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
echo "### cluster"
time $PY -m souporcell_gpu cluster -r "$O/ref.mtx" -a "$O/alt.mtx" -b "$O/barcodes.tsv" -k 4 --min_ref 10 --min_alt 10 \
    --restarts 200 --seed 4 --device cpu --restart_log "$T/restarts.tsv" -o "$T/clusters_tmp.tsv"
echo "### troublet"
time $PY -m souporcell_gpu troublet -r "$O/ref.mtx" -a "$O/alt.mtx" -c "$T/clusters_tmp.tsv" --device cpu -o "$T/clusters.tsv"
awk 'NR>1 {n[$2]++} END {for (s in n) print "  " s, n[s]}' "$T/clusters.tsv"
echo "End: $(date)"
