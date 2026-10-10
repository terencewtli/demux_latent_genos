#!/bin/bash
#$ -N A09b_soup_cluster_cputest
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09b_soup_cluster_cputest.$JOB_ID
#$ -j y
# CPU smoke test of souporcell_gpu cluster + troublet on real counts (ys3a GEX, 10 restarts; submit with -hold_jid A09a_soup_counts), so the one GPU job
# (A09b_soup_cluster_gpu) is not the first time this code sees real data. Writes to $SCR/cputest only.
set -euo pipefail
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos/A09
O=$SCR/souporcell_v2/ys3a_gex
T=$SCR/cputest
PY=/u/home/t/terencew/project-cluo/miniconda3/envs/soupgpu/bin/python
export PYTHONPATH=/u/project/cluo/terencew/claude/exploration/engineer/code/souporcell_gpu/src
export OMP_NUM_THREADS=${NSLOTS:-1}
mkdir -p "$T"
echo "Start: $(date)  host=$(hostname -s)"
time $PY -m souporcell_gpu cluster -r "$O/ref.mtx" -a "$O/alt.mtx" -b "$O/barcodes.tsv" -k 4 --min_ref 10 --min_alt 10 \
    --restarts 10 --device cpu --restart_log "$T/ys3a_gex.restarts.tsv" -o "$T/ys3a_gex.clusters_tmp.tsv"
time $PY -m souporcell_gpu troublet -r "$O/ref.mtx" -a "$O/alt.mtx" -c "$T/ys3a_gex.clusters_tmp.tsv" --device cpu \
    -o "$T/ys3a_gex.clusters.tsv"
awk 'NR>1 {n[$2" "$3]++} END {for (s in n) print "  " s, n[s]}' "$T/ys3a_gex.clusters.tsv" | sort | head -12
echo "End: $(date)"
