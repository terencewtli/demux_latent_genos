#!/bin/bash
#$ -N A09b_soup_cluster_gpu
#$ -cwd
#$ -l gpu,cuda=1,h_data=4G,h_rt=2:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09b_soup_cluster_gpu.$JOB_ID
#$ -j y
# souporcell_gpu clustering (200 restarts, k = 4) + troublet for all 20 A09 runs, one GPU job (one queue wait;
# ~1-2 min per run on a 2080 Ti). Skips runs already done or whose counts are missing.
# Submit: qsub -hold_jid A09a_soup_counts scripts/qsub/A09b_soup_cluster_gpu.sh
set -euo pipefail
PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos/A09/souporcell_v2
FASTA=/u/project/cluo/terencew/reference/hg38_igvf/IGVF_hg38/fasta/genome.fa
COMMON=/u/project/cluo/terencew/reference/cellsnp_db/hg38.AF5e2.possorted.chr_prefix.vcf
PY=/u/home/t/terencew/project-cluo/miniconda3/envs/soupgpu/bin/python
export PYTHONPATH=/u/project/cluo/terencew/claude/exploration/engineer/code/souporcell_gpu/src
echo "Start: $(date)  host=$(hostname -s)"; nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
FAIL=0
while read -r POOL MOD; do
    O=$SCR/${POOL}_${MOD}
    [ -e "$O/vartrix.done" ] || { echo "SKIP $POOL $MOD: no counts"; continue; }
    BC=$O/barcodes.tsv
    time $PY -m souporcell_gpu pipeline -i NA -b "$BC" -f "$FASTA" -k 4 -o "$O" --common_variants "$COMMON" \
        --restarts 200 --steps cluster,troublet --device cuda || { echo "FAILED $POOL $MOD"; FAIL=$((FAIL + 1)); continue; }
    awk 'NR>1 {n[$2" "$3]++} END {for (s in n) print "  " s, n[s]}' "$O/clusters.tsv" | sort | head -12
done < "$PROJ/txt/A09_runs.txt"
echo "End: $(date)  failed runs: $FAIL"
[ "$FAIL" -eq 0 ]
