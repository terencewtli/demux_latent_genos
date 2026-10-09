#!/bin/bash
#$ -N A08a_select_eqtls
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 1
#$ -t 1-14:1
#$ -tc 14
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A08a_select_eqtls.$JOB_ID.$TASK_ID
#$ -j y

# OneK1K published cis-eQTLs (onek1k.org S3 full tables, one per cell type) -> significant rows, using the paper's
# own criterion (powellgenomicslab/onek1k_phase1 round*.run_spearman_rank_test.R):
#   Spearman rho per (gene, SNP), cis = gene +/- 1 Mb, imputed SNPs with MAF > 0.05; qvalue() per cell type x chr;
#   significant = localFDR < 0.05 (the table's FDR column); rounds 1-5 = conditional analysis on the previous lead.
# Step 1 (awk, ~2 min per table): keep FDR < 0.05 rows -> results/A08_eqtl_anchor/sig/<ct>.sig.tsv.gz
# Step 2 (A08a_select_eqtls.py, after all 14 tasks): leads + loci table. Run once by hand or with -hold_jid.
# Submit with (from latent_genos/): qsub scripts/qsub/A08a_select_eqtls.sh

source ~/.bashrc
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
IN=$PROJ/reference/onek1k/published_eqtl
OUT=$PROJ/results/A08_eqtl_anchor/sig
CTS=(cd4et cd4nc cd4sox4 cd8nc cd8et cd8s100b nk nkr bmem bin plasma monoc mononc dc)
# ID=1
ID=$SGE_TASK_ID
CT=${CTS[$((ID - 1))]}
F=$IN/${CT}_eqtl_table.tsv.gz
mkdir -p "$OUT"
echo "Start: $(date)  host=$(hostname -s)  ct=$CT"
gzip -t "$F" || { echo "ERROR: $F is incomplete or corrupt" >&2; exit 1; }
# columns: 1 CELL_ID 3 RSID 4 SNPID 5 GENE 6 GENE_ID 7 CHR 8 POS 9 A1 10 A2 11 A2_FREQ_ONEK1K 13 RHO 15 P 16 Q 17 FDR 18 RSQUARE 19 GENOTYPED 20 ROUND
time zcat "$F" | awk -F'\t' 'NR == 1 || $17 < 0.05' | gzip > "$OUT/$CT.sig.tsv.gz.tmp"
mv "$OUT/$CT.sig.tsv.gz.tmp" "$OUT/$CT.sig.tsv.gz"
echo "rows: $(zcat "$OUT/$CT.sig.tsv.gz" | tail -n +2 | wc -l)"
echo "End: $(date)  ct=$CT"
