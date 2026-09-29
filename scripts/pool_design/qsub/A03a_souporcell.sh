#!/bin/bash
#$ -N A03a_souporcell
#$ -cwd
#$ -l h_data=2G,h_rt=24:00:00
#$ -pe shared 10
#$ -t 1-12:1
#$ -tc 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A03a_souporcell.$JOB_ID.$TASK_ID
#$ -j y

# Genotype-free souporcell on one (tree, pool, modality) row of
# txt/latent_genos_souporcell_tasks.txt, for the latent_genos project
# (github/demux_latent_genos): the cluster_genotypes.vcf it writes is the
# "latent genotype" we score against the true 1000G genotypes.
#
# Settings mirror the Li 2025 benchmark templates
# (demux_benchmark/template_demux/souporcell_nogenos_{gex,atac}.sh):
# --skip_remap with the 1000G AF>5% common-variant list, 200 restarts, --no_umi
# for ATAC. Two changes: the fasta is the cellranger-arc reference the BAMs were
# aligned to, and threads = NSLOTS (the template asked for -t 50 on 8 slots).
#
# Row format: <tree> <pool> <k> <gex|atac>
# Output: <tree>/<pool>/demux/souporcell/<modality>/<sample>/
# Submit with: qsub scripts/latent_genos/qsub/A03a_souporcell.sh

source ~/.bashrc

conda activate souporcell

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
SOUPORCELL=/u/project/cluo/terencew/programs/souporcell/souporcell_pipeline.py
FASTA=/u/project/cluo/terencew/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa
COMMON=/u/project/cluo/terencew/reference/cellsnp_db/hg38.AF5e2.possorted.chr_prefix.vcf

LIST=$PROJDIR/txt/latent_genos_souporcell_tasks.txt
# ID=1
ID=$SGE_TASK_ID
read -r TREE POOL K MOD <<< "$(sed -n "${ID}p" "$LIST")"

CR_OUTS=$PROJDIR/$TREE/$POOL/cr_arc/$SAMPLE/outs
BAM=$CR_OUTS/${MOD}_possorted_bam.bam
BARCODES=$CR_OUTS/filtered_feature_bc_matrix/barcodes.tsv.gz
OUT=$PROJDIR/$TREE/$POOL/demux/souporcell/$MOD/$SAMPLE

echo "Start: $(date)  host=$(hostname -s)  tree=$TREE pool=$POOL k=$K mod=$MOD slots=${NSLOTS:-1}"

if [ ! -f "$BAM" ] || [ ! -f "$BARCODES" ]; then
    echo "$(date): cr_arc output not found for $TREE/$POOL ($MOD), skipping"
    echo "  expected: $BAM"
    echo "  expected: $BARCODES"
    exit 0
fi

if [ -s "$OUT/clusters.tsv" ] && [ -s "$OUT/cluster_genotypes.vcf" ]; then
    echo "$(date): souporcell output exists, skipping ($OUT)"
    exit 0
fi

mkdir -p "$OUT"

EXTRA=()
if [ "$MOD" = "atac" ]; then
    EXTRA=(--no_umi True)
fi

time $SOUPORCELL -i "$BAM" -b "$BARCODES" -f "$FASTA" \
    -t "${NSLOTS:-1}" -o "$OUT" -k "$K" \
    --common_variants "$COMMON" \
    --restarts 200 --skip_remap True "${EXTRA[@]}"

echo "End: $(date)  tree=$TREE pool=$POOL mod=$MOD"
