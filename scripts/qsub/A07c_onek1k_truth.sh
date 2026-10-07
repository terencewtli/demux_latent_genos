#!/bin/bash
#$ -N A07c_onek1k_truth
#$ -cwd
#$ -l h_data=8G,h_rt=24:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07c_onek1k_truth.$JOB_ID
#$ -j y

# OneK1K pilot truth genotypes + small metadata.
#   1. Zenodo 7619796 OneK1K.noGP.vcf.gz (12.7 GB; Minimac4 v1.0.2 output, GRCh37 contigs "1".."22",
#      1,098 samples OneK1K_<n>; coordinates are GRCh38)
#   2. truth = array-genotyped sites only (INFO/TYPED or TYPED_ONLY), biallelic SNVs, chr1-22.
#      Minimac4 conditions on the observed genotype at typed sites, so GT there ~ array call.
#      Imputed sites are NOT used as truth (their own error + shared-panel bias with GLIMPSE2).
#   3. rename contigs 1 -> chr1 (the VCF is already GRCh38 despite "1".."22" contig names; no liftover).
#   4. GEO per-pool donor lists (*GenotypeSamples.txt, IDs "682_683") and per-cell donor barcodes
#      (*Individual_Barcodes.csv) for the pilot pools; 10x v2 whitelist 737K-august-2016.
# Output: $SCR/A07/truth/onek1k.typed.b38.vcf.gz ; latent_genos/reference/{onek1k_geo,10x}/
# Submit with (from latent_genos/): qsub scripts/qsub/A07c_onek1k_truth.sh

source ~/.bashrc
set -euo pipefail
module load bcftools/1.11 htslib/1.12

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos/A07/truth
REF=/u/project/cluo/terencew/reference/hg38_igvf
GEO=$PROJ/reference/onek1k_geo
WLDIR=$PROJ/reference/10x
N=${NSLOTS:-1}
mkdir -p "$SCR" "$GEO" "$WLDIR"
echo "Start: $(date)  host=$(hostname -s)"

### 4 (small, first). whitelist + GEO metadata
if [ ! -s "$WLDIR/737K-august-2016.txt" ]; then
    curl -fsL https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-august-2016.txt.gz | gunzip > "$WLDIR/737K-august-2016.txt.part"
    mv "$WLDIR/737K-august-2016.txt.part" "$WLDIR/737K-august-2016.txt"
fi
echo "whitelist: $(wc -l < "$WLDIR/737K-august-2016.txt") barcodes"
while read -r POOL GSM NDON; do
    base=https://ftp.ncbi.nlm.nih.gov/geo/samples/${GSM:0:7}nnn/$GSM/suppl
    for f in $(curl -fsL "$base/" | grep -oE "${GSM}_[^\"]+(GenotypeSamples\.txt|Individual_Barcodes\.csv)\.gz" | sort -u); do
        [ -s "$GEO/$f" ] || curl -fsL -o "$GEO/$f" "$base/$f"
    done
    echo "pool $POOL ($GSM): $(zcat "$GEO/${GSM}"_*GenotypeSamples.txt.gz | grep -c .) donors listed (expected $NDON), $(zcat "$GEO/${GSM}"_*Individual_Barcodes.csv.gz | tail -n +2 | wc -l) barcodes"
done < "$PROJ/txt/onek1k_pilot_pools.txt"

### 1. download
VCF=$SCR/OneK1K.noGP.vcf.gz
if [ ! -s "$VCF.csi" ]; then
    time curl -fsL -C - -o "$VCF.part" "https://zenodo.org/records/7619796/files/OneK1K.noGP.vcf.gz?download=1"
    mv "$VCF.part" "$VCF"
    curl -fsL -o "$VCF.csi" "https://zenodo.org/records/7619796/files/OneK1K.noGP.vcf.gz.csi?download=1"
fi

### 2-3. typed SNVs, contigs renamed 1 -> chr1. NO liftover: despite the "1".."22" contig names the
### VCF is already GRCh38 (REF = hg38 base at 400/400 sampled sites, hg19 ~25%; checked 2026-10-06).
### The first run (job 15073119/15073954) lifted hg38 -> "hg38" again with Picard; its outputs
### onek1k.typed.hg38.vcf.gz / lifted*.vcf.gz are WRONG and unused (user to delete).
OUT=$SCR/onek1k.typed.b38.vcf.gz
if [ ! -s "$OUT.tbi" ]; then
    seq 1 22 | awk '{print $1"\tchr"$1}' > "$SCR/rename.txt"
    time bcftools view --threads "$N" -r "$(seq -s, 1 22)" -i 'INFO/TYPED=1 || INFO/TYPED_ONLY=1' \
        -m2 -M2 -v snps "$VCF" -Ou \
        | bcftools annotate --threads "$N" --rename-chrs "$SCR/rename.txt" -x INFO,^FORMAT/GT -Oz -o "$OUT.part.vcf.gz"
    mv "$OUT.part.vcf.gz" "$OUT"
    bcftools index -t -f "$OUT"
fi
### sanity: REF must equal the hg38 base
module load samtools 2>/dev/null || true
CHECK=$(bcftools query -f '%CHROM\t%POS\t%REF\n' "$OUT" | awk 'NR % 1000 == 1' | head -500 \
    | while read -r c p r; do b=$(samtools faidx "$REF/GRCh38.fa" "$c:$p-$p" | tail -1 | tr a-z A-Z); [ "$b" = "$r" ] && echo 1 || echo 0; done | awk '{s+=$1} END {print s "/" NR}')
echo "$(date): typed GRCh38 SNVs $(bcftools index -n "$OUT"), samples $(bcftools query -l "$OUT" | wc -l), REF==hg38 base at $CHECK sampled sites"
echo "End: $(date)"
