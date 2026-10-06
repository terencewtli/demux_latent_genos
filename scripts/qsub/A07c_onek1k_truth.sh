#!/bin/bash
#$ -N A07c_onek1k_truth
#$ -cwd
#$ -l h_data=8G,h_rt=24:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07c_onek1k_truth.$JOB_ID
#$ -j y

# OneK1K pilot truth genotypes + small metadata.
#   1. Zenodo 7619796 OneK1K.noGP.vcf.gz (12.7 GB; Minimac4 v1.0.2 output, GRCh37 contigs "1".."22",
#      1,098 samples OneK1K_<n>)
#   2. truth = array-genotyped sites only (INFO/TYPED or TYPED_ONLY), biallelic SNVs, chr1-22.
#      Minimac4 conditions on the observed genotype at typed sites, so GT there ~ array call.
#      Imputed sites are NOT used as truth (their own error + shared-panel bias with GLIMPSE2).
#   3. rename contigs 1 -> chr1, Picard LiftoverVcf to GRCh38 (hg19ToHg38 chain, hg38_igvf GRCh38.fa,
#      RECOVER_SWAPPED_REF_ALT), keep chr1-22, index.
#   4. GEO per-pool donor lists (*GenotypeSamples.txt, IDs "682_683") and per-cell donor barcodes
#      (*Individual_Barcodes.csv) for the pilot pools; 10x v2 whitelist 737K-august-2016.
# Output: $SCR/A07/truth/onek1k.typed.hg38.vcf.gz ; latent_genos/reference/{onek1k_geo,10x}/
# Submit with (from latent_genos/): qsub scripts/qsub/A07c_onek1k_truth.sh

source ~/.bashrc
set -euo pipefail
module load bcftools/1.11 htslib/1.12 picard_tools/2.25.0

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

### 2-3. typed SNVs -> hg38
OUT=$SCR/onek1k.typed.hg38.vcf.gz
if [ ! -s "$OUT.tbi" ]; then
    seq 1 22 | awk '{print $1"\tchr"$1}' > "$SCR/rename.txt"
    time bcftools view --threads "$N" -r "$(seq -s, 1 22)" -i 'INFO/TYPED=1 || INFO/TYPED_ONLY=1' \
        -m2 -M2 -v snps "$VCF" -Ou \
        | bcftools annotate --threads "$N" --rename-chrs "$SCR/rename.txt" -x INFO,^FORMAT/GT -Oz -o "$SCR/onek1k.typed.hg19.vcf.gz"
    bcftools index -t -f "$SCR/onek1k.typed.hg19.vcf.gz"
    echo "$(date): typed hg19 SNVs $(bcftools index -n "$SCR/onek1k.typed.hg19.vcf.gz")"
    time java -Xmx24g -jar "$PICARD" LiftoverVcf I="$SCR/onek1k.typed.hg19.vcf.gz" O="$SCR/lifted.vcf.gz" \
        CHAIN="$REF/hg19ToHg38.over.chain.gz" REJECT="$SCR/lifted.reject.vcf.gz" R="$REF/GRCh38.fa" \
        RECOVER_SWAPPED_REF_ALT=true WARN_ON_MISSING_CONTIG=true MAX_RECORDS_IN_RAM=200000 TMP_DIR="${TMPDIR:-/tmp}"
    time bcftools view --threads "$N" -t "$(seq -s, 1 22 | sed 's/[0-9]*/chr&/g')" -m2 -M2 -v snps "$SCR/lifted.vcf.gz" -Ou \
        | bcftools sort -m 4G -T "${TMPDIR:-/tmp}" -Oz -o "$OUT.part"
    mv "$OUT.part" "$OUT"
    bcftools index -t -f "$OUT"
fi
echo "$(date): typed hg38 SNVs $(bcftools index -n "$OUT"), rejected $(zcat "$SCR/lifted.reject.vcf.gz" | grep -vc '^#'), samples $(bcftools query -l "$OUT" | wc -l)"
echo "End: $(date)"
