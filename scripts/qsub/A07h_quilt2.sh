#!/bin/bash
#$ -N A07h_quilt2
#$ -cwd
#$ -l h_data=8G,h_rt=12:00:00
#$ -pe shared 8
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07h_quilt2.$JOB_ID.$TASK_ID
#$ -j y

# A07h arm (iii): QUILT2 read-aware imputation on the per-group BAMs (A07h_split_bams.py), the scTAPAS method
# (Kockelbergh 2026: method diploid, nGen 100, buffer 500 kb, 5 Mb chunks). Panel: 1000G 30x, the same source as
# GLIMPSE2 (A07e); PANEL_SUBSET=EUR restricts it to the 633 EUR samples if memory requires (scTAPAS used 160 GB per
# chr6 chunk with the 2,504-sample phase 3 panel).
# Three stages, submitted separately (STAGE env var, -t set at submission):
#   STAGE=chunks   (no array)  writes $QREF/<chr>/chunks.txt (5 Mb windows) and the QUILT-format genetic map
#   STAGE=prep     -t 1-<n_chunks over chr6+chr22>      QUILT2_prepare_reference.R per chunk (shared by all pools)
#   STAGE=impute   -t 1-<10 x n_chunks>                  QUILT2.R per (pool x mode) x chunk
#   STAGE=concat   -t 1-20                               concat + index per pool x mode x chrom
# Output: results/A07_onek1k/impute_quilt/<mode>/pool<P>/chrN.imputed.bcf (FORMAT DS, as the other arms)
# Order: chunks -> prep -> (A07h_split_bams) -> impute -> concat -> A07h_score.

source ~/.bashrc
conda activate quilt
set -euo pipefail

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
SRC=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom
PED=/u/project/cluo/terencew/demux_benchmark/pool_design/csv/1000G/meta/20130606_g1k_3202_samples_ped_population.txt
QREF=$SCR/quilt_ref${PANEL_SUBSET:+_$PANEL_SUBSET}
CHROMS=(chr6 chr22)
CHUNK=5000000
BUFFER=500000
N=${NSLOTS:-1}
STAGE=${STAGE:?set STAGE=chunks|prep|impute|concat}
ID=${SGE_TASK_ID:-1}
[ "$ID" = "undefined" ] && ID=1
echo "Start: $(date)  host=$(hostname -s)  stage=$STAGE id=$ID panel=${PANEL_SUBSET:-all} slots=$N"

all_chunks() { for C in "${CHROMS[@]}"; do sed "s/^/$C /" "$QREF/$C/chunks.txt"; done; }

case $STAGE in
chunks)
    for C in "${CHROMS[@]}"; do
        mkdir -p "$QREF/$C"
        PANEL=$QREF/$C/panel.vcf.gz
        if [ ! -s "$PANEL.tbi" ]; then
            SAMP=()
            if [ -n "${PANEL_SUBSET:-}" ]; then
                awk -v p="$PANEL_SUBSET" 'NR>1 && ($6==p || $7==p) {print $2}' "$PED" > "$QREF/samples.txt"
                SAMP=(-S "$QREF/samples.txt")
            fi
            time bcftools view "${SAMP[@]}" -m2 -M2 -v snps -c1 "$SRC/1000G.$C.vcf.gz" -Oz -o "$PANEL"
            tabix -f "$PANEL"
        fi
        # QUILT map: position COMBINED_rate(cM/Mb) Genetic_Map(cM), from the GLIMPSE b38 map (pos chr cM)
        zcat "$SCR/glimpse_maps/b38/$C.b38.gmap.gz" | awk 'NR>1 {if (n) {r=(p>pp)?($3-pc)/(($1-pp)/1e6):0; print pp, r, pc}
            pp=$1; pc=$3; p=$1; n=1} END {print pp, 0, pc}' | awk 'BEGIN{print "position COMBINED_rate.cM.Mb. Genetic_Map.cM."} {print}' \
            > "$QREF/$C/map.txt"
        FIRST=$(bcftools query -f '%POS\n' "$PANEL" | head -1)
        LAST=$(bcftools query -f '%POS\n' "$PANEL" | tail -1)
        awk -v a="$FIRST" -v b="$LAST" -v w="$CHUNK" 'BEGIN{s=int(a/w)*w+1; while (s<=b) {e=s+w-1; print s, e; s=e+1}}' \
            > "$QREF/$C/chunks.txt"
        echo "$C: $(wc -l < "$QREF/$C/chunks.txt") chunks, panel $(bcftools query -l "$PANEL" | wc -l) samples"
    done
    all_chunks | wc -l | awk '{print "total chunks (prep -t 1-" $1 "; impute -t 1-" 10*$1 ")"}'
    ;;
prep)
    read -r C S E <<< "$(all_chunks | sed -n "${ID}p")"
    OUT=$QREF/$C/prep_${S}_${E}
    if [ -e "$OUT/done" ]; then echo "exists"; exit 0; fi
    mkdir -p "$OUT"
    time QUILT2_prepare_reference.R --outputdir="$OUT" --chr="$C" --regionStart="$S" --regionEnd="$E" \
        --buffer="$BUFFER" --nGen=100 --reference_vcf_file="$QREF/$C/panel.vcf.gz" --genetic_map_file="$QREF/$C/map.txt"
    touch "$OUT/done"
    ;;
impute)
    NC=$(all_chunks | wc -l)
    G=$(( (ID - 1) / NC )); K=$(( (ID - 1) % NC + 1 ))
    ROW=$(( G / 2 + 1 )); MODES=(soup oracle); MODE=${MODES[$(( G % 2 ))]}
    read -r POOL GSM KK <<< "$(sed -n "${ROW}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
    read -r C S E <<< "$(all_chunks | sed -n "${K}p")"
    BAMDIR=$SCR/A07/h_bams/pool$POOL/$MODE
    OUT=$SCR/A07/quilt/$MODE/pool$POOL/$C
    F=$OUT/quilt.${C}.${S}.${E}.vcf.gz
    if [ -s "$F" ]; then echo "exists: $F"; exit 0; fi
    mkdir -p "$OUT"
    REFD=$(ls "$QREF/$C/prep_${S}_${E}"/RData/QUILT_prepared_reference.*.RData)
    echo "pool=$POOL mode=$MODE $C:$S-$E groups=$(wc -l < "$BAMDIR/names.txt")"
    time QUILT2.R --prepared_reference_filename="$REFD" --bamlist="$BAMDIR/quilt_bamlist.txt" \
        --sampleNames_file="$BAMDIR/names.txt" --chr="$C" --regionStart="$S" --regionEnd="$E" --buffer="$BUFFER" \
        --method=diploid --nGen=100 --nCores="$N" --output_filename="$F"
    ;;
concat)
    I=$((ID - 1)); CI=$((I % 2)); MI=$(((I / 2) % 2)); ROW=$((I / 4 + 1))
    C=${CHROMS[$CI]}; MODES=(soup oracle); MODE=${MODES[$MI]}
    read -r POOL GSM KK <<< "$(sed -n "${ROW}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
    IN=$SCR/A07/quilt/$MODE/pool$POOL/$C
    OUT=$PROJ/results/A07_onek1k/impute_quilt/$MODE/pool$POOL
    mkdir -p "$OUT"
    LIST=$IN/files.txt
    while read -r S E; do echo "$IN/quilt.${C}.${S}.${E}.vcf.gz"; done < "$QREF/$C/chunks.txt" > "$LIST"
    while read -r f; do [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done < "$LIST"
    time bcftools concat -f "$LIST" -Ob -o "$OUT/$C.imputed.bcf"
    bcftools index -f "$OUT/$C.imputed.bcf"
    echo "$(date): $(bcftools index -n "$OUT/$C.imputed.bcf") sites"
    ;;
esac
echo "End: $(date)"
