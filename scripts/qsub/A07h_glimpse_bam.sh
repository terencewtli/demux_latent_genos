#!/bin/bash
#$ -N A07h_glimpse_bam
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
# GLIMPSE2 static binaries need AVX2 (Illegal instruction on older nodes, 2026-10-01)
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-20:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07h_glimpse_bam.$JOB_ID.$TASK_ID
#$ -j y

# A07h arm (ii): GLIMPSE2 directly on per-group BAMs (A07h_split_bams.py), so genotype likelihoods are computed
# at EVERY 1000G panel site covered by reads, not only souporcell's site list (A07f, arm i). Same panel (A07e),
# same burn-in / main iterations and the same segmented ligation as A07f.
# Task -> (pool row, mode, chrom): 5 pools x {soup, oracle} x {chr6, chr22}
#   ID-1 = ((row-1) * 2 + mode) * 2 + chrom
# Output: latent_genos/results/A07_onek1k/impute_bam/<mode>/pool<P>/chrN.imputed.bcf
# Requires A07h_split (BAMs) and A07e. Overrides for tests: BAMDIR, OUTDIR, ONLY_CHUNKS (e.g. 0,1).

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
GBIN=$SCR/bin
FASTA=/u/project/cluo/terencew/reference/hg38_igvf/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
GOPTS=${GOPTS:---burnin 2 --main 5}
# ID=1
ID=$SGE_TASK_ID
I=$((ID - 1))
CI=$((I % 2)); MI=$(((I / 2) % 2)); ROW=$((I / 4 + 1))
CHRS=(chr6 chr22); MODES=(soup oracle)
CHR=${CHRS[$CI]}; MODE=${MODES[$MI]}
read -r POOL GSM K <<< "$(sed -n "${ROW}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
BAMDIR=${BAMDIR:-$SCR/A07/h_bams/pool$POOL/$MODE}
OUT=${OUTDIR:-$PROJ/results/A07_onek1k/impute_bam/$MODE/pool$POOL}
REF=$SCR/glimpse_ref/onek1k_full_1000G/$CHR
TMP=$SCR/A07/tmp/impute_bam.$MODE.pool$POOL.$CHR
N=${NSLOTS:-1}

echo "Start: $(date)  host=$(hostname -s)  pool=$POOL mode=$MODE chr=$CHR gopts=[$GOPTS] slots=$N"
[ -s "$BAMDIR/stats.tsv" ] || { echo "ERROR: BAMs missing ($BAMDIR)" >&2; exit 1; }
[ -e "$REF/split.done" ] || { echo "ERROR: A07e panel missing: $REF" >&2; exit 1; }
mkdir -p "$OUT" "$TMP"
IMPUTED=$OUT/$CHR.imputed.bcf
if [ -s "$IMPUTED.csi" ]; then echo "$(date): $IMPUTED exists, skipping"; exit 0; fi
echo "groups: $(wc -l < "$BAMDIR/bamlist.txt")"

### 1. phase + impute per chunk, GLs from the BAMs at panel sites
while read -r CIDX CCHR IRG ORG REST; do
    if [ -n "${ONLY_CHUNKS:-}" ] && ! echo ",$ONLY_CHUNKS," | grep -q ",$CIDX,"; then continue; fi
    S=${IRG#*:}
    BIN=$REF/bin/ref_${CCHR}_${S%-*}_${S#*-}.bin
    [ -s "$BIN" ] || { echo "ERROR: missing $BIN" >&2; exit 1; }
    [ -s "$TMP/chunk_$CIDX.bcf.csi" ] && continue
    time "$GBIN/GLIMPSE2_phase_static" --bam-list "$BAMDIR/bamlist.txt" --fasta "$FASTA" --reference "$BIN" \
        --threads "$N" $GOPTS --output "$TMP/chunk_$CIDX.bcf" --log "$TMP/chunk_$CIDX.log" > /dev/null
done < "$REF/chunks.txt"
if [ -n "${ONLY_CHUNKS:-}" ]; then echo "$(date): test chunks done in $TMP"; exit 0; fi

### 2. ligate in segments (as A07f / A04g), trim each to its chunks' output regions, concat
awk -v d="$TMP" -v c="$CHR" 'BEGIN{seg=0} {s=$3; sub(/.*:/, "", s); split(s, r, "-"); st[NR]=r[1]; en[NR]=r[2];
        o=$4; sub(/.*:/, "", o); split(o, q, "-"); os[NR]=q[1]; oe[NR]=q[2]; f[NR]=d "/chunk_" $1 ".bcf"}
    END{for (i=1; i<=NR; i++) {if (i>2 && st[i] <= en[i-2]) seg++; sg[i]=seg; print f[i] > (d "/ligate_seg" sprintf("%03d", seg) ".txt")}
        for (i=1; i<=NR; i++) {if (!(sg[i] in a)) a[sg[i]]=os[i]; b[sg[i]]=oe[i]}
        for (k in a) print c ":" a[k] "-" b[k] > (d "/ligate_seg" sprintf("%03d", k) ".region")}' "$REF/chunks.txt"
: > "$OUT/$CHR.ligate.log"
: > "$TMP/segments.txt"
for L in "$TMP"/ligate_seg*.txt; do
    "$GBIN/GLIMPSE2_ligate_static" --input "$L" --output "${L%.txt}.full.bcf" --threads "$N" >> "$OUT/$CHR.ligate.log"
    bcftools index -f "${L%.txt}.full.bcf"
    bcftools view -r "$(cat "${L%.txt}.region")" "${L%.txt}.full.bcf" -Ob -o "${L%.txt}.bcf"
    bcftools index -f "${L%.txt}.bcf"
    echo "${L%.txt}.bcf" >> "$TMP/segments.txt"
done
if [ "$(wc -l < "$TMP/segments.txt")" -eq 1 ]; then
    cp "$(cat "$TMP/segments.txt")" "$TMP/imputed.bcf"
else
    bcftools concat --naive -f "$TMP/segments.txt" -Ob -o "$TMP/imputed.bcf"
fi
mv "$TMP/imputed.bcf" "$IMPUTED"
bcftools index -f "$IMPUTED"
cat "$TMP"/chunk_*.log > "$OUT/$CHR.phase.log"
echo "$(date): imputed $(bcftools index -n "$IMPUTED") sites  (tmp files left in $TMP)"
echo "End: $(date)  pool=$POOL mode=$MODE chr=$CHR"
