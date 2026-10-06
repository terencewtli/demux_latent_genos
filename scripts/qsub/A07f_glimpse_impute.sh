#!/bin/bash
#$ -N A07f_glimpse_impute
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
# GLIMPSE2 static binaries need AVX2 (Illegal instruction on older nodes, 2026-10-01)
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-110:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07f_glimpse_impute.$JOB_ID.$TASK_ID
#$ -j y

# OneK1K pilot: GLIMPSE2 imputation of one (pool x chromosome) from souporcell cluster GLs, against
# the full 1000G panel (A07e). Same method as the simulation headline arm, A04g ARM=ambient with
# GOPTS="--burnin 2 --main 5" (ambient_b2m5: GEX typed r2 0.94 in simulation):
#   1. lib/souporcell_to_gl.py --ambient -> sites matching the 1000G 30x panel exactly
#   2. GLIMPSE2_phase per chunk, --impute-reference-only-variants
#   3. GLIMPSE2_ligate in segments (A04g's "Three files overlapping" workaround), concat
# Scoring is separate (A07g) because truth is the OneK1K array, not a pool VCF.
# Task ID -> row of txt/onek1k_pilot_pools.txt x chr1-22: row = (ID-1)/22 + 1, chr = (ID-1)%22 + 1
# Output: latent_genos/results/A07_onek1k/impute/pool<P>/chrN.{imputed.bcf,target.bcf,*.log}
# Requires A07d, A07e. Submit with: qsub -hold_jid A07d_souporcell,A07e_glimpse_ref_full scripts/qsub/A07f_glimpse_impute.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
LIB=/u/project/cluo/terencew/claude/project_ideas/pool_design/scripts/latent_genos/lib
PANEL=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
GBIN=$SCR/bin
GOPTS=${GOPTS:---burnin 2 --main 5}
# ID=22
ID=$SGE_TASK_ID
ROW=$(( (ID - 1) / 22 + 1 ))
C=$(( (ID - 1) % 22 + 1 ))
CHR=chr$C
read -r POOL GSM K <<< "$(sed -n "${ROW}p" "$PROJ/txt/onek1k_pilot_pools.txt")"
SDIR=$SCR/A07/souporcell/pool$POOL
SOUP=$SDIR/cluster_genotypes.vcf
REF=$SCR/glimpse_ref/onek1k_full_1000G/$CHR
OUT=$PROJ/results/A07_onek1k/impute/pool$POOL
TMP=$SCR/A07/tmp/impute.pool$POOL.$CHR
N=${NSLOTS:-1}

echo "Start: $(date)  host=$(hostname -s)  pool=$POOL k=$K chr=$CHR gopts=[$GOPTS] slots=$N"
[ -s "$SOUP" ] && [ -s "$SDIR/ambient_rna.txt" ] || { echo "ERROR: souporcell output missing ($SDIR)" >&2; exit 1; }
[ -e "$REF/split.done" ] || { echo "ERROR: A07e panel missing: $REF" >&2; exit 1; }
mkdir -p "$OUT" "$TMP"
IMPUTED=$OUT/$CHR.imputed.bcf
if [ -s "$IMPUTED.csi" ]; then echo "$(date): $IMPUTED exists, skipping"; exit 0; fi

### 1. GL target, allele-matched to the panel
time python "$LIB/souporcell_to_gl.py" "$SOUP" --chrom "$CHR" --error 0.01 --ambient "$SDIR/ambient_rna.txt" \
    2> "$OUT/$CHR.convert.log" | bcftools view -Ob -o "$TMP/raw.bcf"
cat "$OUT/$CHR.convert.log"
bcftools index -f "$TMP/raw.bcf"
TGT=$OUT/$CHR.target.bcf
bcftools isec -c none -n =2 -w 1 "$TMP/raw.bcf" "$PANEL/1000G_30x.$CHR.sites.vcf.gz" -Ob -o "$TGT"
bcftools index -f "$TGT"
echo "$(date): target sites raw=$(bcftools index -n "$TMP/raw.bcf") matched=$(bcftools index -n "$TGT")"

### 2. phase + impute per chunk
while read -r CIDX CCHR IRG ORG REST; do
    S=${IRG#*:}
    BIN=$REF/bin/ref_${CCHR}_${S%-*}_${S#*-}.bin
    [ -s "$BIN" ] || { echo "ERROR: missing $BIN" >&2; exit 1; }
    [ -s "$TMP/chunk_$CIDX.bcf.csi" ] && continue
    time "$GBIN/GLIMPSE2_phase_static" --input-gl "$TGT" --reference "$BIN" \
        --impute-reference-only-variants --threads "$N" $GOPTS \
        --output "$TMP/chunk_$CIDX.bcf" --log "$TMP/chunk_$CIDX.log" > /dev/null
done < "$REF/chunks.txt"

### 3. ligate in segments (see A04g for why), trim each to its chunks' output regions, concat
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
echo "End: $(date)  pool=$POOL chr=$CHR"
