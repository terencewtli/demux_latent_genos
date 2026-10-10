#!/bin/bash
#$ -N A09e_glimpse_impute
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
# GLIMPSE2 static binaries need AVX2 (Illegal instruction on older nodes, 2026-10-01)
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-682:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A09e_glimpse_impute.$JOB_ID.$TASK_ID
#$ -j y

# A09 IGVF multiome: GLIMPSE2 imputation of one (target x chromosome) from donor-level souporcell GLs
# (scripts/A09e_build_gl.py; clusters renamed to donors by A09d). Same GLIMPSE2 settings and full 1000G panel as
# A07f (ambient_b2m5). Targets = rows of txt/A09_targets.txt (all_{multi,gex,atac}, then ys3?_{multi,gex,atac}; ys3c gex + multi dropped, see A09e_build_gl.py).
# Task ID -> row = (ID-1)/22 + 1, chr = (ID-1)%22 + 1.
# Output: results/A09_multiome/impute/<target>/chrN.{imputed.bcf,target.bcf,*.log}; skips finished chromosomes.
# Requires A09d. Submit with (from latent_genos/): qsub [-t 1-66] scripts/qsub/A09e_glimpse_impute.sh

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

PROJ=/u/project/cluo/terencew/claude/project_ideas/latent_genos
PANEL=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
GBIN=$SCR/bin
GOPTS=${GOPTS:---burnin 2 --main 5}
# ID=22
ID=$SGE_TASK_ID
ROW=$(( (ID - 1) / 22 + 1 ))
C=$(( (ID - 1) % 22 + 1 ))
CHR=chr$C
TARGET=$(sed -n "${ROW}p" "$PROJ/txt/A09_targets.txt")
REF=$SCR/glimpse_ref/onek1k_full_1000G/$CHR
OUT=$PROJ/results/A09_multiome/impute/$TARGET
TMP=$SCR/A09/tmp/impute.$TARGET.$CHR
N=${NSLOTS:-1}

echo "Start: $(date)  host=$(hostname -s)  target=$TARGET chr=$CHR gopts=[$GOPTS] slots=$N"
[ -s "$PROJ/results/A09_multiome/match/cluster_map.tsv" ] || { echo "ERROR: A09d cluster_map.tsv missing" >&2; exit 1; }
[ -e "$REF/split.done" ] || { echo "ERROR: panel missing: $REF" >&2; exit 1; }
mkdir -p "$OUT" "$TMP"
IMPUTED=$OUT/$CHR.imputed.bcf
if [ -s "$IMPUTED.csi" ]; then echo "$(date): $IMPUTED exists, skipping"; exit 0; fi

### 1. donor-level GL target, allele-matched to the panel
time python "$PROJ/scripts/A09e_build_gl.py" --target "$TARGET" --chrom "$CHR" \
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
    \cp -f "$(cat "$TMP/segments.txt")" "$TMP/imputed.bcf"
else
    bcftools concat --naive -f "$TMP/segments.txt" -Ob -o "$TMP/imputed.bcf"
fi
mv "$TMP/imputed.bcf" "$IMPUTED"
bcftools index -f "$IMPUTED"
cat "$TMP"/chunk_*.log > "$OUT/$CHR.phase.log"
echo "$(date): imputed $(bcftools index -n "$IMPUTED") sites  (tmp files left in $TMP)"
echo "End: $(date)  target=$TARGET chr=$CHR"
