#!/bin/bash
#$ -N A07j_impute
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
# GLIMPSE2 static binaries need AVX2 (Illegal instruction on older nodes, 2026-10-01)
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-8:1
#$ -tc 8
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A07j_impute.$JOB_ID.$TASK_ID
#$ -j y

# A07j spacing vs per-site information: GLIMPSE2 on the simulated targets of scripts/A07j_targets.py (54 OneK1K
# pilot donors), same panel and options as A07f (ambient_b2m5: --burnin 2 --main 5). Steps 2-3 are A07f's verbatim.
# Task -> arm x chrom: arms perfect_soup, perfect_random, weak_soup, weak_random; chroms chr20, chr22.
# Output: <scratch>/A07j/impute/<arm>.<chr>.imputed.bcf. Requires A07j_targets.py for both chromosomes.

source ~/.bashrc
conda activate allcools
set -euo pipefail
module load bcftools

SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
GBIN=$SCR/bin
GOPTS=${GOPTS:---burnin 2 --main 5}
ARMS=(perfect_soup perfect_random weak_soup weak_random)
CHRS=(chr20 chr22)
# ID=1
ID=$SGE_TASK_ID
ARM=${ARMS[$(( (ID - 1) % 4 ))]}
CHR=${CHRS[$(( (ID - 1) / 4 ))]}
REF=$SCR/glimpse_ref/onek1k_full_1000G/$CHR
TGT=$SCR/A07j/targets/$ARM.$CHR.bcf
OUT=$SCR/A07j/impute
TMP=$SCR/A07j/tmp/$ARM.$CHR
N=${NSLOTS:-1}
echo "Start: $(date)  host=$(hostname -s)  arm=$ARM chr=$CHR gopts=[$GOPTS] slots=$N"
[ -s "$TGT.csi" ] || { echo "ERROR: target missing: $TGT" >&2; exit 1; }
mkdir -p "$OUT" "$TMP"
IMPUTED=$OUT/$ARM.$CHR.imputed.bcf
if [ -s "$IMPUTED.csi" ]; then echo "$(date): $IMPUTED exists, skipping"; exit 0; fi

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
: > "$OUT/$ARM.$CHR.ligate.log"
: > "$TMP/segments.txt"
for L in "$TMP"/ligate_seg*.txt; do
    "$GBIN/GLIMPSE2_ligate_static" --input "$L" --output "${L%.txt}.full.bcf" --threads "$N" >> "$OUT/$ARM.$CHR.ligate.log"
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
cat "$TMP"/chunk_*.log > "$OUT/$ARM.$CHR.phase.log"
echo "$(date): imputed $(bcftools index -n "$IMPUTED") sites  (tmp files left in $TMP)"
echo "End: $(date)  arm=$ARM chr=$CHR"
