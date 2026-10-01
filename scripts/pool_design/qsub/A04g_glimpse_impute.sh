#!/bin/bash
#$ -N A04g_glimpse_impute
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
# GLIMPSE2 static binaries need AVX2 (Illegal instruction on older nodes, 2026-10-01)
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-264:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04g_glimpse_impute.$JOB_ID.$TASK_ID
#$ -j y

# GLIMPSE2 phasing + imputation of one (souporcell run x chromosome) from
# genotype likelihoods, against the pool's leave-pool-out binary panel (A04f).
# The GL-aware arm parallel to A04c (Eagle + Minimac4 on hard calls).
#
# Task ID -> row of txt/latent_genos_souporcell_tasks.txt (12 rows) x chr1-22:
#   row = (ID-1)/22 + 1, chr = (ID-1)%22 + 1   (same scheme as A04c)
#   e.g. ambisim_final EUR_only__greedy_maxmin__rep1 gex chr20 = 240
# chr1 needs the A04f panel built from the re-downloaded chr1 (A04h) and the
# rebuilt A04a chr1 sites file.
#
# GL model (ARM, set with qsub -v ARM=...):
#   binom    (default) PL from cluster AO/RO, binomial, base error 0.01
#   ambient  same, plus souporcell's ambient RNA fraction mixed in at the pooled
#            site alt fraction (lib/souporcell_to_gl.py --ambient)
# Optional GLIMPSE2_phase settings for cost tests: qsub -v ARM=binom,GOPTS="--burnin 2 --main 5",SUFFIX=_b2m5
#   (commas inside -v values need care; for interactive tests export GOPTS/SUFFIX).
#   Output then goes to <ARM><SUFFIX>/.
#
# Steps:
#   1. chrN.target.bcf: lib/souporcell_to_gl.py -> sites whose CHROM/POS/REF/ALT
#      exactly match the 1000G 30x panel (A04a sites file, as in A04b)
#   2. GLIMPSE2_phase per A04f chunk, --impute-reference-only-variants (the target
#      has only the ~2-5% of panel sites souporcell covers; everything else is
#      imputed from the panel)
#   3. GLIMPSE2_ligate -> chrN.imputed.bcf (GT:DS:GP, every panel site)
#   4. lib/score_latent_imputation.py (clusters -> donors from A03b assign.tsv):
#      typed-site r2 vs souporcell's naive calls, untyped r2 by 1000G MAF bin.
#      Truth = the pool VCF (<tree>/vcf/, what ambisim simulated from; MAF >= 5%
#      SNVs only, chr1 ends at 190.67 Mb); panel AF/AC from the A04a sites file
#
# Output: latent_genos/results/A04g_glimpse_impute/<ARM>/<tree>__<pool>__<mod>/chrN.*
# Requires A03a, A04a (sites), A04f (panel); scoring also needs A03b.
# Submit with: qsub scripts/latent_genos/qsub/A04g_glimpse_impute.sh

source ~/.bashrc

conda activate allcools

set -euo pipefail

module load bcftools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
LIB=$PROJDIR/scripts/latent_genos/lib
PANEL=/u/project/cluo/terencew/reference/topmed/local_1000G_30x
SCR=/u/project/cluo_scratch/terencew/claude/latent_genos
GBIN=$SCR/bin
RES=/u/project/cluo/terencew/claude/project_ideas/latent_genos/results
A03B=$RES/A03b_souporcell_dosage_r2
ARM=${ARM:-binom}
GOPTS=${GOPTS:-}
SUFFIX=${SUFFIX:-}

LIST=$PROJDIR/txt/latent_genos_souporcell_tasks.txt
# ID=240
ID=$SGE_TASK_ID
ROW=$(( (ID - 1) / 22 + 1 ))
C=$(( (ID - 1) % 22 + 1 ))
CHR=chr$C
read -r TREE POOL K MOD <<< "$(sed -n "${ROW}p" "$LIST")"

TAG=${TREE}__${POOL}__${MOD}
SDIR=$PROJDIR/$TREE/$POOL/demux/souporcell/$MOD/$SAMPLE
SOUP=$SDIR/cluster_genotypes.vcf
REF=$SCR/glimpse_ref/${TREE}__${POOL}/$CHR
OUT=$RES/A04g_glimpse_impute/$ARM$SUFFIX/$TAG
TMP=$SCR/A04g/$ARM$SUFFIX.$TAG.$CHR
N=${NSLOTS:-1}

echo "Start: $(date)  host=$(hostname -s)  arm=$ARM$SUFFIX gopts=[$GOPTS] tree=$TREE pool=$POOL mod=$MOD chr=$CHR slots=$N"

if [ ! -e "$SDIR/consensus.done" ] || [ ! -s "$SOUP" ]; then
    echo "$(date): souporcell not finished ($SDIR), skipping"
    exit 0
fi
[ -e "$REF/split.done" ] || { echo "ERROR: A04f panel missing: $REF" >&2; exit 1; }
mkdir -p "$OUT" "$TMP"

IMPUTED=$OUT/$CHR.imputed.bcf
if [ -s "$IMPUTED.csi" ]; then
    echo "$(date): $IMPUTED exists, skipping imputation"
else
    ### 1. GL target, allele-matched to the panel
    AMB=()
    [ "$ARM" = "ambient" ] && AMB=(--ambient "$SDIR/ambient_rna.txt")
    time python "$LIB/souporcell_to_gl.py" "$SOUP" --chrom "$CHR" --error 0.01 ${AMB[@]+"${AMB[@]}"} \
        2> "$OUT/$CHR.convert.log" | bcftools view -Ob -o "$TMP/raw.bcf"
    cat "$OUT/$CHR.convert.log"
    bcftools index -f "$TMP/raw.bcf"
    TGT=$OUT/$CHR.target.bcf
    bcftools isec -c none -n =2 -w 1 "$TMP/raw.bcf" "$PANEL/1000G_30x.$CHR.sites.vcf.gz" -Ob -o "$TGT"
    bcftools index -f "$TGT"
    echo "$(date): target sites raw=$(bcftools index -n "$TMP/raw.bcf") matched=$(bcftools index -n "$TGT")"

    ### 2. phase + impute per chunk
    : > "$TMP/ligate.txt"
    while read -r CIDX CCHR IRG ORG REST; do
        S=${IRG#*:}
        BIN=$REF/bin/ref_${CCHR}_${S%-*}_${S#*-}.bin
        [ -s "$BIN" ] || { echo "ERROR: missing $BIN" >&2; exit 1; }
        echo "$TMP/chunk_$CIDX.bcf" >> "$TMP/ligate.txt"
        if [ -s "$TMP/chunk_$CIDX.bcf.csi" ]; then
            echo "$(date): chunk $CIDX done, reusing"
            continue
        fi
        time "$GBIN/GLIMPSE2_phase_static" --input-gl "$TGT" --reference "$BIN" \
            --impute-reference-only-variants --threads "$N" $GOPTS \
            --output "$TMP/chunk_$CIDX.bcf" --log "$TMP/chunk_$CIDX.log" > /dev/null
    done < "$REF/chunks.txt"

    ### 3. ligate, in segments. GLIMPSE2_ligate refuses when three chunks overlap one position ("Three files
    ### overlapping"), which --sequential chunking produces around centromeres / large gaps (chr1 ~125 Mb, chr22
    ### ~27.8 Mb; 2026-10-01). Start a new segment wherever chunk i-1 and chunk i+1 overlap, ligate each segment,
    ### then concatenate (ligate output is trimmed to the chunks' output regions, so segments do not overlap; phase
    ### across a segment break is arbitrary, dosages are unaffected).
    ### Ligate output still spans the chunks' INPUT (buffered) regions, so each segment is trimmed to the output
    ### regions of its first..last chunk (column 4 of chunks.txt) before the concat.
    awk -v d="$TMP" -v c="$CHR" 'BEGIN{seg=0} {s=$3; sub(/.*:/, "", s); split(s, r, "-"); st[NR]=r[1]; en[NR]=r[2];
            o=$4; sub(/.*:/, "", o); split(o, q, "-"); os[NR]=q[1]; oe[NR]=q[2]; f[NR]=d "/chunk_" $1 ".bcf"}
        END{for (i=1; i<=NR; i++) {if (i>2 && st[i] <= en[i-2]) seg++; sg[i]=seg; print f[i] > (d "/ligate_seg" sprintf("%03d", seg) ".txt")}
            for (i=1; i<=NR; i++) {if (!(sg[i] in a)) a[sg[i]]=os[i]; b[sg[i]]=oe[i]}
            for (k in a) print c ":" a[k] "-" b[k] > (d "/ligate_seg" sprintf("%03d", k) ".region")}' "$REF/chunks.txt"
    : > "$OUT/$CHR.ligate.log"
    : > "$TMP/segments.txt"
    for L in "$TMP"/ligate_seg*.txt; do
        time "$GBIN/GLIMPSE2_ligate_static" --input "$L" --output "${L%.txt}.full.bcf" --threads "$N" >> "$OUT/$CHR.ligate.log"
        bcftools index -f "${L%.txt}.full.bcf"
        bcftools view -r "$(cat "${L%.txt}.region")" "${L%.txt}.full.bcf" -Ob -o "${L%.txt}.bcf"
        bcftools index -f "${L%.txt}.bcf"
        echo "$(basename "$L" .txt) $(cat "${L%.txt}.region") $(bcftools index -n "${L%.txt}.bcf") sites" >> "$OUT/$CHR.ligate.log"
        echo "${L%.txt}.bcf" >> "$TMP/segments.txt"
    done
    echo "$(date): ligated $(wc -l < "$TMP/segments.txt") segment(s)"
    if [ "$(wc -l < "$TMP/segments.txt")" -eq 1 ]; then
        mv "$(cat "$TMP/segments.txt")" "$TMP/imputed.bcf"
    else
        bcftools concat --naive -f "$TMP/segments.txt" -Ob -o "$TMP/imputed.bcf"
    fi
    mv "$TMP/imputed.bcf" "$IMPUTED"
    bcftools index -f "$IMPUTED"
    cat "$TMP"/chunk_*.log > "$OUT/$CHR.phase.log"
    echo "$(date): imputed $(bcftools index -n "$IMPUTED") sites"
    rm -f "$TMP"/chunk_* "$TMP"/ligate_seg* "$TMP/segments.txt" "$TMP/ligate.txt" "$TMP/raw.bcf" "$TMP/raw.bcf.csi"
fi

### 4. score
ASSIGN=$A03B/$TAG/assign.tsv
if [ ! -s "$ASSIGN" ]; then
    echo "$(date): no A03b assign.tsv for $TAG, skipping scoring (run A03b)"
elif [ -s "$OUT/$CHR.score.donor.tsv" ]; then
    echo "$(date): scores exist, skipping"
else
    time python "$LIB/score_latent_imputation.py" --dose "$IMPUTED" --typed "$OUT/$CHR.target.bcf" \
        --soup "$SOUP" --assign "$ASSIGN" --truth "$PROJDIR/$TREE/vcf/$POOL.vcf.gz" \
        --panel-sites "$PANEL/1000G_30x.$CHR.sites.vcf.gz" --region "$CHR" --out "$OUT/$CHR.score"
fi

rmdir "$TMP" 2>/dev/null || true
echo "End: $(date)  arm=$ARM tree=$TREE pool=$POOL mod=$MOD chr=$CHR"
