#!/bin/bash
#$ -N A06b_level_a_finemap
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 2
# R user library has AVX2 builds (data.table crashes with 'illegal operand' elsewhere)
#$ -l arch=intel-gold*|intel-E5-2650|intel-6736p
#$ -t 1-100
#$ -tc 50
#$ -o /u/project/cluo/terencew/claude/project_ideas/latent_genos/logs/A06b_level_a_finemap.$JOB_ID.$TASK_ID
#$ -j y

# A06b: Level A fine-mapping smoke test (docs/A06_FINEMAP_PLAN.md); task = region row of <indir>/regions.tsv.
# Usage: qsub [-t ...] scripts/qsub/A06b_level_a_finemap.sh <indir> <outdir> [n_seeds=3] [level=A|B]
source /u/local/Modules/default/init/bash
module load R/4.1.0
export LD_LIBRARY_PATH=/u/local/compilers/gcc/10.2.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=${NSLOTS:-2} OPENBLAS_NUM_THREADS=${NSLOTS:-2} MKL_NUM_THREADS=${NSLOTS:-2}
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/latent_genos
INDIR=$1
OUTDIR=$2
NSEEDS=${3:-3}
LEVEL=${4:-A}
# ID=1
ID=$SGE_TASK_ID
REGION=$(awk -F'\t' -v i=$((ID + 1)) 'NR==i{print $NF}' $INDIR/regions.tsv)
echo "$(date): A06b $REGION on $(hostname -s)"
time Rscript $PROJDIR/scripts/A06b_level_a_finemap.R $INDIR $REGION $OUTDIR $NSEEDS $LEVEL
echo "$(date): done"
