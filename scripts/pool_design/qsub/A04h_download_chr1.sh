#!/bin/bash
#$ -N A04h_download_chr1
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A04h_download_chr1.$JOB_ID
#$ -j y

# Re-download the NYGC 1000G 30x phased panel for chr1 (the by_chrom copy in demux_benchmark fails gzip -t
# with a crc error; reads stop at ~190.7 Mb). Same file on EBI: 2,390,515,911 bytes. EBI gives no md5, so
# integrity = byte size + full gzip -t + tabix index. The demux_benchmark copy is left untouched.
set -euo pipefail
OUT=/u/project/cluo/terencew/claude/project_ideas/latent_genos/reference/1000G_30x
U=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
F=1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
cd $OUT
echo "$(date): download on $(hostname -s)"
time wget -q -c $U/$F
time wget -q -c $U/$F.tbi
test "$(stat -c %s $F)" = 2390515911 && echo "size ok"
time gzip -t $F && echo "gzip ok"
md5sum $F > $F.md5
cat $F.md5
module load htslib
echo "last chr1 record: $(tabix $F chr1:248000000-249000000 | tail -1 | cut -f1-2)"
echo "$(date): done"
