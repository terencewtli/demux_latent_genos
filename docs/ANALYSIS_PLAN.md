# Analysis plan

This is a working draft. Expect it to change once the first souporcell runs land.

## Step 0: souporcell on a few pools (running)

`A03a_souporcell.sh` runs on 6 pools × {GEX, ATAC} = 12 tasks
(`txt/latent_genos_souporcell_tasks.txt`):

- `ambisim_n16/`: EUR_only × {random, greedy_maxmin} × rep{1,2}, with k=16
- `ambisim_final/`: EUR_only × {random, greedy_maxmin} × rep1, with k=8. Same
  universe and strategy as the n=16 pools, but different donors. `ambisim/` (the
  n=8 pools the n16 set was matched to) no longer has its BAMs.

The settings are the Li 2025 benchmark settings: `--skip_remap` against the 1000G
AF>5% common-variant list, 200 restarts, and `--no_umi` for ATAC. The output is
`<tree>/<pool>/demux/souporcell/<mod>/<sample>/{clusters.tsv,cluster_genotypes.vcf,ambient_rna.txt}`.

## Step 1: score the latent genotypes

1. **Cluster → donor matching.** For each cluster and each true donor, compute the
   correlation of alt-allele dosage over the sites they share. Match with the
   Hungarian algorithm. Also record the second-best match margin; a merged or
   split cluster shows up as a small margin or a many-to-one match. Score cell
   assignment accuracy against `drop_data_rand.txt` at the same time, so that
   genotype quality and demux quality can be related within a pool.
2. **Genotype metrics**, computed per donor and pooled by allele-frequency bin:
   - **Dosage r²**: squared Pearson correlation of true dosage against expected
     dosage (from GP/GL if souporcell writes them, otherwise from GT). This is the
     standard imputation metric, so it compares directly before and after
     imputation.
   - **Non-reference concordance** (NRC). Plain concordance is inflated by
     hom-ref sites, so report NRC as well.
   - **Genotype error by true class** (hom-ref / het / hom-alt). Ambient
     contamination should push hom calls toward het; check whether it does.
   - **Site yield**: the number and fraction of true variant sites with any
     latent call. Imputation's second job is to fill these in.
3. **Stratify** by per-site read depth in the cluster (souporcell's AO+RO), by
   modality, by n=8 vs n=16, and by pool strategy. The depth → r² curve is the
   most useful single plot, because it connects directly to the coverage
   breakdown in Fig. S9.

## Step 2: imputation

- **Protocol.** Match the TOPMed Imputation Server: Eagle2 pre-phasing, then
  Minimac4, and report its Rsq. Chunk-level QC on the server drops chunks with low
  call rate or low overlap with the reference panel, and sparse latent genotypes
  may fail it. Record how much of the genome survives.
- **Panel (decided 2026-09-29).** Local Minimac4 on the NYGC 1000G 30x phased
  panel (3,202 samples, SHAPEIT4). It is already on disk:
  `demux_benchmark/pool_design/vcf/1000G/by_chrom/1000G.chr*.vcf.gz`, from
  `ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/`.
  Build a separate panel for each pool that leaves out the pool's donors **and
  their 1000G relatives**. The 3,202 samples include 602 trios, so dropping only
  the donors still leaks their haplotypes through parents and children; this
  matters most for the adversarial_family pools. Eagle 2.4.1 genetic map:
  `programs/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz`.
- **Scripts.** `A04a` builds one msav per chromosome from the full panel. `A04b`
  converts souporcell output: the GO field holds natural-log likelihoods in the
  order 0/0, 1/1, 0/1, which are rewritten as log10 GL in VCF order; BACKGROUND
  loci are dropped; sites must exactly match a panel site; then the two QC sets
  (`server`: call rate >= 0.9 and not monomorphic; `all`: not monomorphic). It
  also writes the leave-out lists (donors + parents, children and siblings from
  the 3,202 pedigree). `A04c` runs Eagle against a per-pool leave-out reference
  at the target sites, then minimac4 on the full msav with
  `--sample-ids-file keep.txt`.
- **TOPMed server leakage.** The r3 panel is 133,597 TOPMed WGS samples from
  NHLBI cohort studies. I could not confirm from the docs whether any 1000G
  samples are in it. Empirical check before trusting server results: submit the
  donors' *true* genotypes thinned to array density, then look at r² for variants
  that are rare in 1000G. Near-perfect recovery of donor-private variants would
  mean the donors are in the panel.
- **Genotype-likelihood-aware alternative.** Latent genotypes are really
  low-coverage allele counts, not hard calls. GLIMPSE2 or Beagle 4.1 in GL mode
  would use the counts directly and may do much better than Minimac4 on hard
  calls. Worth running side by side.
- Score the imputed genotypes with the same Step 1 metrics, split into sites that
  souporcell called and sites that imputation added.

## Step 3: close the loop

Run demuxlet (and the Alvarez 2025 method) with (i) true, (ii) raw latent, and
(iii) imputed latent genotypes. Compare singlet accuracy and doublet calls against
souporcell's own `clusters.tsv`, reusing the pool_design scoring notebooks.

**Expectation.** Imputation can only help clusters that are already separated
correctly. If souporcell merges two donors, the latent genotype is a mixture, and
no amount of imputation recovers the missing donor. So the realistic gain is in
low-UMI cells and doublets near the decision boundary, not in fixing
cluster-level failures. Step 1's matching margins tell us which regime each pool
is in.

## Other ideas

- **Merge the modalities.** The GEX and ATAC clusters come from the same cells, so
  their latent genotypes can be matched and combined: exonic sites from GEX, and
  promoter/enhancer sites from ATAC. This increases site yield before imputation.
- **Downsample** cells or reads within one pool to map genotype quality against
  depth per donor. This is cheaper than new simulations and speaks to Fig. S9.
- **Iterate.** Impute, re-demultiplex with genotypes, re-pseudobulk each donor,
  and re-call. Check whether it converges.
- **Beyond demultiplexing.** Kinship and ancestry recovery from latent genotypes,
  and whether they are good enough for cis-QTL mapping in pooled designs.
- **Monopogen.** The analysis plan is still to be decided. The natural first
  version is Monopogen on per-cluster pseudobulk BAMs, compared against
  souporcell's latent genotypes.
- **Low-pass WGS.** Simulate 0.5–1x WGS for the pool donors, impute with
  GLIMPSE2, and ask whether it is enough for genotype-based demux.
