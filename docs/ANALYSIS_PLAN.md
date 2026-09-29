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
- **Panel choice.** This is an open decision.
  - (a) The TOPMed server itself. This is the most faithful option. It means
    uploading 1000G genotypes, which are public, so that's acceptable. Check
    whether any 1000G samples are in the TOPMed panel; if they are, that would
    leak truth into the imputation.
  - (b) A local Minimac4 run with the 1000G 30x panel, **leaving out every donor
    in the pool**. Without the leave-out, truth leaks into the panel. This option
    is fully reproducible and scriptable.
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
