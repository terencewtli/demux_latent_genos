# Ideas: what to do with latent genotypes, and how to test it

Brainstorm, 2026-10-02. This repo is now the umbrella for latent-genotype and demultiplexing-improvement ideas.
Nothing here has been run. Status of the actual analyses: `md/PROGRESS.md`; findings: `md/RESULTS.md`.

## Where we are

- souporcell latent genotypes are noisy (naive dosage r² ~0.69). GLIMPSE2 on ambient-aware genotype likelihoods
  lifts this to **0.92–0.94 at typed sites** and **0.56 (GEX) / 0.74 (ATAC) at untyped common sites** (A04g, A04l).
- souporcell already assigns singlets almost perfectly in our 8- and 16-donor simulations. The only measurable
  demux gap is **doublet detection in 16-donor GEX pools** (0.80–0.84 vs 0.95 with true genotypes; A05b).
- So "better demultiplexing" has little headroom *in these simulations*. The open question is **utility**: are
  imputed latent genotypes good enough for genetics?

---

## Part 1. Critique of the proposed directions

### 1a. "Compare cell-type eQTL power with imputed vs true genotypes": the right question, but the wrong simulator

- **ambisim cannot answer it as built.** It samples reads from a reference dataset, so expression does not depend
  on donor genotype (already noted in ANALYSIS_PLAN, real-data section). The donor set (~16 fibroblast-derived
  donors, arbitrary) is also far too small: cis-eQTL mapping needs ~100+ donors to have power worth comparing.
- **Much of the answer is analytic, and we already have the inputs.** With non-differential genotype error, the
  association test's non-centrality scales with the dosage r² between imputed and true genotype. So the effective
  sample size is ≈ n × r² at the eQTL SNP. At our r²:
  - typed sites (0.92–0.94): ~6–8% loss of effective n;
  - untyped GEX sites (~0.56): ~45% loss.
  This predicts most of the power comparison without any simulation. A simulation is only worth building if it
  tests what the analytic argument assumes away.
- **What the analytic argument assumes away is the interesting part: the errors are not independent of expression.**
  Latent genotypes are inferred *from the same RNA reads* that measure expression.
  - **Depth follows expression.** Genotypes at SNPs in highly expressed genes are well called; genotypes for
    lowly expressed genes are mostly imputed. Error therefore varies with the phenotype being mapped.
  - **ASE biases genotype calls.** In a donor heterozygous for a strong cis-eQTL, the allele with higher
    expression dominates the reads at exonic SNPs in LD. Those SNPs drift toward being called homozygous for that
    allele, exactly at the genes with the strongest effects. Direction of the bias on effect sizes and false
    positives is not obvious; it could deflate (hets look hom, effect diluted) or create spurious structure.
  - **Ambient RNA biases toward the pool's average expression-weighted allele frequency.** That is again
    correlated with expression.
  - This "circularity" (genotype and phenotype from the same molecules) is, as far as I know, untested for
    pooled sc-eQTL designs. The bulk analogue (calling genotypes from RNA-seq for eQTL) is used, but its bias at
    ASE genes is rarely quantified. **I'd make this the central scientific question**; it also makes the
    simulator worth building, because only a genotype-aware expression simulator can test it.
- **Imputation helps exactly here.** The lead eQTL SNP is often not exonic. An imputed lead SNP gets its dosage
  from LD with many SNPs (including non-expressed ones in ATAC), so it should be less biased than a directly
  observed exonic SNP. That is a testable prediction: bias at observed exonic SNPs > bias at imputed non-exonic
  SNPs.

### 1b. "Is low-pass WGS enough to assign latent clusters to donors?": almost certainly yes, so ask a harder question

- **Assignment is easy.** Matching a cluster's latent genotype to the right donor needs a few hundred informative
  SNPs. At 0.1x, ~10% of ~5M common SNPs carry one read, i.e. hundreds of thousands of sites. The answer at 1x
  is a foregone yes. The publishable result is **the floor**: does 0.01x work, and does 0.001x? That is a
  cheap sweep: thin real high-coverage BAMs (1000G 30x), no read simulation.
- **For eQTL mapping, donor identity is not needed at all.** If each donor appears in one pool, the latent +
  imputed genotype *is* that donor's genotype for association purposes. External genotypes only buy:
  1. linking clusters to metadata (sex, age, disease, other assays);
  2. tracking the same donor across pools or time points;
  3. catching sample swaps and contamination;
  4. **better genotypes**. This is the one worth testing, because low-pass reads are unbiased by expression and
     cover the whole genome.
- **So the stronger version of the idea is joint imputation.** Feed GLIMPSE2 *both* the cluster's single-cell
  allele counts and the donor's low-pass WGS reads as one set of genotype likelihoods. Predictions:
  - beats either source alone;
  - removes most of the ASE-driven bias of 1a, because WGS reads are allele-balanced;
  - the gain over latent-only is largest at untyped and lowly expressed sites, which is where latent genotypes are
    weakest (GEX untyped 0.56).
- **Cost framing needs care.** "$30 per donor" is mostly library prep, not sequencing. At ≤ 0.5x, sequencing is a
  few dollars, and tagmentation-based low-input kits are ~$10–30 per sample. A genotyping array is ~$30–60. So the
  honest comparison is low-pass vs array vs nothing:
  - low-pass wins over arrays on non-European ancestry and on imputation of rarer variants (Li et al. 2021
    *Genome Res*);
  - **"nothing" (latent + imputation only) is the real competitor** and should be the baseline in every figure.
  - The message to experimentalists becomes: "for eQTL you may need nothing; for identity you need very little;
    here is the depth floor and the price".

### 1c. Scaling up: 16 arbitrary fibroblast donors is not a genetics study

Agreed. Options, cheapest first:
1. **Real data** (best): a large pooled sc dataset with genotypes. See Part 3.
2. **Count-level simulation**: skip reads and alignment entirely (Part 2, tier 1). This scales to 1,000 donors.
3. **Read-level ambisim** with expression made genotype-dependent: needed only to validate tier 1 on a handful of
   pools. Too expensive (cellranger per pool; disks are ~99% full) to run at eQTL scale.

### 1d. Is demultiplexing really "done"?

Not quite, and I would keep one demux thread alive:
- Our simulations use a fixed 20% ambient and well-separated EUR donors. Li 2025 Fig. S9 shows souporcell
  collapsing at low coverage past a handful of donors. So the regime where latent → impute → genotype-based demux
  could matter is **low coverage × many donors**, which we haven't simulated. A small grid (coverage × n donors) is
  cheap with the count-level simulator.
- The n16 doublet gap is real and measurable now (A05 already set up).

---

## Part 2. A simulation design that would actually test the genetics (three tiers)

**Tier 0: error injection, no reads (days).**
- Take real genotypes for many donors (1000G 3,202 or HGDP).
- Corrupt them with the *measured* error model from our runs: depth-stratified confusion matrices from A03b /
  A04g, as a function of per-site depth.
- Simulate pseudobulk cell-type expression from the true genotypes (linear model with known cis effects, or
  splatPop for per-cell counts).
- Map eQTL with true vs corrupted-then-imputed genotypes.
- Answers: power vs n, by site class.
- Cannot answer: circularity (1a), because the errors here are independent of expression. Use it as the null
  against which tier 1's bias is measured.

**Tier 1: count-level joint simulator (the main investment, ~2–3 weeks).**
- Per donor: genotype (1000G / HGDP) → per-gene expression with cis-eQTL effects, cell-type-specific (splatPop
  already does this).
- For each heterozygous exonic SNP, split that gene's reads between alleles by the ASE implied by the eQTL
  genotype at the causal variant (via haplotype phase).
- Add ambient (pool-average allele counts × expression) and doublets.
- Output per-cell **ref/alt count matrices**. souporcell's clustering and vireo both run from count matrices, so
  no BAMs and no alignment are needed. This makes 1,000 donors × 10–16 per pool affordable.
- Then: cluster → latent GLs → GLIMPSE2 → eQTL. Compare against true genotypes. Bonus: run Tier 0 on the same
  donors to separate random-error loss from circularity bias.
- Validate the count simulator against ambisim on 2–3 pools, i.e. the same donors and depths should give similar
  latent-genotype error curves.

**Tier 2: read-level (only if a reviewer insists).** Patch ambisim so per-donor expression and allele ratios
come from the Tier-1 model. Run a few pools to show tier 1 does not miss alignment or mapping-bias effects
(reference bias at het SNPs is a known ASE artefact that only read-level simulation captures).

**Low-pass arm (runs on any tier).** Thin real 30x BAMs of the same donors to 0.001–1x (`samtools view -s`).
Run GLIMPSE2 on low-pass alone, latent alone, and both jointly. Measure:
- cluster → donor assignment accuracy;
- dosage r²;
- eQTL power and bias.
This needs no read simulation because 1000G 30x BAMs exist. Use HGDP donors if the panel needs to be leakage-free.

---

## Part 3. Real data: what exists (to verify before planning around it)

| dataset | why | access (verify) |
|---|---|---|
| **OneK1K** (Yazar 2022 *Science*; ~982 donors, 75 pools of ~12–14, PBMC 10x) | The ideal test: pooled design, large n, published cell-type eQTLs to compare against, array genotypes as truth. Run souporcell per pool → impute → re-map eQTL → compare with published effects. | Processed expression is open (CELLxGENE / GEO GSE196830); tools (cellink) fetch per-donor genotypes. **Raw reads / BAMs: unclear, possibly controlled.** That is the first thing to check, because without reads there are no latent genotypes. |
| HipSci-derived pooled iPSC 10x (e.g. Cuomo 2020; Jerber 2021, ~200 lines in pools) | Pooled, differentiation time courses, sc-eQTL published; HipSci has an open-access genotype subset | Raw-data access mixed (ENA / EGA); check which lines are open |
| User's 4-donor fibroblast → iPSC (open) | ASE validation (already chosen) | Open |
| User's 51-line pools | Right scale | Protected; excluded by decision 2026-09-30 |
| **Open-access 10x multiome QTL cohort** (wanted, 2026-10-08) | A06 Level B: at simulated coverage, ATAC / multiome imputed genotypes keep most eGene and fine-mapping power (h² 0.10, untyped causal: eGene 0.72 vs GEX 0.37; causal in CS 0.44 vs 0.18). A real pooled multiome cohort with genotypes would test this. | None found yet; the only open multiome in hand is the user's 4-donor fibroblast → iPSC set (too small for QTLs). Search pending. |

If OneK1K raw reads are obtainable, it beats every simulation for the eQTL question. Tier 1 would then mainly
provide the mechanism (circularity) and the low-pass / design sweeps.

---

## Part 4. Other ideas for the umbrella

Ranked by my sense of value per effort.

1. **Depth-aware hybrid genotype** (cheap, already motivated): take the observed call where cluster depth > 10 and
   the GLIMPSE2 posterior elsewhere. A04l shows GLIMPSE2 over-smooths high-depth sites. A few hours.
2. **GEX + ATAC latent-genotype merge** for multiome: same cells, complementary sites (exons vs peaks). Should lift
   untyped GEX accuracy toward ATAC's. Cheap with existing runs.
3. **Genotype QC from latent genotypes, as a tool:**
   - sex check (chrX het rate / chrY);
   - ancestry PCs;
   - kinship between clusters, which catches related donors or a donor split across two clusters;
   - cross-pool identity: the same donor in two pools, or the same donor across time points.
   This is practical, needs no truth, and experimentalists would use it. It also guards against the failure mode
   where souporcell merges two donors: the merged "genotype" has excess heterozygosity and impossible kinship.
4. **Demux → impute → re-demux loop** aimed at the n16 doublet gap and at low-coverage × many-donor pools.
   (Step 3 is already set up.)
5. **Cross-pool designs without any genotyping:** put each donor in 2 pools with different partners. Latent
   genotypes can be matched across pools (identity from consistency, not truth), which doubles cells per donor
   and enables batch-effect correction. Connects to the pool_design repo.
6. **mtDNA variants as an orthogonal signal** in ATAC: mitochondrial heteroplasmy separates donors (and clones)
   independently of nuclear latent genotypes. It is useful as a tie-breaker for merged clusters.
7. **Privacy risk assessment.** Public scRNA BAMs plus our pipeline give genome-wide imputed genotypes of the
   donors (r² ~0.9 at typed sites). That is a re-identification risk beyond the eQTL-linking attacks discussed in
   the literature (Cho 2024 *Patterns*). It is worth a short note, carefully framed as risk quantification. **Do
   not** run it on real public data without an ethics conversation first.
8. **Monopogen comparison** (already planned): its LD refinement is the closest published method. The fair test is
   Monopogen on per-cluster pseudobulk vs souporcell latent + GLIMPSE2.

9. **ASE as an imputation check (idea only, 2026-10-08; user: not implementing, too many problems).** Donors imputed het at a known eQTL lead should show phase-consistent allelic imbalance at exonic hets of that gene in donor pseudobulk; compare imputed + GLIMPSE2-phased vs array genotypes. Problems: 3′ reads see few exonic hets, per-cell dropout and bursting, phasing error between lead and exonic SNP, and circularity when the same reads call the het. Note that strongly monoallelic sites are exactly where naive souporcell calls read homozygous.

## Suggested order

0. **NEXT (decided 2026-10-02): A06 fine-mapping / colocalization go/no-go test**, `docs/A06_FINEMAP_PLAN.md`.
   It is a focused version of Tier 0 that uses real imputation error (emulated latent genotype likelihoods +
   GLIMPSE2) and susieR / coloc. It decides whether there is a story before any simulator is built.
1. Check OneK1K raw-read access (an hour). It changes everything downstream.
2. Tier 0 plus the analytic n × r² prediction (days). This gives the baseline story (largely covered by A06).
3. Low-pass floor sweep on thinned 1000G 30x BAMs, plus joint GLIMPSE2 (a week). This is the experimentalist
   message.
4. Tier 1 count-level simulator with ASE coupling (2–3 weeks). This is the circularity question, the main science.
5. Cheap side items whenever there is idle time: hybrid genotype, GEX + ATAC merge, QC tool.

## Open questions for the user

- Is the target a methods paper ("you don't need genotyping") or a resource/tool (QC + joint imputation)? The
  first needs real-data eQTL (OneK1K-scale). The second can stand on simulations.
- How much simulator engineering is acceptable? Tier 1 is the real cost. Tier 0 + the low-pass sweep could be a
  short paper on its own.
- For the circularity question, is any real dataset with both pooled scRNA and independent WGS large enough? If
  not, Tier 1 is the only route.
