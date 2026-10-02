# NEXT STEP (A06): does latent-genotype error damage eQTL fine-mapping and colocalization?

Planned 2026-10-02. Not run yet. Motivation and context: `docs/IDEAS.md` (Part 1a and the go/no-go discussion).

## The question

Imputed latent genotypes are accurate at sites souporcell observes (dosage r² ~0.93) and much less accurate at
sites filled in by imputation (untyped common sites: r² ~0.56 GEX, ~0.74 ATAC; A04g / A04l).
- **eGene discovery** probably survives, because a gene is detected if *any* SNP tags the causal variant, and the
  gene's own exonic SNPs are the well-genotyped ones.
- **Fine-mapping and colocalization** may not. They compare association strength across SNPs in LD, so uneven
  accuracy could tilt the signal toward well-genotyped exonic SNPs, moving lead SNPs and breaking credible sets.

This is a **cheap go/no-go test** before building any read-level or count-level eQTL simulator.

**Deliberately out of scope:** the "circularity" bias, where genotype errors depend on expression through
allele-specific expression. Here the errors are realistic in depth and LD structure but independent of the
simulated phenotype. If A06 finds damage, circularity can only add to it. If A06 finds none, circularity becomes
the remaining open question (IDEAS Part 2, tier 1).

## Design

### Donors and genotypes (no read simulation)
- **Source.** 1000G NYGC 30x phased genotypes already on disk in pool_design:
  - `demux_benchmark/pool_design/vcf/1000G/by_chrom/1000G.chr*.vcf.gz`;
  - metadata in `demux_benchmark/pool_design/csv/1000G/meta/` (633 EUR; 525 with no parents in the set).
- **Test cohort.** n = 300 unrelated EUR (pedigree-filtered), plus an n = 150 subsample for a sample-size contrast.
  This is the scale of a typical sc-eQTL study (OneK1K-like designs are larger; that is a later scale-up).
- **Imputation panel.** All other 1000G samples minus the test donors and their relatives (reuse
  `lib/leave_pool_out.py` logic). That leaves ~225 EUR plus all other ancestries.
  - This is fewer EUR haplotypes than the A04 runs had (those left out only one pool's donors), so imputation here
    will be somewhat *worse* than measured.
  - The emulated r² will be reported next to A04's to keep that visible.
- **Regions.** 100 protein-coding genes on chr16–22, each with a ±500 kb cis window (~100 Mb in total), restricted
  to biallelic SNPs with MAF ≥ 1% in the test cohort.

### Genotype arms
Each test donor gets the per-site read depth profile of a **real souporcell cluster**:
- the depth vector sampled from the 80 GEX or 80 ATAC clusters in the A03a runs;
- this preserves the real pattern of which sites have reads (expressed exons for GEX, peaks for ATAC).

Allele counts are then drawn from the donor's true genotype, with 20% ambient at the cohort allele frequency (as in
the ambisim pools). Genotype likelihoods come from `lib/souporcell_to_gl.py --ambient`, and imputation from
GLIMPSE2 against the leave-out panel, using the same settings as the best A04g arm (ambient, b2m5).

| arm | what it represents |
|---|---|
| `truth` | 1000G genotypes |
| `gex_imp` | GEX latent genotype + GLIMPSE2 (main test arm) |
| `atac_imp` | ATAC latent genotype + GLIMPSE2 |
| `multiome_imp` | GEX + ATAC allele counts summed per site, then GLIMPSE2 |
| `gex_naive` | GEX hard calls at observed sites only, no imputation (the "do nothing" floor) |
| `lowpass_0.5x` *(optional)* | 0.5x reads at every site (Poisson depth, no ambient) + GLIMPSE2. Previews the low-pass WGS question with the same machinery. |

**Level A smoke test (run first, hours).**
- Before any GLIMPSE2, inject independent per-site noise into the truth dosages, calibrated to the measured r² by
  site class (observed ~0.93; untyped GEX ~0.56 / ATAC ~0.74).
- It checks the pipeline end to end and gives an upper bound on how much an uneven accuracy alone matters.
- Real imputation error is correlated along LD blocks, which is why Level B (the arms above) is the real test.

### eQTL simulator: a donor-level pseudobulk linear model (custom, ~100 lines of R)
- **Phenotype.** y = Σ β_k G_k + ε per gene, with G the true dosages, k ∈ causal set, and ε Gaussian. The β are
  scaled to a cis-heritability h² ∈ {0.05, 0.10, 0.20}.
- **Causal configurations.**
  - 1 causal (70% of genes) or 2 causals (30%).
  - Causal-variant class drawn in three scenarios: (i) any cis SNP; (ii) only untyped / noncoding SNPs (souporcell
    never observes them; the realistic worst case, since most causal eQTL variants are noncoding); (iii) only
    observed exonic SNPs (best case).
- **Why not splatPop.** The question is about genotype error, and most sc-eQTL pipelines map donor-level pseudobulk.
  splatPop's per-cell count model only adds noise that is equivalent to lowering h², which the h² grid already
  covers. splatPop stays the choice for tier-1 work, where per-cell counts and ASE matter.
- **Replicates.** 3 seeds × 100 genes × 3 h² × 3 causal scenarios × 2 cohort sizes.

### Fine-mapping and colocalization
- **Fine-mapping.** susieR 0.14.2 on individual-level data, `susie(X, y, L = 10)`, X = dosage matrix of each arm
  (installed in `module load R/4.1.0`, user library; smoke-tested 2026-10-02 — X must be stored as double).
- **eGene call.** At least one 95% credible set, plus a cis min-p / Bonferroni comparison.
- **Colocalization.** coloc 5.2.3.
  - Simulate a GWAS from the *true* genotypes (large effective n, summary statistics with reference LD).
  - Half the genes share the eQTL causal variant (colocalization should be supported); half have a nearby,
    different causal in moderate LD, r² 0.3–0.8 (it should not be).
  - Run `coloc.abf`, and `coloc.susie` where both traits have credible sets. Each eQTL arm is compared against the
    same GWAS.

## Metrics: each arm vs `truth`, split by h², causal scenario, cohort size and causal-SNP class

1. eGene recall.
2. **Credible-set coverage**: the fraction of true causal variants inside a 95% credible set.
3. Credible-set size (number of SNPs).
4. Causal-variant PIP (mean, and the fraction with PIP > 0.5).
5. **Lead-SNP shift**: how often the top-PIP SNP is an observed (exonic / peak) site when the causal is untyped,
   vs the same rate under truth.
6. Coloc: posterior for a shared variant (PP.H4) when the causal is shared (sensitivity), and when it is not
   (false colocalization rate).
7. Effect-size bias at the causal variant (dosage posterior means vs hard calls).

## Go / no-go

- **Story (go to tier 1 / real data):** credible-set coverage of untyped noncoding causals drops by ≥ 15 points
  vs truth (e.g. 0.90 → 0.75), **or** lead SNPs shift toward observed sites at ≥ 2× the truth rate, **or**
  colocalization sensitivity or false-colocalization rate changes by ≥ 10 points, at realistic h² and n.
- **No story:** losses within a few points. The paper message is then positive ("latent genotypes are sufficient
  for eQTL discovery and fine-mapping"), with the multiome / low-pass arms showing what closes any residual gap.

## Cost estimate

**Compute.** GLIMPSE2 cost is scaled from the chr20 test: ~89 CPU-min for one ~16-donor run of chr20 (64 Mb) at
default iterations, ~2.5× faster with b2m5. That gives roughly 0.035–0.09 CPU-min per donor-Mb.

| step | CPU-h |
|---|---|
| Level A smoke test (noise injection + susie) | < 5 |
| Region panels (leave-out, 100 Mb) | ~1 |
| Allele-count emulation + GLs (all arms) | ~2 |
| GLIMPSE2: 300 donors × 100 Mb × 3 imputed arms (+1 if low-pass) | ~50–150 (+20–50) |
| susieR: ~16k fits at ~3 s (n = 300, ~2–3k SNPs per window) | ~15 |
| coloc | < 1 |
| **total** | **~75–225 CPU-h** |

At `-tc 40`, that is roughly half a day of wall time including queue waits. Disk is a few GB in scratch (region-
restricted BCFs).

**Analyst time.** About 3–4 working days:
- day 1: cohort, regions, emulator, Level A smoke test (gives a first answer);
- day 2: GLIMPSE2 arms and the susie / coloc pipeline running;
- days 3–4: metrics, figures, write-up, go/no-go call.

## Scripts (planned; naming continues A05)

| script | does |
|---|---|
| `A06a_finemap_cohort.py` | test donors (EUR, pedigree-filtered), leave-out panel lists, 100 gene regions + cis SNP sets |
| `A06b_emulate_latent_gl.py` | per-donor depth profiles from real souporcell clusters → allele counts with ambient → GLs (reuses `souporcell_to_gl.py`) |
| `A06c_glimpse_regions.sh` | qsub: GLIMPSE2 per arm × region chunk against the leave-out panel |
| `A06d_simulate_eqtl.R` | phenotypes (h² grid, causal scenarios, seeds) + simulated GWAS summary statistics |
| `A06e_finemap_coloc.R` | susieR + coloc per gene × arm; writes one long result table |
| `A06f_score.py` | metrics above → `results/A06_finemap/` |

## Known limitations (state them in any write-up)

- Errors are independent of the phenotype (no circularity / ASE coupling); this is by design.
- The EUR-only cohort and the reduced EUR panel will make imputation slightly pessimistic.
- Single-ancestry LD; fine-mapping behaviour differs in multi-ancestry cohorts.
- A pseudobulk phenotype does not model cell-type composition or per-cell noise beyond h².
