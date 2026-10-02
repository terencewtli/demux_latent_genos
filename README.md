# Latent genotypes from genotype-free demultiplexing

> **Known data issue (2026-09-29): chr1 past 190,673,377 is missing from every ambisim pool truth VCF and from the
> simulated reads** (corrupt 1000G chr1 source; details in the pool_design repo's `docs/pipeline_issues.md`).
> Decision: for the existing simulations the truncated `vcf/` files are the correct truth, because the reads were
> simulated from them. Scoring against them simply has no truth past chr1:190.7 Mb (~2% of the genome), so those sites
> drop out. Imputed dosages there are left unscored. Clean `vcf_fixed/` files exist for any re-simulation. This is
> an exploration, so we proceed with it as is.


Genotype-free demultiplexers such as souporcell cluster cells in a multiplexed
single-cell pool. Along the way they infer a *latent genotype* for each cluster.
How close are those latent genotypes to the donors' true genotypes? And can
reference-panel imputation (TOPMed-style Minimac4) clean them up enough to feed a
genotype-based demultiplexer?

This is open-ended exploration, not a paper (yet). Start narrow: first find out
whether souporcell gets the latent genotypes right at all.

## Motivation

In the Li 2025 multiome demultiplexing benchmark, souporcell was one of the
strongest methods. It still degrades as the number of donors grows (Fig. S1), and
in the lower-coverage version of the experiment it breaks down quickly past a
handful of donors (Fig. S9). Genotype-based methods (demuxlet; the Alvarez 2025
method) hold up better, but they need genotypes that many experiments don't have.

## Questions

1. **Baseline.** How accurate are souporcell's latent genotypes against the true
   genotypes? Break this down by modality (GEX vs ATAC), donor count (8 vs 16),
   per-site coverage, allele frequency, and pool composition.
2. **Imputation.** Does running imputation on the latent genotypes, with the same
   settings as the TOPMed Imputation Server (Eagle2 phasing + Minimac4), improve
   accuracy on average? Imputation can fix wrong calls and fill in sites that had
   no reads.
3. **Closing the loop.** If a genotype-based demultiplexer is given the
   souporcell-derived genotypes (raw or imputed), does it beat souporcell's own
   assignments? How close does it get to the same method given true genotypes?
4. **Design guidance (optimistic).** Which pool sizes and coverages actually need
   external genotypes? Where does "souporcell → impute → genotype-based demux" get
   you as far as having the genotypes?

## Planned / noted for later

- **Monopogen** (single-cell germline SNV calling with LD-based refinement). We
  should run it. The analysis plan is still to be written; see `md/PROGRESS.md`.
- **Low-pass WGS** (~0.5–1x) plus imputation, as the cheapest external genotype
  source to compare against latent genotypes.

## Data

Simulated multiome pools from the pool-design project
(`project_ideas/pool_design`, repo `demux_pool_design`). Reads from 1000 Genomes
donors are simulated with ambisim at 20% ambient contamination and aligned with
cellranger-arc, so the truth genotypes and the true cell-to-donor assignments are
both known. Most of the ~130 pools have 8 donors. The `ambisim_n16/` tree has
16-donor pools, which is where souporcell should start to struggle, so those come
first.

## Scope (2026-10-02)

This repo is now the umbrella for all latent-genotype ideas and for ideas on improving demultiplexing: genetics
from latent genotypes (eQTL / ASE), low-pass WGS as a cheap identity and genotype source, simulators, and
genotype-QC tools. The brainstorm and critique are in `docs/IDEAS.md`.

> **Next step: A06, the eQTL fine-mapping / colocalization go/no-go test** (`docs/A06_FINEMAP_PLAN.md`).

## Layout

| path | contents |
|---|---|
| `md/PROGRESS.md` | What has run and what's next |
| `md/RESULTS.md` | Validated findings only |
| `md/JOURNAL.md` | Dated log |
| `md/LITERATURE.md` | Lit review: has this been done? |
| `docs/ANALYSIS_PLAN.md` | Metrics, the imputation protocol, and pitfalls |
| `docs/IDEAS.md` | Brainstorm: genetics with latent genotypes, low-pass WGS, simulator tiers, other ideas |
| `docs/A06_FINEMAP_PLAN.md` | **Next step**: fine-mapping / coloc go/no-go test (design, metrics, cost) |
| `docs/CONVENTIONS.md` | Paths, environments, script naming |
| `scripts/pool_design/` | Scripts that live in and run from `project_ideas/pool_design/scripts/latent_genos/` |
| `txt/` | Task lists |

## References

- Li et al. 2025, bioRxiv 10.1101/2025.02.06.636969: multiome demultiplexing
  benchmark (Figs. S1, S9)
- Alvarez et al. 2025: genotype-based demultiplexing (collaborator)
- Heaton et al. 2020, *Nat Methods*: souporcell
- Dou et al. 2024, *Nat Biotechnol*: Monopogen
