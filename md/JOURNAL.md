# Journal

## 2026-09-29
- Project started. Scope: accuracy of souporcell latent genotypes, then whether TOPMed-style imputation improves them. Open-ended exploration.
- Wrote `A03a_souporcell.sh` in pool_design (`scripts/latent_genos/qsub/`) and submitted it on 6 pools. The n=16 pools have priority.
- Noted Monopogen for later.

## 2026-09-29 ~22:30: naive dosage r² of souporcell latent genotypes (A03b)

Script `scripts/A03b_souporcell_dosage_r2.py` (working copy in `latent_genos/scripts/`; outputs in `latent_genos/results/A03b_souporcell_dosage_r2/`, summaries mirrored under `results/`). Six finished runs; the rest are still in A03a.

- **Method:**
  - Souporcell `cluster_genotypes.vcf` vs the pool's truth VCF, biallelic SNVs matched by chrom / pos / ref / alt (2.5–3% of souporcell SNVs are not in the truth common-biallelic VCF; no allele swaps).
  - Clusters are assigned to donors by Hungarian matching on −r.
  - Pearson r² per matched donor, from GT hard-call dosage and from GL expected dosage.
- **Mean r² across donors: 0.687 GT / 0.676 GL** (0.690 / 0.679 without the merged cluster below).

| pool | sites / donor | mean r² GT / GL | min assignment margin |
|---|---|---|---|
| n8 greedy ATAC | ~317k | 0.731 / 0.765 | 0.49 |
| n8 greedy GEX | ~182k | 0.722 / 0.733 | 0.52 |
| n8 random GEX | ~180k | 0.719 / 0.734 | 0.51 |
| n16 greedy rep1 GEX | ~165k | 0.672 / 0.649 | 0.44 |
| n16 greedy rep2 GEX | ~164k | 0.674 / 0.648 | 0.43 |
| n16 random rep1 GEX | ~170k | 0.655 / 0.625 | 0.006 (c3 merges HG00237 + HG00353; c10 has no singlets) |

- **Depth is the main driver:**

| mean r² (GT) | 1–2 reads | 3–5 | 6–10 | > 10 |
|---|---|---|---|---|
| ATAC | 0.51 | 0.63 | 0.77 | 0.88 |
| GEX | 0.60 | 0.70 | 0.77 | 0.85 |

  - GL dosage beats GT only at ≥ 3 reads. Souporcell stores integer-rounded GLs.
- **Allele frequency:**

| | EUR MAF < 5% | 5–20% | > 20% |
|---|---|---|---|
| r² (GEX / ATAC) | 0.49 / 0.44 | 0.73 / 0.75 | 0.63 / 0.72 |
| concordance | 0.90 | | 0.73–0.83 |

  - Within the pool, sites with MAF < 10% (typically one het carrier) have r² 0.83–0.87; sites with MAF > 25% have 0.61–0.69.
  - **Reading:**
    - Rare variants: r² is low while concordance is high, because the truth is almost all hom-ref and a few false alt calls dominate r².
    - Common, het-rich sites: r² drops because a het needs both alleles seen, and at 1–2 reads it is often called hom (het undercalling).
- Next: imputation (A04b / A04c, Minimac4 leave-pool-out) should mostly fix the het undercalling at low depth. A genotype-likelihood-aware imputer (GLIMPSE2 or Beagle GL mode, with GLs recomputed from the cluster AO / RO counts) is the natural comparison arm, since Minimac4 takes hard calls.
