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

## 2026-09-29 ~23:00: chr1 source corrupt; chr1 re-download; A04b / A04c submitted for 6 finished runs; GLIMPSE2 arm in progress

- **The 1000G 30x chr1 source is corrupt.** `/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom/1000G.chr1.vcf.gz` (April 2022) fails `gzip -t` with a crc / length error. It has the same byte size as the EBI file (2,390,515,911), so the bytes were damaged after download.
  - The A04a chr1 task failed ("Number of columns … 2221 vs 3202"). chr2–22 msav + sites built fine.
  - **Knock-on:** pool truth VCFs built from these files stop at chr1:190,673,377, so they lack the last ~58 Mb of chr1 (~2% of the genome). For A03b this only means those sites were not scored. Check whether anything upstream of the ambisim pools (donor selection, read simulation) used this chr1 file.
  - The demux_benchmark file is left untouched.
- **Fix:** `qsub/A04h_download_chr1.sh` re-downloads chr1 from EBI into `latent_genos/reference/1000G_30x/` (size + full gzip -t + md5 recorded). A04a now accepts `IN_OVERRIDE` (qsub -v), and the chr1 A04a task is held on the download.
- **A04b** (rows 1, 5, 7, 9, 11, 12; held on the chr1 A04a, since A04b needs all 22 sites files) → **A04c** (those rows × chr1–22, held on A04b). All at 8 GB × 4 slots.
- **GLIMPSE2 arm (being built; agent):**
  - GLs recomputed from the cluster AO / RO counts (binomial, e = 0.01).
  - Per-pool leave-out reference panels (GLIMPSE2 cannot drop samples at run time) → chunk / split_reference → phase / ligate.
  - Scored with `lib/score_imputation.py`, plus typed-site r² comparable to A03b.
  - First test: n8 greedy GEX, chr20.

## 2026-09-29 ~23:20: chr1 truncation confirmed in every pool truth VCF

- All 130 pool VCFs (`pool_design/ambisim_final/vcf`, 126; `ambisim_n16/vcf`, 4) have **zero chr1 records past 190,673,377**, the exact point where the corrupt 1000G chr1 source stops reading. Table: `results/chr1_truncation_pool_vcfs.tsv` (chr1 > 190.7 Mb record count, plus chr2 first-58-Mb count as reference, ~159k).
- chr22 ends normally (50,807,702) in the files checked; nothing else looks damaged.
- **Impact:**
  - ~58 Mb of chr1 (~2% of the genome, ~150k variants per pool VCF) is missing from truth.
  - Scoring (A03b, and later A04 / GLIMPSE2) simply has no truth there, so no bias.
  - Anything *generated* from these VCFs (donor selection metrics, simulated genotypes / reads) also lacks the chr1 tail. Check with the pool_design owner.
- **Fix:** rebuild the pool VCFs' chr1 from the re-downloaded clean chr1 (`latent_genos/reference/1000G_30x/`, A04h) by subsetting to each pool's donors. Not done; the user decides.

## 2026-09-29 ~23:40: chr1 decision

- ambisim simulates reads from the pool VCF genotypes, so the simulated donors carry no variants past chr1:190,673,377.
- **Truth for all current analyses = the original truncated `vcf/` files** (internally consistent). The chr1 tail has no truth, so it is unscored for A03b, A04 (Minimac4) and GLIMPSE2.
- Do **not** score against `vcf_fixed/`: it has true donor variants that the simulated reads never contained, which would unfairly penalize both souporcell and imputation.
- The imputation panel chr1 (A04a) uses the clean re-download, so imputation near the boundary is not artificially starved.
- `vcf_fixed/` is being built for all 266 pool VCFs (A04i, held on the chr1 download), for future re-simulation only.
- The user will not regenerate the simulations (pool_design is now a simulation library). The issue is flagged in both repos' READMEs.

## 2026-09-30 ~02:30: GLIMPSE2 arm works (chr20 test, run 11 = n8 greedy rep1 GEX)

- **Built** (`scripts/pool_design/`; working copies in `pool_design/scripts/latent_genos/`):
  - `lib/souporcell_to_gl.py`: PL recomputed from each cluster's AO / RO, binomial with e = 0.01, plus an optional ambient-aware model. Souporcell estimates 28% ambient RNA for this run.
  - `qsub/A04f_glimpse_ref.sh`: per-pool leave-pool-out binary panels (13 samples excluded, 3,189 kept).
  - `qsub/A04g_glimpse_impute.sh`: chunk → phase (`--impute-reference-only-variants`) → ligate.
  - `lib/score_latent_imputation.py`.
- **Truth** = the original pool VCF (MAF ≥ 5% SNVs only), so untyped scoring covers the 5–50% bin only.

| chr20, mean over 8 donors | naive souporcell GT | GLIMPSE2 (binomial PL) | Minimac4 (hard calls) |
|---|---|---|---|
| typed sites (4,245) r² | 0.716 | **0.911** | 0.878 |
| typed, 1–2 reads | 0.58 | 0.92 | |
| typed, 3–5 reads | 0.68 | 0.93 | |
| typed, > 10 reads | 0.875 | 0.857 | |
| untyped MAF 5–50% (154,590 sites), aggregate r² | – | **0.524** | 0.453 |

- **Reading:** imputation fixes the low-depth het undercalling, as expected. GLIMPSE2 with proper GLs beats Minimac4 on hard calls. At > 10 reads it is slightly below naive, probably ambient reads that the binomial model ignores. The ambient arm is written but untested.
- **Cost** (1 core, chr20 → × ~39 for chr2–22):
  - Panels: ~10 CPU-h and ~51 GB scratch per pool (~306 GB for 6 pools; scratch has 479 GB free, so build and consume pool by pool).
  - Imputation: ~57 CPU-h per n8 run and ~115 per n16 run (~1,100 CPU-h for 12 runs).
  - Minimac4 is ~12× faster.
- **Bug found:** `lib/score_imputation.py` calls `bcftools query -m2 -v snps`. In bcftools 1.11, query's `-v` means `--vcf-list`, and the exit status is not checked, so it **silently scores 0 sites**. The A04e (TOPMed leakage calibration) scores are therefore invalid; fix and rerun. The new scorer avoids this.
- GLIMPSE2 binaries run from executable copies in `cluo_scratch/.../latent_genos/bin`; b38 maps in `.../glimpse_maps/b38`.
