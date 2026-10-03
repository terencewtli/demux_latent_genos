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

## 2026-09-30 ~10:00: ambient-aware GLs are the best arm (chr20 test)

Same test as before (run 11 = n8 greedy rep1 GEX, chr20), scored with `score_latent_imputation.py`. Mean over 8 donors:

| | naive souporcell GT | GLIMPSE2 binomial | **GLIMPSE2 ambient** | binomial, fewer iters | Minimac4* |
|---|---|---|---|---|---|
| typed sites | 0.716 | 0.911 | **0.952** | 0.910 | 0.878 |
| 1–2 reads | 0.580 | 0.916 | **0.952** | 0.911 | 0.803 |
| 3–5 reads | 0.679 | 0.929 | **0.962** | 0.928 | 0.857 |
| 6–10 reads | 0.791 | 0.937 | **0.965** | 0.936 | 0.908 |
| > 10 reads | 0.875 | 0.857 | 0.926 | 0.858 | 0.971 |
| untyped, MAF 5–50% | – | 0.524 | **0.538** | 0.515 | 0.453 |

\* Minimac4 used A04b's monomorphic filter (4,112 vs 4,245 sites); its naive baseline on that set is 0.726.

- The binomial arm's dip above 10 reads was **ambient RNA** (28% for this run); the ambient model recovers most of it. Minimac4 is still best at that depth only.
- **The bug that blocked the first ambient attempt was not the model**: `GLIMPSE2_ligate` v2.0.0 has no `--log`. A04g now writes ligate output to a file and reuses finished chunks on rerun.
- **Cost, and where it goes** (binomial, 1 thread, 89 min for chr20): 26 chunks × ~63k panel variants; ~95% of time is HMM imputation, 21 iterations at ~5.5 s, ~1,900 reference haplotypes. 4 threads → 33.5 min (2.65×), 1.8 CPU-h, 15 GB peak.
- **Speedups:** `--burnin 2 --main 5` measured at 2.5× for a loss of 0.001 typed / 0.009 untyped r² — **recommended**. Untested: `--Kpbwt 1000` (~1.7–1.9×), panel restricted to MAF ≥ 1% SNVs (~2–2.5×, ~4× less disk), SNVs only (~1.15×).
- **Full-run estimate** (12 runs × chr2–22): ~1,150–1,400 CPU-h at defaults, **~450–560 CPU-h with fewer iterations**; panels add ~63 CPU-h and ~306 GB scratch (479 GB free). At `-tc 40` (160 cores) that is roughly 4–6 h of compute, 6–10 h including queue waits.
- **Not launched** — awaiting the go-ahead.
- Still open: `score_imputation.py` scores 0 sites under bcftools 1.11 (`query -v` = `--vcf-list`), so A04e's leakage scores are invalid and need a rerun; chr1 waits on the rebuilt A04a sites file.

## 2026-09-30 ~15:20: full GLIMPSE2 run launched (user go-ahead)

- **Arm:** ambient-aware GLs with `--burnin 2 --main 5` (`ARM=ambient`, `SUFFIX=_b2m5`) → `results/A04g_glimpse_impute/ambient_b2m5/`.
  All 12 runs × chr1–22.
- **A03b** (14983421, new wrapper `scripts/qsub/A03b_souporcell_dosage_r2.sh`): fills in assign.tsv for the 6 runs
  that lacked it (the 5 ATAC runs + n16 random rep2 GEX). Scoring needs it.
- **A04f panels,** one job per pool (`-tc 3` each, heavy I/O on the 1000G source): 14983422 / 24 / 26 / 29 / 31 / 33.
- **A04g,** one job per pool covering its GEX + ATAC rows (`-tc 7` each; held on its pool's A04f + A03b):
  14983423 / 25 / 28 / 30 / 32 / 34.
- **Scratch:** ~306 GB of panels needed, 2.9 TB free; nothing is deleted.
- **chr1:** the panel uses the clean A04h download (A04a chr1 sites rebuilt 09-30 03:45). Truth still ends at 190.67 Mb.
- **Estimate:** ~450–560 CPU-h, roughly 6–10 h wall time.
- **Still open:** the `score_imputation.py` bcftools-1.11 bug (A04e leakage scores invalid); the TOPMed leakage download.

## 2026-09-30 ~22:45: TOPMed leakage results downloaded; score_imputation.py fixed; GLIMPSE2 progress

- **TOPMed job-20260929-163357-077 downloaded** with imputationbot into `results/topmed_leakage/topmed/`; zip md5 OK.
  Server QC: 13,350 / 14,374 sites kept (89.87% reference overlap), 0 strand flips, 1,017 monomorphic, 1,506 typed-only.
- **`lib/score_imputation.py` fixed.** bcftools 1.11 `query` has no `-v` / `-m` / `-M` (`-v` is `--vcf-list`), and
  the exit status was never checked, so the old version silently scored 0 sites.
  - Site filters now go through `bcftools view`, piped to `query`.
  - Any non-zero exit raises an error.
  - The A04e local scores (all-zero) are left in place but are invalid.
- **`qsub/A04j_leakage_score.sh`** (14983558, queued) scores all three chr20 imputations with the fixed scorer into
  `results/topmed_leakage/scores_v2/`: TOPMed, local with_donors, local loo.
- **GLIMPSE2 full run** (ambient_b2m5): 60 chromosome imputations done, 59 scored so far. About 26 A04g tasks
  are running; 2 pools' panels are still queued.

## 2026-10-01 ~09:00: TOPMed leakage CONFIRMED; GLIMPSE2 full run at 215 / 264 (two bugs fixed); interim accuracy

**TOPMed leakage test** (A04j, fixed scorer; chr20, 69 donors, 1 in 10 common SNVs as input):

| group (not in the input) | TOPMed r3 | local panel WITH donors | local panel, donors + relatives left out |
|---|---|---|---|
| MAF < 0.1% | **0.983** | 0.876 | 0.536 |
| 0.1–0.5% | 0.994 | 0.937 | 0.805 |
| 1–5% | 0.998 | 0.966 | 0.910 |
| 5–50% | 0.999 | 0.990 | 0.962 |
| private to the 69 donors (in 1000G) | 0.999 | 0.986 | 0.953 |
| **private singletons** | **0.976** (mean per site 0.997) | 0.781 | 0.349 |

- Variants carried only by one donor within 3,202 1000G samples are imputed almost perfectly by TOPMed: better
  than our panel that *contains* the donors, far above leave-out.
- **The 1000G donors are in (or duplicated in) the TOPMed r3 panel.** The TOPMed arm is not usable for these
  simulations. Local leave-pool-out stays the main arm.
- A clean TOPMed comparison would need re-simulated pools from donors certainly absent from TOPMed (e.g. HGDP), and
  even that would need its own leakage test.

**GLIMPSE2 full run** (ambient_b2m5): 215 / 264 chromosome × run imputed and scored. Two failures, both fixed:
1. **`GLIMPSE2_ligate`: "Three files overlapping"** on chr1 (~125 Mb centromere) and chr22 (~27.8 Mb) in all 12
   runs.
   - `--sequential` chunking makes chunk i−1 and i+1 input regions overlap there.
   - `--recursive` fixes chr22 but not chr1 (tested on the sites files).
   - **Fix in A04g:** split the ligate list wherever chunk i−1 and i+1 overlap, ligate each segment, trim it to its
     chunks' output regions (ligate keeps the buffers), then `bcftools concat --naive`.
   - Verified on final greedy rep1 GEX chr22: 2 segments (1–29.86 Mb, 29.86–50.8 Mb); typed r² 0.946.
2. **"Illegal instruction"** in GLIMPSE2_split_reference (static AVX2 build) on older nodes: 5 panels in the two
   ambisim_final pools (greedy chr3 / 11 / 14, random chr15 / 17). A04f / A04g now request
   `-l arch=intel-gold*|intel-E5-2650|intel-6736p`.
- Resubmitted: A04f 14988524 / 25 / 26 (5 panels) → A04g 14988527 (all 264; finished tasks skip, failed ones reuse
  their completed chunks).

**Interim accuracy** (mean over donors and scored chromosomes, 14–20 chromosomes per run):

| | naive souporcell GT | GLIMPSE2 ambient, typed sites | untyped MAF 5–50% |
|---|---|---|---|
| ATAC runs (6) | 0.66–0.73 | 0.91–0.93 | **0.73–0.75** |
| GEX runs (6) | 0.65–0.72 | 0.91–0.95 | **0.56** |
| all 12 | 0.685 | 0.927 | 0.650 |

- Typed-site r² is consistent with the chr20 test (0.95).
- **ATAC imputes untyped sites much better than GEX** (0.74 vs 0.56), presumably because ATAC sites are spread
  genome-wide while GEX is restricted to expressed genes. Combining GEX + ATAC from the same multiome nuclei is an
  obvious next arm.

## 2026-10-01 ~09:40: Minimac4 arm fixed and extended to all 12 runs

- **State before:** A04c had run on 6 runs (those with A03b assignments on 09-29), 20 / 22 chromosomes each, in both
  QC sets, with no genome-wide scoring.
  - **chr1:** "VCF parse error" in Eagle. A04c built its leave-out phasing reference from the corrupt by_chrom chr1.
    It now uses the clean A04h download, as A04f does. The chr1 msav was already rebuilt from it on 09-30.
  - **chr15:** Minimac4 "not enough target variants". The first 20 Mb chunk (the acrocentric p-arm) has 0 / 18,679
    typed sites. A04c now passes `--min-ratio-behavior skip`, so such chunks are skipped instead of failing.
- **Resubmitted:** A04b 14988792 (all 12 rows; existing targets skip) → A04c 14988794 (264 tasks; existing doses
  skip) → **A04k 14988795** (new `qsub/A04k_score_minimac4.sh`).
- A04k scores both QC sets with the same scorer, truth and panel sites as GLIMPSE2, so the arms are directly
  comparable. Output: `results/A04k_minimac4_score/<qc>/<tag>/`.
- **The method comparison now rests on Minimac4 + GLIMPSE2** (local 1000G leave-pool-out). TOPMed would need donors
  outside its panel.

## 2026-10-01 ~19:00: GLIMPSE2 ambient_b2m5 complete genome-wide (12 runs × chr1–22)

All 264 run × chromosome imputations are done and scored (A04g 14988527, no errors). Summary tables:
`results/A04g_glimpse_impute/summary_{run,depth,donor}.ambient_b2m5.tsv`. Means are weighted by sites over
chromosomes. "Untyped" is scored only at MAF ≥ 5%: the panel-site scoring has 0 sites in the lower-MAF bins.

| | naive souporcell GT | GLIMPSE2, typed sites | untyped MAF ≥ 5% | concordance naive → GLIMPSE2 (called sites) |
|---|---|---|---|---|
| ATAC (6 runs; ~330k typed sites / donor) | 0.683 | **0.918** | **0.744** | 0.825 → 0.953 |
| GEX (6 runs; ~195k typed sites / donor) | 0.685 | **0.938** | **0.561** | 0.798 → 0.965 |

- **Pool size:** n16 pools have a lower naive r² (0.66–0.67) than n8 (0.72–0.73). After imputation the gap mostly
  closes (typed r² 0.914–0.942 vs 0.923–0.946).
- **Untyped:** ATAC is better than GEX (0.74 vs 0.56) in every run, as in the interim numbers. ATAC has more typed
  sites, spread genome-wide; GEX sites cluster in expressed exons.
- **Depth** (typed r², GLIMPSE2 / naive):

| reads | ATAC | GEX |
|---|---|---|
| 0 (imputed from neighbours only) | 0.852 / – | 0.916 / – |
| 1–2 | 0.887 / 0.510 | 0.937 / 0.604 |
| 3–5 | 0.923 / 0.647 | 0.948 / 0.699 |
| 6–10 | **0.949** / 0.782 | 0.946 / 0.765 |
| > 10 | 0.937 / 0.895 | 0.918 / 0.849 |

  - Imputation gains most at 1–5 reads (+0.28 to +0.38).
  - Above 10 reads, r² still dips in both modalities, so the ambient model removes only part of the high-depth
    penalty. On chr20 Minimac4 was better at > 10 reads.
- **The merged cluster** (n16 random rep1 GEX, c3 = HG00237 + HG00353) is the only outlier: typed 0.58, untyped
  0.47. Imputation cannot repair a mixture, as expected. Every other donor in every run is at typed ≥ 0.90.
- **Open:**
  - rare-variant (MAF < 5%) accuracy on the real runs is unscored. The leakage test shows that is where panels
    differ;
  - Minimac4 A04c is at 169 / 264, then A04k scores it with the same scorer for the head-to-head.
- **Next:** Step 3 (re-demux with true / latent / imputed genotypes), set up in the following entry.

## 2026-10-01 ~19:30: Step 3 set up (re-demux with recovered genotypes); baselines say the gap is GEX doublets

**Design** (`scripts/pool_design/qsub/A05a_demuxlet_latent.sh`, 12 tasks, one per run):
- Three demuxlet arms on the **existing** pool_design pileups, so no new pileup is needed:
  - `latent_GT`: souporcell hard calls;
  - `imputed_GT` and `imputed_GP`: GLIMPSE2 ambient_b2m5, as hard calls and as posteriors (`--field GP`).
- Samples are the souporcell clusters, mapped to donors with A03b `assign.tsv` at scoring time.
- **Sites:** every arm is restricted (exact alleles) to the pileup sites. These are 1000G global AF ≥ 5% biallelic
  SNPs, chosen without regard to who is in the pool (14% of chr22 sites are monomorphic in the n16 random rep1 pool).
  So no arm gets truth-selected sites, and true_GT is an apples-to-apples ceiling.
- Coverage on chr22 for n16 random rep1 GEX: latent covers 4,075 / 102,467 sites; imputed covers all 102,467.
- souporcell leaves its FORMAT tags undeclared, which bcftools 1.11 cannot subset, so the latent VCF is built in awk
  (GT only).
- **ambisim_final greedy_maxmin rep1 has no pileup.** It was never demuxlet'd in pool_design, so it has no true_GT
  baseline. Its 2 tasks exit; it is 10 / 12 runs unless a pileup is made.
- Test: task 1 (n16 random rep1 GEX, the merged-cluster run) is job 15000307. Submit the other 11 after it checks
  out.

**Scoring** (`scripts/A05b_score_demux.py` → `results/A05b_score_demux/{summary,by_depth}.tsv`). Methods are
true_GT (pool_design demuxlet), souporcell `clusters.tsv`, and the three arms. The baseline numbers are already in.

| run (10 with a baseline) | best-guess accuracy, souporcell / true_GT | doublet recall, souporcell / true_GT | singlet recall, souporcell / true_GT |
|---|---|---|---|
| ATAC (5) | 0.999–1.000 / 0.998–1.000 | 0.989–1.000 / 0.989–0.998 | 0.999–1.000 / 0.89–0.92 |
| GEX n8 random | 0.999 / 0.992 | 0.993 / 0.961 | 0.999 / 0.832 |
| GEX n16 (4) | 0.937–1.000 / 0.982–0.986 | **0.80–0.84** / 0.95–0.96 | 0.937–1.000 / 0.77–0.81 |

- **demuxlet with true genotypes over-calls doublets.** It calls ~20% of true GEX singlets DBL, yet its best guess is
  right for 98–99% of them. It has no ambient model, and ambisim adds ambient. So singlet recall mostly measures
  doublet calling. A05b adds `best_guess_acc` (best single-donor guess right, whatever the droplet call) to measure
  genotype informativeness.
- **souporcell already assigns singlets essentially perfectly** (≥ 0.999) in every run except n16 random rep1 GEX
  (0.937, the merged cluster). There is almost no room for imputed genotypes to improve singlet assignment in these
  simulations.
- **The measurable gap is doublet detection in n16 GEX.** souporcell gets 0.80–0.84 vs 0.95–0.96 for demuxlet on
  true genotypes. Step 3's real question is how much of that gap demuxlet on imputed genotypes recovers, and at
  what cost in singlet over-calling.
- The merged cluster cannot be rescued by imputation (as the plan predicted). Check how demuxlet on its
  latent / imputed genotypes splits those cells.

## 2026-10-01 ~23:50: Minimac4 vs GLIMPSE2, genome-wide (A04l)

Minimac4 (A04c Eagle + Minimac4 on souporcell hard calls, QC sets `all` and `server`) finished 264 / 264 and was
scored by A04k with the same scorer, truth and panel as GLIMPSE2. Tables: `results/A04l_compare/{run,depth}.tsv`.
Means are over runs.

| | typed sites / donor | naive r² (own site set) | **typed r²** | **untyped MAF ≥ 5%** | concordance (called) |
|---|---|---|---|---|---|
| ATAC, GLIMPSE2 | 330k | 0.683 | **0.918** | **0.744** | 0.953 |
| ATAC, Minimac4 all | 318k | 0.692 | 0.821 | 0.632 | 0.893 |
| ATAC, Minimac4 server | 242k | 0.736 | 0.862 | 0.638 | 0.912 |
| GEX, GLIMPSE2 | 195k | 0.685 | **0.938** | **0.561** | 0.965 |
| GEX, Minimac4 all | 190k | 0.690 | 0.789 | 0.464 | 0.873 |
| GEX, Minimac4 server | 92k | 0.753 | 0.867 | 0.471 | 0.917 |

- **GLIMPSE2 wins in every run.**
  - Typed sites: +0.05–0.19.
  - Untyped sites: +0.07–0.13.
- **The gap is largest in the n16 pools.** Minimac4 drops (all: 0.73–0.79 GEX, 0.79 ATAC) while GLIMPSE2 holds
  (0.91–0.94). Hard calls from noisier, lower-depth clusters hurt Minimac4; GLIMPSE2 uses the genotype
  likelihoods.
- The `server` QC set keeps fewer, higher-quality sites (GEX 92k vs 190k). It raises Minimac4's typed r² but
  barely moves untyped accuracy.
- **Depth** (typed r², GLIMPSE2 / Minimac4 all):

| reads | ATAC | GEX |
|---|---|---|
| 0 | 0.85 / 0.58 | 0.92 / 0.62 |
| 1–2 | 0.89 / 0.71 | 0.94 / 0.72 |
| 3–5 | 0.92 / 0.80 | 0.95 / 0.81 |
| 6–10 | 0.95 / 0.89 | 0.95 / 0.87 |
| > 10 | 0.94 / **0.97** | 0.92 / **0.94** |

  - GLIMPSE2's advantage is all at low depth.
  - **Above 10 reads Minimac4 is better**, and so are souporcell's naive calls (0.85–0.89). GLIMPSE2 slightly
    over-smooths well-covered sites (the ambient dip).
  - An obvious hybrid: keep the observed genotype where depth > 10 and the GLIMPSE2 posterior elsewhere. It is
    untested.
- **Decision for the real-data arm:** GLIMPSE2 (ambient GLs, b2m5) is the default imputation method. Minimac4 is
  the comparison.
- Not done: the exact-site matched untyped comparison (`untyped_matched.tsv`). It was too slow on the login node
  (millions of sites × 264 files) and was stopped. Run it as a qsub job if needed; the per-donor aggregate above
  already scores the same truth sites.

## 2026-10-02 ~20:30: A06 Level A smoke test (fine-mapping / coloc with noise-injected genotypes)

Plan: `docs/A06_FINEMAP_PLAN.md`. Scripts `scripts/A06a_finemap_cohort.py`, `A06b_level_a_finemap.R` (qsub
`scripts/qsub/A06b_level_a_finemap.sh`, job 15019949, 100 regions × 3 seeds) and `A06c_score_finemap.py`. Tables in
`results/A06_finemap/level_a/`.

**Setup:**
- 300 unrelated 1000G EUR founders; 100 protein-coding genes on chr16–22, TSS ± 500 kb; MAF ≥ 5% (729–3,091 SNPs per
  region).
- Each donor carries the real per-site depth profile of a souporcell cluster (GEX and ATAC from the same real pool
  donor; 78 profiles). Sites with reads: GEX 1.5–6%, ATAC 4–10%.
- Arms are noise-injected at the measured accuracy per depth class:
  - imputed arms: posterior-mean (Berkson) noise;
  - naive arm: GEX hard calls at observed sites only, mean-imputed elsewhere.
- Realized mean per-SNP r²: GEX 0.66, ATAC 0.82, multiome 0.83, naive 0.45.

**Results: untyped (noncoding) causal SNP** (rates over 300 phenotypes; truth / GEX-imp / ATAC-imp / multiome /
naive):

| h² | eGene (Bonferroni) | causal in 95% CS | median CS size | lead SNP r² with causal |
|---|---|---|---|---|
| 0.05 | 0.49 / 0.50 / 0.58 / 0.57 / 0.08 | 0.24 / 0.02 / 0.04 / 0.06 / 0 | 52 / 6 / 25 / 20 / 2 | 0.61 / 0.51 / 0.60 / 0.59 / 0.19 |
| 0.10 | 0.96 / 0.92 / 0.96 / 0.95 / 0.27 | 0.69 / 0.09 / 0.26 / 0.27 / 0 | 39 / 6 / 12 / 13 / 2 | 0.92 / 0.81 / 0.89 / 0.90 / 0.27 |
| 0.20 | 1.00 / 1.00 / 1.00 / 1.00 / 0.36 | 0.98 / 0.22 / 0.44 / 0.43 / 0 | 19 / 3 / 5 / 5 / 2 | 0.96 / 0.92 / 0.95 / 0.94 / 0.29 |

- **Colocalisation, shared causal (PP.H4 > 0.8):** truth 0.77; GEX 0.53, ATAC 0.66, multiome 0.67, naive 0.17.
- **False colocalisation (distinct causal):** 0.09–0.13 in every arm, vs 0.13 under truth.
- **Go/no-go (A06c):** every arm passes "GO".
  - credible-set coverage drop: 0.30 (multiome) to 0.46 (naive);
  - lead SNP on an observed site: 3.7× truth for GEX;
  - colocalisation sensitivity change: 0.10–0.60.

**Reading, with the key caveat:**
1. **eGene discovery is essentially unaffected** for imputed arms (within a few points of truth at every h²), as
   predicted. Skipping imputation (naive) is the real loss: eGene 0.27 vs 0.96 at h² = 0.10.
2. **Fine-mapping looks badly hurt, in a specific way: credible sets become small and confident but miss the causal
   variant.** The lead SNP is still in high LD with the causal (r² 0.81–0.95 vs 0.92–0.96), yet the 95% set shrinks
   from 19–52 SNPs to 3–25 and excludes it. That is miscalibration, not lost signal.
3. **This pattern is very likely exaggerated by Level A's independent per-SNP noise.**
   - SNPs in near-perfect LD get *different* noise here, so susie can "tell them apart" by whichever is least noisy,
     producing falsely small credible sets.
   - Real imputation errors come from copying panel haplotypes, so tightly linked SNPs share their errors and LD is
     largely preserved.
   - Level A therefore cannot separate a real fine-mapping problem from an artefact of the noise model.
4. **Decision:** do not draw conclusions from Level A. Run **Level B** (emulated GLs + GLIMPSE2; correlated, realistic
   imputation error) on the same cohort, regions and phenotype seeds, so the arms are directly comparable. Level A
   becomes the "independent error" reference that shows how much LD-coherent error matters.
5. Multiome and ATAC are consistently better than GEX for fine-mapping and colocalisation. That ordering is likely
   to survive Level B.

**Next:** implement A06d (emulate allele counts + GLs from the same depth profiles, GLIMPSE2 against a leave-out panel,
region-restricted) and rerun A06b / A06c with the Level B dosages.
