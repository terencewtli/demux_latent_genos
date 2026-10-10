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

## 2026-10-03 ~02:50: A06 Level B built and submitted (GLIMPSE2 on emulated latent-genotype GLs)

- **A06d** (`scripts/A06d_emulate_targets.py`):
  - Per test donor, at every site in its assigned real cluster's souporcell site list, alt reads ~ Binomial(real
    cluster depth, p_g) with the A04g ambient model.
  - ρ is the run's own souporcell ambient fraction: 0.21 ATAC, 0.31–0.32 GEX.
  - GLs come from `lib/souporcell_to_gl.py` functions (identical to A04g). Multiome = GEX + ATAC counts summed.
  - Regions or arms with no listed biallelic 1000G SNV (e.g. r009, pericentromeric chr16) get a `.none` marker and
    no-information dosages (2 × AF).
- **A06e** (`scripts/qsub/A06e_glimpse_regions.sh`):
  - Region-restricted leave-out panel: 1000G 30x minus the 300 test donors and 88 relatives, 2,814 kept (245 EUR).
  - GLIMPSE2_split_reference with input TSS ± 1 Mb and output TSS ± 500 kb (the binary is named after the INPUT
    region; fixed after the first test).
  - GLIMPSE2_phase b2m5 per arm, then `scripts/A06e_extract_ds.py`, which matches samples by name and SNPs by
    POS:REF:ALT.
  - Test on r002: 192 / 166 / 326 target sites (GEX / ATAC / multiome), 0 cohort SNPs missing from the output,
    ~3 min per arm on 4 cores. Estimate for the full run: ~60 CPU-h.
- **A06b** now takes a `level` argument. `B` reads `<region>/imp_<arm>.tsv.gz` for the imputed arms. Phenotypes, the
  GWAS and the seeds are identical to Level A, because they are reseeded per phenotype, so A and B are directly
  comparable.
- **Submitted as one dependency chain:** A06d 15022575 → A06e 15022576 (100 tasks) → A06b level B 15022577 → A06c
  15022579.
  - Outputs: `cluo_scratch/.../latent_genos/A06/level_b_out/`;
  - scored tables: `results/A06_finemap/level_b/`.
- **Next session:** compare `results/A06_finemap/level_b/{summary,paired,coloc,gonogo}.tsv` with `level_a/`.
  - The question is whether the small, causal-missing credible sets survive LD-coherent imputation error.
  - Check `logs/A06e_glimpse_regions.15022576.*` for failed regions first.

## 2026-10-03 ~17:00: Level B vs Level A read; the "high r², low power" paradox was a reporting bug; OneK1K access

- **Bug in A06c:** the summary's `realized_r2` uses `.first()` per arm, so it reports region r000 only (GEX 0.77). The
  real mean over the 100 regions:

| arm | Level A | Level B |
|---|---|---|
| GEX | 0.62 | **0.37** |
| ATAC | 0.81 | 0.72 |
| multiome | 0.81 | 0.74 |

  Per-region table: `results/A06_finemap/level_b/realized_r2_by_region.tsv` (new file; existing outputs untouched).
  **Fixed 2026-10-03:** A06c averages over regions and also reports the median (`realized_r2_median`). Reran
  A06c for both levels (jobs 15029737 / 15029742); only those columns of `summary.tsv` change.
  - Level B mean / median: GEX 0.369 / 0.409, ATAC 0.722 / 0.773, multiome 0.743 / 0.795.
  - Level A: GEX 0.622, ATAC 0.805, multiome 0.813.
  - **Naive r² fixed as well.** A06b dropped SNPs whose arm column is constant: naive mean-imputed sites and Level B
    `.none` regions gave NaN correlations that `na.rm = TRUE` removed. A06b now counts them as r² = 0, over SNPs
    polymorphic in the truth.
  - A06b skips regions whose output already exists, so the reruns went to new directories,
    `$scratch/latent_genos/A06/level_{a,b}_out_r2fix` (jobs 15031448 / 15031453, A06c 15031454 / 15031455). The
    old `level_{a,b}_out` are untouched; the user can remove them.
  - **Corrected realized r², mean / median across regions:**

| arm | Level A | Level B |
|---|---|---|
| GEX imputed | 0.622 / 0.619 | 0.365 / 0.384 |
| ATAC imputed | 0.805 / 0.805 | 0.715 / 0.773 |
| multiome | 0.813 / 0.811 | 0.735 / 0.790 |
| GEX naive | **0.034** / 0.027 (was 0.465) | 0.034 / 0.027 |

  - Every fine-mapping metric is identical to the previous run (deterministic reseeding); only `realized_r2`
    changed.
  - `level_b/realized_r2_by_region.tsv` was computed from the pre-fix outputs; only the imputed arms of region r009
    differ (NaN there, 0 now).
- **Level B GEX quality is set by how many GEX-covered sites a ±500 kb window has.** Real cluster site lists are
  used; r041 / r074 have 4 GEX sites vs 60–158 ATAC.
  - eGene loss vs per-region r²: r = −0.69 (GEX).

| GEX target sites | regions | r² B | eGene A | eGene B | causal in CS, A | causal in CS, B |
|---|---|---|---|---|---|---|
| ≤ 10 | 5 | 0.09 | 0.76 | 0.11 | 0.07 | 0.06 |
| 11–50 | 9 | 0.18 | 0.75 | 0.36 | 0.10 | 0.09 |
| 51–200 | 48 | 0.34 | 0.81 | 0.47 | 0.11 | 0.19 |
| > 200 | 38 | 0.48 | 0.80 | 0.57 | 0.13 | 0.20 |

- **Reading:**
  - Level A assumed genome-average accuracy everywhere. Level B shows that GEX latent genotypes are only as good as
    the expressed genes near the locus. That is a real property of GEX reads, not a simulator artefact.
  - Both levels hurt fine-mapping vs truth, by different mechanisms:
    - A: small, confident, wrong credible sets (an independent-noise artefact);
    - B: inflated credible sets and lost eGenes where coverage is sparse.
  - ATAC and multiome are much less affected.
  - **Caveat:** A06d uses only the sites in souporcell's selected site list. A real pipeline could pile up all panel
    SNPs covered by reads (cellsnp-lite), giving more sites, so B may understate the information.
    Worth a sensitivity arm.
- **OneK1K is usable for real data:**
  - Raw 10x reads are public in SRA (GSE196830 → PRJNA807386, 75 runs, "pub", released 2022-04-10).
  - Genotypes are open access on Zenodo 7619796 ("OneK1K pseudobulk eQTL dataset", K. Alasoo, CC-BY-4.0):
    `OneK1K.noGP.vcf.gz` (12.7 GB), Minimac4-imputed from arrays, 1,098 samples, GRCh37 contigs, with
    TYPED / IMPUTED flags and R2.
  - The truth for scoring is therefore array-typed sites (TYPED / TYPED_ONLY); imputed truth carries its own error
    (R2). onek1k.org only links summary statistics and points to a contact form for individual-level data.

## 2026-10-06: A07 OneK1K real-data pilot (5 pools, GEX only) submitted

- **Why:** real-data check of the GEX latent-genotype + GLIMPSE2 arm, with the A06 question in mind
  (real eQTL fine-mapping with ~980 donors if the pilot is good). OneK1K is 10x 3' v2 scRNA-seq (no
  ATAC) of PBMCs: 75 pools in SRA (pools 40 and 66 absent), 9-18 donors per pool, 20 NovaSeq lanes per pool,
  ~58 B reads / ~2.9 TB SRA total.
- **Truth:** Zenodo 7619796 `OneK1K.noGP.vcf.gz` (Minimac4 v1.0.2, GRCh37, 1,098 samples `OneK1K_<n>`).
  Only array-typed sites (INFO TYPED / TYPED_ONLY) are used as truth, lifted to hg38 with Picard. Imputed sites
  are excluded: they carry their own error, and their errors are shared with ours because both imputations
  lean on reference-panel LD.
- **GEO extras** (per GSM): `*GenotypeSamples.txt` (pool donor list, IDs `682_683`) and `*Individual_Barcodes.csv`
  (published per-cell donor calls). Which half of the GEO ID equals `OneK1K_<n>` is unknown, so A07g assigns
  clusters to donors genetically against all 1,098 samples and reports which half matches. Per-cell
  concordance with the published calls comes for free.
- **Pools** (span donor count at similar depth, 680-871 M reads per pool): 70 (9 donors), 1 (12), 55 (14),
  11 (17), 19 (18). `txt/onek1k_pilot_pools.txt`, `txt/onek1k_pilot_runs.txt` (100 runs).
- **Pipeline** (`scripts/qsub/A07*`, logs in `latent_genos/logs/`, scratch `$SCR/A07/`):
  - A07a: prefetch + fasterq-dump per run. _1 = I1 8 bp, _2 = R1 26 bp, _3 = R2 98 bp, checked on
    SRR18028385. Only R1/R2 .gz are kept (~64 GB per pool).
  - A07b: STARsolo, CellRanger-like v2 settings, hg38_igvf STAR index, `~/bin/STAR` 2.7.10a (the
    `mapping` env has kb, not STAR).
  - A07c: truth + GEO metadata + 737K-august-2016 whitelist (teichlab mirror).
  - A07d: souporcell with the same settings as the simulation GEX arm (A03a): skip_remap, common variants,
    200 restarts; k = GEO donor count.
  - A07e: full (non-leave-out) 1000G 30x GLIMPSE2 panel. OneK1K donors are not in 1000G, so no leave-out.
  - A07f: GLIMPSE2 ambient_b2m5, the same arm as the simulation headline (GEX typed r² 0.94).
  - A07g: assignment + GEO checks + `score_latent_imputation.py` per chromosome, then a genome-wide
    n-weighted per-donor summary.
- **Jobs (2026-10-06):** A07a 15073118 · A07c 15073119 · A07e 15073122 · A07b 15073123 · A07d 15073125 ·
  A07f 15073126 · A07g 15073127 (holds chained).
- **Expected cost:** ~150-250 CPU-h per pool, ~1k CPU-h total. Scratch: ~320 GB fastq + ~300 GB BAMs.
  The fastqs can be removed after A07b; the user deletes them.
- **Go/no-go:** per-donor typed-site r² of imputed GEX genotypes. Around 0.4 or higher on real PBMC pools → plan
  the 75-pool run plus real eQTL fine-mapping (the A06 question on real data).

### 2026-10-06 ~17:50: correction, the OneK1K VCF is already GRCh38 (no liftover)
- Picard rejected 245,159 of 492,638 typed SNVs (241,535 MismatchedRefAllele). A check of REF against both
  references at the *same* coordinates: hg38 400/400, hg19 ~25%. The `1`..`22` contig names had misled me into
  assuming GRCh37, and the 10-06 entry above is wrong on this point.
- A07c now only renames contigs (1 → chr1) and checks REF == hg38 base. Output `onek1k.typed.b38.vcf.gz`; A07g
  reads it. The double-lifted `onek1k.typed.hg38.vcf.gz`, `lifted.vcf.gz` and `lifted.reject.vcf.gz` are wrong and
  unused, left for the user to delete. The first A07c failure (job 15073119) was a wrong chain path.
- Rerun A07c job: see PROGRESS; A07g's hold was updated to it.

## 2026-10-06 (late): Kockelbergh 2026 (scTAPAS) half-scoops the eQTL use case; does it transfer to latent genotypes?

**Paper** (bioRxiv 2026.09.30.755703, Luo lab Oxford; PDF in `reference/kockelbergh_2026/`). COMBAT 5' 10x v1.1
PBMC scRNA-seq, 140 samples / 124 donors (102 with arrays used).
- Pipeline:
  1. Pools are demultiplexed with souporcell (`--skip_remap`, common variants = TOPMed-imputed GSA sites), and
     clusters are matched to array genotypes by Pearson r. **This step uses known genotypes.**
  2. Each donor's chr6 reads → one BAM.
  3. QUILT2 (read-aware, 1000G 30x panel; nGen 100; ~4 h per chunk, 160 GB).
  4. QC: INFO > 0.8, MAF > 5%, HWE → 76,833 chr6 SNVs.
- Accuracy: per-variant dosage R² median 0.80 vs GSA / TOPMed. INFO_SCORE tracks concordance (R² 0.70-0.83).
  Accuracy falls with fewer cells (log cells R² 0.70; < 500 cells fails) and for MAF < 5%. Variants in
  blood-expressed genes pass QC far more often (β 2.6).
- Downstream: cell-type cis-eQTL (TensorQTL, 6 cell types) recovers 68.6% of array cell-type-eGene pairs with 18%
  of the variants; effect sizes agree (R² 0.94); scTAPAS-only eGenes are replicated in OneK1K (80%).
  Ancestry assignment matches arrays 100%. HLA imputation from MHC SNVs (two-field R² 0.59; beats arcasHLA for
  class I, loses for class II); HLA-DRB1*03 → TRAV12-2 in CD4 T cells.
- Their stated limits: 5' only, so **3' untested**; chr6 only; common variants only; ancestry representation;
  linear-reference bias.

**What that takes from us, and what it doesn't**
- Taken: "imputed scRNA genotypes are good enough for eQTL discovery" on real data. Our A06 Level A/B framing
  ("is it useful downstream?") is now partly answered for **common-variant eGene discovery from 5' scRNA**.
- Not taken:
  1. Genotype-free end to end. They label donors with arrays before imputing. In a cohort with no genotypes,
     clusters are anonymous, so demux errors (doublets, mis-assigned cells, ambient) flow straight into the per-donor
     BAMs. That is our setting; whether souporcell latent assignment is clean enough is the transfer question.
  2. Input representation. They impute from all reads (QUILT2 on BAM, every panel site covered). Our A04g / A06d
     impute from souporcell's cluster AO/RO at **souporcell's common-variant site list only** (the A06 Level B
     caveat). Their result suggests we have been throwing away most of the information.
  3. 3' chemistry (their open limitation) is exactly OneK1K, our A07 pilot.
  4. ATAC / multiome; per-donor accuracy curves vs depth, donors and ambient; fine-mapping and coloc quality
     (A06 found these, not eGene discovery, are what degrade); re-demultiplexing; genome-wide.

**Does it transfer to souporcell latent genotypes? Expectation:** mostly yes, because souporcell singlets are
already ~perfect in our simulations (step 3 baselines). The per-cluster BAM is then nearly the per-donor BAM they
use, and labels are irrelevant to imputation. The risks are doublet / ambient leakage and n16 GEX, where we saw
weaker clustering. What *won't* transfer automatically is our current input (souporcell sites only); that is the
part to change.

**Plan A07h** (OneK1K pilot, after A07d):
- Split the STARsolo BAM by souporcell singlet cluster → per-cluster BAMs.
- Impute three ways against the same full 1000G 30x panel and score at the same OneK1K array-typed sites:
  - (i) current A07f: GLIMPSE2 from souporcell-site GLs;
  - (ii) GLIMPSE2 `--bam-list` on the per-cluster BAM (GLs at every panel site covered by reads);
  - (iii) QUILT2 on the per-cluster BAM (their method; bioconda r-quilt, ~160 GB per chunk, so chr22 / chr6 first).
- Plus (iv) an "oracle-labelled" version: cells assigned with OneK1K's published GEO per-cell donor labels instead
  of souporcell clusters. (iv) vs (ii)/(iii) isolates the cost of genotype-free demultiplexing, which is the
  transfer question.
- Readouts: per-donor and per-variant r² by MAF and INFO; then (if time) eGene recovery on the pilot donors
  against OneK1K's published eQTLs.

## 2026-10-07: A07h (read-based imputation arms) and A07i (cis-eQTL step) built; A07d runtime

- **A07d souporcell:** `qalter` to h_rt 16 h (the 48 h request kept it queued). Running since 10:18 on all 5 pools.
  Settings unchanged: `--skip_remap` with the 1000G AF > 5% common-variant list. Expected 6–14 h per pool.
- **A07h, built** (chr6 + chr22 first):
  - `A07h_split_bams.py`: per-group BAMs, primary reads with MAPQ 255 and a valid CB / UB, one read per
    (CB, UB, start, strand). Two labellings:
    - `soup`: souporcell singlets;
    - `oracle`: GEO per-cell labels.
  - **Gotcha [verified]:** GEO `Individual ID` is `<a>_<b>`, and both OneK1K_<a> and OneK1K_<b> are VCF samples
    (9 / 9 in pool 70). The naming convention is therefore taken from A07g `assign.tsv`; oracle mode needs A07g.
  - **Test, pool 70 chr22 oracle [verified]:** 15.2M reads scanned, 7.6M kept, 1.6M UMI duplicates, 2.3M from
    unlabelled cells; 0.54–1.06M reads per donor.
  - `qsub/A07h_glimpse_bam.sh`: GLIMPSE2 `--bam-list --fasta` against the A07e panel, b2m5, A07f ligation. One chunk
    takes 30 s; GLIMPSE reports 0.38× mean coverage and ~95% of panel sites with no read.
    - Test chunk, chr22:25.0–27.8 Mb: per-donor r² over all 857 array sites is 0.37, and INFO ≈ 1 everywhere
      (uninformative with 9 samples).
    - This is not comparable with scTAPAS's per-variant r² after an INFO filter, hence the per-variant scorer below.
  - `qsub/A07h_quilt2.sh`: QUILT2, scTAPAS settings (diploid, nGen 100, 5 Mb chunks, 500 kb buffer), in stages
    chunks → prep → impute → concat. `PANEL_SUBSET=EUR` is available if memory requires. Not yet tested: the
    `quilt` conda env is installing.
  - `A07h_score.py`: per-variant r² across all pilot donors (pooled over pools) by MAF and by expressed-gene
    membership, plus per-donor r², for the arms A07f / bam_glimpse_{soup,oracle} / bam_quilt_{soup,oracle}.
- **A07i, built:**
  - `A07i_pseudobulk.py`: donor × cell-type pseudobulk from marker-score lineages. Pool 70: CD4T 50%, CD8T 17%,
    NK 12%, B 11%, Mono 6%, other 3.5%.
  - `A07i_eqtl.py`: tensorQTL `map_cis` (nperm 1000, 1 Mb, BH on pval_beta).
    - Phenotypes: log CPM, inverse-normal; covariates: 5 expression PCs + pool.
    - Every arm is compared with OneK1K's array-imputed genotypes (R2 ≥ 0.8) on identical donors and phenotypes:
      eGene recovery, slope r², lead-variant LD.
    - Then the fully genotype-free version (soup labels + soup arms).
    - The `tensorqtl` env is installing.
- **Submitted:** A07h_split 15092042 and A07i_pseudobulk 15092047 (both held on A07g); A07h_glimpse_bam 15092043
  (held on the split). Still to submit: QUILT2 stages, A07h_score, A07i_eqtl.
- **Context, for the record:**
  - The 0.94 GEX r² is at souporcell-typed sites only. Untyped common sites are 0.56 GEX / 0.74 ATAC.
  - Level B mean window r² was GEX 0.37, which is why GEX eGene discovery and fine-mapping both dropped in the
    simulations.
  - A07h tests whether all-reads imputation (scTAPAS-style) closes that gap on real 3' data.

## 2026-10-08: A07d souporcell done (4/5 pools); pre-imputation accuracy (A07g_onek1k_preimpute); reruns

- **Souporcell:** pools 70 / 1 / 55 / 11 finished. The singlet rate is about 87% everywhere, doublets 12-13%.
  Pool 19 hit the 24 h limit in troublet (job 15073125.5).
- **Pre-imputation** (`scripts/A07g_onek1k_preimpute.py`; raw souporcell GT vs array GT at souporcell-called array sites).
  Of ~220-250k souporcell SNVs per pool, ~8-9k are on the array, ~6k called per donor. Clean clusters:

  | pool | clean clusters | r2_gt | concordance |
  |---|---|---|---|
  | 70 | 7/9 | 0.70 | 0.79 |
  | 1 | 12/12 | 0.67 | 0.79 |
  | 55 | 14/14 | 0.72 | 0.80 |
  | 11 | 14/17 | 0.73 | 0.81 |

  Het recall is ~0.48 and hom-alt recall ~0.96: hets read as homozygous at low depth, which is what imputation should
  fix. Per-cell agreement with GEO labels is 0.9985 (pool 1) and 0.9995 (pool 55).
- **Pool 70 failed clustering** ("Non-finite gradient"): donors 211 + 221 merged into one cluster, plus one
  all-donor junk cluster.
- **Pool 11 is clean but k was wrong.** GEO lists 17 donors but only 14 have cells, which left 3 tiny junk clusters
  (Hungarian assigned them to the duplicated `_2` truth samples). Its 0.77 GEO-cell concordance is an ID-convention
  artefact of the 846/847/848 run of IDs plus those duplicates, not bad clustering.
- **Truth VCF:** 20 duplicated samples (`OneK1K_<n>_2`), which is why margin_donor = 0 for those donors.
- **Decision (user):** k = donors with cells in GEO Individual_Barcodes (pool 11: 14, pool 19: 13; pool 19 also has a
  donor with only 67 cells). `txt/onek1k_pilot_pools.txt` was updated; the GEO-list k is kept in
  `onek1k_pilot_pools.geo_k.txt`.
  - Rerun: `qsub/A07d_souporcell_rerun.sh` 15105391 (pools 70 seed 2, 11 k=14, 19 k=13) into
    `souporcell/pool<P>_v2`, reusing the first run's vartrix matrices; the clustering binary is called with --seed.
  - A07f / A07g take `SFX=_v2` / `--suffix _v2`. v2 imputation 15105395 / 15105396 and preimpute 15105393 / 15105394
    are held on the rerun.
- **Left for the user to qdel** (the classifier blocked qdel):
  - 15105237: pool 19 k=18 resume;
  - 15073126.89-110: held, pool 19 v1;
  - the old downstream chain 15073127 / 15092042 / 15092043 / 15092047, held on v1 A07f, which needs resubmitting
    against v2 for pools 70 / 11 / 19. A07h / A07i still read `score/pool<P>` and need the suffix.


## 2026-10-08 (evening): souporcell site spacing vs the array; OneK1K depth vs the simulations; OneK1K eQTL download

- **Site spacing** (`scripts/A07g_site_spacing.py`; `results/A07_onek1k/coverage/site_spacing_pool1_vs_array.txt`).
  Autosomal SNVs only.

  | | souporcell pool 1 | OneK1K array (typed) |
  |---|---|---|
  | sites | 249,165 | 492,638 |
  | median gap | 1.4 kb | 3.0 kb |
  | 90th / 99th percentile gap | 25 / 152 kb | 13 / 33 kb |
  | gaps > 500 kb | 100 | 44 |
  | genome in gaps > 250 kb | 389 Mb (~13%) | 146 Mb (~5%, mostly centromeres) |
  | 100 kb / 1 Mb windows with a site | 79% / 94% | 92% / 94% |

  - Souporcell sites are clumped into expressed 3′ gene ends: the median gap is smaller than the array's, but the
    tail is much longer.
  - About 240 Mb more of the genome than the array sits in big gaps, where imputation falls back to the prior.
  - The array is a designed tag-SNP backbone; souporcell sites are chosen by expression, not by LD coverage.
- **Depth per donor** (`scripts/A07g_cluster_depth.py`; `results/A07_onek1k/coverage/coverage_compare.tsv`).
  Medians over clusters with ≥ 100 singlet cells; AO + RO from `cluster_genotypes.vcf`.

  | | cells / donor | sites with reads | ≥ 5 reads | ≥ 10 reads | median depth (covered) | covered sites at 1–2 reads |
  |---|---|---|---|---|---|---|
  | sim n8 GEX | ~1,017 | 188k | 97k | 45k | 5 | 22% |
  | sim n16 GEX | ~520 | 177k | 53k | 23k | 3 | 46% |
  | sim n8 ATAC | ~1,017 | 326k | 236k | 135k | 8 | 10% |
  | sim n16 ATAC | ~510 | 321k | 149k | 71k | 4 | 30% |
  | OneK1K (pools 1/55/11/70) | 1,140–1,510 | 158–170k | 43–47k | 20–23k | 2 | 51–57% |

  - OneK1K has more cells per donor than any simulation, but per-site depth is **at or below the sim n16 GEX runs**
    (3′ v2 chemistry, PBMC). Its mean depth is higher (10–12) because of a few highly expressed genes.
  - n16 GEX GLIMPSE2 typed r² was 0.91–0.94 and untyped common 0.56, so depth alone doesn't make the simulations
    optimistic for OneK1K. Naive r² agrees: sim n16 GEX 0.65–0.67 vs OneK1K pre-imputation 0.67–0.73.
  - What the simulations leave out (ambisim draws reads from the truth genotypes): allele-specific / monoallelic
    expression, RNA editing, reference mapping bias, and real error structure. All of these push real data down.
  - Expectation for A07g: r² ≥ 0.9 at the ~6–9k array sites souporcell covers, and roughly the sim untyped value
    (~0.5–0.6) across all 492k array sites.
- **OneK1K eQTL tables** (onek1k.org S3, 14 cell types + "All") are downloading to
  `latent_genos/reference/onek1k/published_eqtl/` (background; ~1.7 GB per cell type, probably full summary statistics).
  The integrity check is still to do.

## 2026-10-08 (night): A08 eQTL-anchored evaluation started; first real imputed accuracy (chr22)

- **Downloads:**
  - OneK1K full cis-eQTL tables (onek1k.org S3; 14 cell types plus "All") → `latent_genos/reference/onek1k/published_eqtl/`.
    Coordinates are **GRCh37** (alleles match hg19 at 200/200 sampled SNVs vs hg38 98/200).
  - TenK10K multiome caQTL release (Hugging Face `anglixue/TenK10K_multiome`) → `reference/tenk10k_multiome/`:
    significant summary (243,225 rows), per-cell-type sig/independent tables, SuSiE, coloc, peaks; **GRCh38**
    (REF matches hg38 at 223/223). The full parquet tables (~1.1 TB) were not downloaded. TenK10K raw single-cell data and
    genomes are EGA managed access (EGAS50000001653 / 1654).
- **A08a** (`scripts/A08a_select_eqtls.py`, qsub 15109046): the paper's criterion, taken from powellgenomicslab/onek1k_phase1
  code. Spearman, gene ± 1 Mb, MAF > 5%, qvalue per cell type × chr, significant = **localFDR < 0.05** (the table's `FDR`
  column), conditional rounds 1–5. Leads = min P per (cell type, gene, round): **26,597, the paper's exact count**.
  Loci for fine-mapping: round-1 leads with |z| ≥ 5 in cd4nc / cd8et / nk / bin / monoc = 6,327.
- **A08b** (`scripts/A08b_extract_dosages.py`): per chromosome, matched latent (GLIMPSE2) and OneK1K DS matrices on the
  same donors, plus distance to the nearest souporcell GL site. Bug fixed before any use: `bcftools query -s` returns
  samples in **list order**, not VCF order; the order is now read back from `bcftools view -h -s`.
- **A08c, chr22 test** (pools 1 + 55, 26 donors; scored where OneK1K R2 ≥ 0.8; per-variant r² across donors):

  | distance to nearest souporcell site | n | median r² |
  |---|---|---|
  | 0 (the site itself) | 4,355 | 0.91 |
  | < 5 kb | 31,930 | 0.58 |
  | 5–20 kb | 21,562 | 0.09 |
  | 20–100 kb | 13,472 | 0.03 |

  | distance to nearest OneK1K eQTL lead | n | median r² |
  |---|---|---|
  | lead SNP | 434 | **0.77** |
  | < 10 kb | 14,472 | 0.65 |
  | 10–100 kb | 32,910 | 0.30 |
  | 100 kb–1 Mb | 23,775 | 0.05 |

  - Per-donor r²: typed 0.61, eQTL leads 0.76.
  - LD fidelity around leads (± 250 kb, MAF ≥ 5%): median correlation of r² vectors 0.88, proxy-set Jaccard 0.78.
  - Reading: imputation is very local (accuracy collapses beyond ~20 kb of a read-covered site). But eQTL leads sit
    where the reads are, so they are imputed much better than the genome average: the A06 "bias helps cis-eQTL"
    expectation holds on real data. Caveat: n = 26, so the chance floor for r² is ~0.04.
- **Not done:** A08d SuSiE-RSS R script and the qsubs for A08b/c/d. `A08d_prep_loci.py` is written but untested. Four LD
  arms: OneK1K all (in-sample gold), OneK1K pilot, latent pilot, 1000G EUR founders.
- ASE recorded as an idea only (IDEAS 9); the user does not want it implemented.

## 2026-10-08 (late): A08d SuSiE-RSS + qsubs written and tested on 3 chr22 loci

- `scripts/A08d_susie_rss.R`: susie_rss(z, R, n = 982, L = 10), with R = (1 − 0.1) cor + 0.1 I for every arm;
  `estimate_s_rss` records the LD–z mismatch. `scripts/A08d_summarize.py` compares each arm with the gold `onek1k_all`
  arm on PIP correlation, credible-set recall/precision, lead in CS, CS size and s, split by souporcell site density at
  the lead.
- qsubs: `A08b_extract_dosages.sh` (22 chromosomes; skips existing), `A08c_score_anchor.sh`, `A08d_prep_loci.sh` and
  `A08d_susie_rss.sh` (5 cell types; MAX_LOCI 200, LAMBDA 0.1). None submitted: they wait for the v2 pools (70 / 11 / 19).
- **Test** (cd4nc, 3 chr22 loci, 26 donors). In-sample (1,098) and 1000G EUR LD give the same single credible set
  (PIP r 0.99). Both 26-donor arms are dominated by small-n LD: s ≈ 0.6–0.8, up to 10 spurious credible sets; PIP r is
  0.39 for onek1k_pilot and 0.21 for latent_pilot.
  - Expectation: at ~60 donors the pilot arms will still be small-n-limited, so the imputation effect has to be read as
    latent_pilot vs onek1k_pilot, not vs gold.
  - A stronger λ (0.3–0.5) or restricting to the top-|z| region are the knobs, if needed.
- **The test outputs must be removed before the real run** (26 donors; the qsubs skip existing outputs). The user runs:
  `rm -r results/A08_eqtl_anchor/finemap/cd4nc results/A08_eqtl_anchor/dosage/chr22.*`

## 2026-10-08 (late night): genome-wide scoring of pools 1 + 55; no-imputation eQTL baseline arms

- **A07g scoring for pools 1 + 55 (v1, the two clean pools)**, without waiting for v2: job 15109554 (tasks 2-3).
  It writes `results/A07_onek1k/score/pool{1,55}/`, including `assign.tsv`, which A07i needs.
- **Metric caveat (affects how the simulations read):** the simulation headline "untyped r² 0.56 GEX" is a
  **per-donor r² across variants**. That number has a floor above zero from allele-frequency variance alone (a dosage
  of 2p correlates with the truth across sites). The per-variant r² across donors is the one that matters for
  association. The same scorer writes it (`*.untyped.bins.tsv`, `mean_site_r2`). GEX greedy_maxmin rep1 chr20:
  aggregate 0.54, **per-site mean 0.34**. n16 ATAC chr1 / chr10: 0.78 / 0.75 aggregate, 0.64 / 0.58 per site. The
  real chr22 per-variant numbers (A08c) are therefore in line with the simulations, not worse.
- **Why accuracy collapses beyond ~5-20 kb with GEX** (`coverage/site_spacing_pool1_vs_array.txt`): souporcell sites
  cluster in expressed 3′ ends. The median gap is 1.4 kb (array: 3.0 kb), but p90 is 25 kb, p99 152 kb, and 389 Mb
  lie in gaps > 250 kb (array: 146 Mb). Each site has 1-5 reads per donor, so one site barely constrains the
  haplotype pair. Array imputation works because every typed site is an accurate hard call spread evenly.
- **No-imputation baseline for the eQTL step** (user request): `scripts/A07i_raw_arms.py` + `qsub/A07i_raw_arms.sh`
  write per-cluster DS BCFs in the layout A07i reads, under `results/A07_onek1k/raw/{gt,gl}/pool<P><sfx>/`:
  - `raw_gt`: souporcell `cluster_genotypes.vcf` hard calls at its own sites (non-BACKGROUND SNVs; ./. missing). Tested on
    pool 1 chr22: 5,270 sites, 68% of cluster × site GTs called.
  - `raw_gl`: posterior-mean dosage from the A07f `target.bcf` PLs (the same ambient-model input GLIMPSE2 sees) under
    an HWE prior from the panel's `AF_EUR_unrel`. It has the per-site information + AF prior, but no LD, so imputed − raw_gl
    is what LD-based imputation adds. Pool 1 chr22: 5,135 sites.
  - Submitted: 15109603 (pools 1, 55, v1); 15109605 / 15109606 (pools 70 / 11 / 19, `SFX=_v2`, held on v2 A07f
    15105395 / 15105396).
- **A07i_eqtl.py changes:**
  - arms `raw_gt`, `raw_gl` (also in the soup-labels run of the qsub);
  - `pool_tag()` picks `pool<P>_v2` when `score/pool<P>_v2/assign.tsv` exists. Before this change, A07i would have
    used the bad v1 imputations of pools 70 / 11 once their v1 score dirs existed;
  - missing dosages: keep variants with call rate ≥ 0.5 among the analysed donors and fill the rest with the variant
    mean (`MIN_CALL`). Before, `dropna()` dropped any variant with a missing donor, which would empty the raw arms
    when pools carry different site sets. It is a no-op for the complete imputed arms.
  - The A07i_eqtl qsub hold list now includes A07i_raw_arms.
- **A07g genome-wide results, pools 1 + 55 (v1; job 15109554, ~20 min per pool).** Per-donor r² against OneK1K array
  GT, means over donors (`results/A07_onek1k/score/pool{1,55}/summary_donor.tsv`):

  | pool | donors | souporcell sites on the array per donor | naive GT (raw souporcell) | imputed, same sites | imputed, untyped array sites (~484k) |
  |---|---|---|---|---|---|
  | 1 | 12 | 6,184 | 0.668 | **0.875** | 0.584 |
  | 55 | 14 | 6,758 | 0.722 | **0.922** | 0.607 |

  - This is the clean before/after: same donors, same sites. Imputation lifts covered-site accuracy by about 0.2 r²,
    mainly by recovering hets that low depth reads as homozygous.
  - Untyped sites, per-variant r² across donors (the association-relevant metric; `*.untyped.bins.tsv`
    mean_site_r2): **0.27 (pool 1), 0.28 (pool 55)**; MAF 5-50% 0.27 / 0.28, MAF 1-5% 0.26 / 0.27. The per-donor
    0.58 / 0.61 above is inflated by allele-frequency variance across sites. Chance floor with 12-14 donors ≈ 0.08.
    Simulation equivalent (GEX n8, per-site mean) 0.34, so real data is in line with the simulations.

## 2026-10-09: v2 imputation done; scoring, eQTL and A08 chain submitted

- **A07f v2 (pools 70 / 11 / 19) is complete:** 22 / 22 chromosomes per pool, every log reached `End:`
  (jobs 15105395 / 15105396). The raw arms (A07i_raw_arms 15109603 / 5 / 6) are also done.
- **Fix before submitting (`A07i_pseudobulk.py`, `A07h_split_bams.oracle_labels`):** pseudobulk read
  `score/pool<P>/assign.tsv` and `souporcell/pool<P>` with no `_v2`. That would have failed for pools 70 / 11 / 19
  (no v1 score), or used the bad v1 clusters. It now uses the same `pool_tag()` rule as `A07i_eqtl.py`:
  `pool<P>_v2` when `score/pool<P>_v2/assign.tsv` exists. Checked: v1 behaviour is unchanged (pool 1: 17,036 soup /
  17,677 oracle labels).
- **Submitted:**

  | step | job | hold |
  |---|---|---|
  | A07g score v2 (pools 70; 11 / 19) | 15111325; 15111326 | none |
  | A07i_pseudobulk | 15111327 | A07g v2 |
  | A07i_eqtl (chr6 + chr22) | 15111328 | pseudobulk |
  | A08b dosages | 15111329 | none |
  | A08c score | 15111330 | A08b |
  | A08d_prep_loci | 15111331 | A08b |
  | A08d_susie_rss | 15111332 | A08d_prep + a **user hold** |

  - A07h (BAM-based arms) never ran, so A07i_eqtl runs without the `bam_*` arms.
  - A08d_susie_rss is held by the user hold because it skips loci whose output exists. The 26-donor chr22 test
    outputs in `results/A08_eqtl_anchor/finemap/cd4nc` are still there (the dosage half was already removed). The
    user runs `rm -r results/A08_eqtl_anchor/finemap/cd4nc`, then `qrls 15111332`.

## 2026-10-09 (later): all-pool accuracy; tensorqtl env fixed; A07j spacing-vs-information test

- **All 5 pools scored (A07g; v2 for 70 / 11 / 19).** Means over donors:

  | pool | donors | naive GT | imputed, covered sites | per-variant untyped r² |
  |---|---|---|---|---|
  | 70_v2 | 9 | 0.637 | 0.825 | 0.282 |
  | 1 | 12 | 0.668 | 0.875 | 0.268 |
  | 55 | 14 | 0.722 | 0.922 | 0.279 |
  | 11_v2 | 14 | 0.729 | 0.925 | 0.275 |
  | 19_v2 | 13 | 0.692 | 0.891 | 0.267 |

  "Per-variant untyped r²" is the vote-weighted mean of `mean_site_r2` over MAF bins.
- **A08c with 54 donors** (clean clusters, all 5 pools):

  | where | median r² |
  |---|---|
  | lead SNPs (17.5k) | 0.79 |
  | within 10 kb of a lead | 0.67 |
  | 10–100 kb from a lead | 0.32 |
  | 100 kb–1 Mb from a lead | 0.04 |
  | all scored variants | 0.06 |

  By distance to the nearest souporcell site: 0.89 at the site, 0.43 within 5 kb, 0.04 at 5–20 kb.

  LD fidelity around leads: corr_r2 0.94, proxy Jaccard 0.83.
- **A07i_eqtl failed** (15111328): the `tensorqtl` env never existed, because the A00z mamba 0.15 build failed
  silently, like the akita one.
  - **Rebuilt with pip pins:** `conda create -n tensorqtl python=3.10 pip`, then numpy 1.26.4, pandas 2.1.4,
    scipy 1.11.4, pyarrow 14.0.2, torch 2.1.2+cpu and tensorqtl 1.0.9 (`--no-deps`).
  - **Import-time deps:** pandas-plink 2.3.1 (`--no-deps`), dask / xarray 2024.1.1, qtl, deprecated, tqdm,
    zstandard, cffi.
  - **Known gap:** pandas-plink's own pins are unmet (pandas ≥ 2.2, pandera), and pgenlib is absent. Neither
    matters here: A07i passes DataFrames to `cis` and never reads plink / pgen.
  - **Resubmitted** as 15111471. The pseudobulk (15111327) had finished: 62 donors, median 680 CD4 T cells.
- **A07j: why does accuracy collapse within ~5–20 kb?** (user question: are arrays designed to be evenly spread, or
  does 3′ RNA cover too little?) This is a 2 × 2 simulation on the same 54 donors, chr20 + chr22:

  | | the donor's own souporcell sites | random sites, same count, MAF-matched |
  |---|---|---|
  | perfect calls | `perfect_soup` | `perfect_random` |
  | weak calls | `weak_soup` | `weak_random` |

  - **Weak calls:** binomial reads at the donor's real cluster depth, error 0.01, no ambient.
  - **Inputs:** OneK1K genotypes with R2 ≥ 0.8.
  - **Random sites:** drawn from non-typed candidates.
  - **Scoring:** array-TYPED SNVs never used as input by any arm. chr22 has 6,885 such SNVs.
  - **Real latent run on that set:** median 0.06 overall; 0.39 at < 2 kb, 0.14 at 2–5 kb, 0.03 at 5–20 kb.
    Consistent with A08c.
  - **Reading:**
    - perfect_random vs perfect_soup = the cost of RNA spacing;
    - weak vs perfect = the cost of depth;
    - weak_soup vs latent = the cost of ambient RNA and cluster errors.
  - **Scripts:** `scripts/A07j_targets.py` (both chromosomes built: per donor ~2,300 (chr22) / 2,700 (chr20)
    input sites), `scripts/qsub/A07j_impute.sh` (job 15111440, 8 tasks; A07f's GLIMPSE2 steps verbatim) and
    `scripts/A07j_score.py`.

## 2026-10-09 (later): A07j result: spacing, not depth, limits latent-genotype imputation

Jobs 15111440 (8 tasks, 20–45 min each); `scripts/A07j_score.py` → `results/A07j_spacing/`. 54 donors; chr20 +
chr22; 18,579 held-out array-typed SNVs (MAF ≥ 1%) that no arm used as input. Median per-variant r² across donors:

| arm | all | MAF 1–5% | MAF 20–50% |
|---|---|---|---|
| perfect_random | **0.49** | 0.24 | 0.69 |
| weak_random | **0.31** | 0.10 | 0.50 |
| perfect_soup | **0.07** | 0.03 | 0.11 |
| weak_soup | 0.042 | 0.017 | 0.076 |
| latent (real run) | 0.037 | 0.018 | 0.058 |

- **Spacing is the dominant cost.**
  - Perfect genotype calls at each donor's own souporcell sites give r² 0.07.
  - The same number of perfect calls at random sites gives 0.49, seven times higher.
  - Even weak calls (the real depth) at random sites give 0.31, far above perfect calls at RNA sites.
  - Souporcell sites are clustered in a few 3′ ends, so most of them carry redundant haplotype information.
- **Depth costs less, and only once spacing is good:** perfect → weak is 0.49 → 0.31 with random spacing, but
  0.07 → 0.04 with RNA spacing.
- **Ambient RNA and cluster errors cost almost nothing:** weak_soup 0.042 vs the real latent run 0.037.
- **By distance to the nearest souporcell site:**

  | distance | perfect_soup | perfect_random |
  |---|---|---|
  | < 2 kb | 0.41 | 0.48 |
  | 2–5 kb | 0.19 | 0.51 |
  | 5–20 kb | 0.045 | 0.48 |
  | 20–100 kb | 0.017 | 0.49 |

  Random spacing is flat because its own inputs are elsewhere. By distance to each arm's OWN nearest input,
  perfect_random still decays (0.58 → 0.32 at 5–20 kb → 0.14 at 20–100 kb), but much more slowly than RNA spacing,
  because the neighbouring random sites constrain the haplotype on both sides.
- **Answer to the user's question:** both are true, but spacing dominates. 3′ RNA covers too little of the genome,
  and in clumps. Arrays win mainly by even spacing (plus accurate calls). Deeper 3′ GEX sequencing would barely
  help; spread-out reads would (ATAC, full-length RNA, low-pass WGS). For cis-eQTL near expressed genes, the RNA
  sites are where they are needed (A08c leads r² 0.79).
- **Caveats:**
  - perfect calls at souporcell sites come from OneK1K's R2 ≥ 0.8 imputed genotypes (~96% of those sites are not
    on the array);
  - random sites are MAF-matched and drawn from non-typed candidates;
  - chance floor ~0.02.

## 2026-10-09 (night): A07i cis-eQTL and A08d SuSiE-RSS results

- **A07i cis-eQTL** (15111471; chr6 + chr22; pseudobulk, 62 donors; souporcell labels;
  `results/A07_onek1k/eqtl/soup/compare.tsv`). Array eGenes recovered:

  | cell type | array eGenes | imputed | raw_gt | raw_gl |
  |---|---|---|---|---|
  | all | 67 | **0.75** | 0.57 | 0.57 |
  | CD4T | 48 | 0.69 | 0.63 | 0.56 |
  | CD8T | 27 | 0.70 | 0.70 | 0.59 |
  | NK | 17 | 0.65 | 0.82 | 0.71 |
  | B | 26 | 0.73 | 0.54 | 0.54 |
  | Mono | 12 | 0.83 | 0.75 | 0.75 |

  - The same lead SNP as the array in 9–21% of cases (all: imputed 0.16, raw 0.08–0.11), but the median LD
    between the two leads is 0.89–1.0.
  - Counts are small (12–67 eGenes per cell type), so the per-cell-type differences between arms are noise. Only the
    "all" contrast (imputed vs raw) is worth reading.
- **A08d SuSiE-RSS** (15111332, 923 loci; `scripts/A08d_summarize.py` rerun after the user removed the stale
  10-08 summaries). Median vs the gold arm (in-sample OneK1K LD, 982 donors):

  | LD arm | PIP r | CS precision | CS size | s_rss |
  |---|---|---|---|---|
  | 1000G EUR | 0.99 | 1.0 | 15 | 0.000 |
  | OneK1K array, pilot donors | 0.87 | 1.0 | 10 | 0.060 |
  | latent (imputed), pilot donors | **0.67** | **0.5** | 6 | 0.104 |

  - By souporcell site density at the lead (0–5 / 6–20 / > 20 sites per 100 kb), latent PIP r is 0.20 / 0.60 / 0.69.
    The array pilot gives 0.59 / 0.88 / 0.87.
  - Same donors, so latent vs array pilot (0.67 vs 0.87) is the cost of imputation. Credible-set recall is 1.0 for
    every arm, but the latent arm adds spurious, smaller credible sets.
