# Progress

## Status (2026-10-10): A09 souporcell done and mapped to donors; imputation + scoring running

**Decision (user, 2026-10-09):** GEX-only imputation will not get much further (untyped per-variant r² 0.27-0.28,
below the 0.4 bar; A07j shows marker spacing caps it). The ATAC / multiome test is A09, on the user's own data.

**A09: IGVF fibroblast→iPSC multiome**, 10 pools (ys3a..ys3j) × GEX / ATAC = 20 runs (`txt/A09_runs.txt`); 4 unrelated
EUR PGP donors (C29 / C37 / C38 / C39) in every pool, each at a different time point within a pool
(`.../multiome/txt/demux_time_map.txt`).
- Inputs (read in place, never modified):
  - CR-ARC `/u/project/cluo/terencew/igvf/2023_YR2/multiome/mapping/cr_arc/igvf_ref/default/ys3{a..j}/outs` (skip `old/`)
  - ambimux singlets `.../multiome/csv/demux/ambimux/wgs/ambimux_joint_sings.csv` (tab-separated, singlets only)
  - WGS truth `.../multiome/vcf/wgs/pgp_filt1.rm_missing.reheader.vcf.gz` (per-sample FT / GQ)
  - cell state `.../multiome/csv/clusters/rna_leiden_v7.csv`
- souporcell (A09a-c, `souporcell_gpu pipeline`, k = 4, 200 restarts): all 20 runs done 2026-10-10 01:56, outputs in
  `/u/project/cluo_scratch/terencew/claude/latent_genos/A09/souporcell_v2/<pool>_<mod>/`. Ambient RNA 14-69%.
- **A09d** cluster → donor map (WGS genotype concordance, Hungarian; `results/A09_multiome/match/`): 19/20 runs map
  cleanly and agree with the ambimux majority. ATAC concordance with ambimux singlets ≥ 0.987; GEX 0.89-0.98. GEX vs
  ATAC per-barcode agreement 0.73-0.97. **ys3c GEX fails**: one cluster merges C29 + C39 and one is junk → excluded
  from imputation.
- **A09e** GLIMPSE2 (A07f settings, full 1000G panel) on donor-level GLs summed over runs (`A09e_build_gl.py`):
  31 targets × 22 chr (`txt/A09_targets.txt`: all_{multi,gex,atac} + per pool ys3?_{multi,gex,atac}, ys3c gex/multi
  dropped); job **15121200**. all_* targets finished.
- **A09g** scoring vs WGS (FT = PASS, GQ ≥ 20; aggregate r² by typed / untyped and panel MAF bin; naive PL-argmax
  baseline at typed sites): jobs **15122720** (all_*), **15122721** (per pool, held on A09e).
  First number, all_multi chr22: typed r² 0.979 (naive 0.806), **untyped r² 0.910** (MAF bins 0.74 at < 1% to 0.92),
  non-ref concordance 0.97 / 0.82.
- **A09f** souporcell benchmark on ys3j ATAC (22k cells × 282k loci, k = 4, 200 restarts, seed 4): GPU (2080 Ti)
  1 min 19 s cluster + 49 s troublet; souporcell_gpu CPU (8 threads) 14 min 57 s + 10 min 18 s; same optimum
  (log-lik −22,008,342 vs −22,008,498), ARI 0.999 on shared singlets. Rust souporcell 2.4 (8 threads) still
  running (15121108).
- Next: per-pool scores → (a) GEX vs ATAC vs multiome per library, (b) gain from pooling a donor across libraries,
  (c) consistency of a donor's imputed genotypes across pools / time points; compare with the simulation headline
  (typed 0.92 / 0.94) and the OneK1K GEX result (untyped per-variant 0.27, a different metric: with 4 donors only
  aggregate r² is defined).

**OneK1K pilot (complete):**

The pilot was: 5 OneK1K pools (70, 1, 55, 11, 19), GEX only. The pipeline:
1. STARsolo;
2. souporcell, with the A03a settings;
3. GLIMPSE2 ambient_b2m5, against the full 1000G panel;
4. scoring against the OneK1K array.

Pools 70 / 11 / 19 use the v2 reruns. The latest results are in JOURNAL 10-08 to 10-09 (night).

**Decision rule (set 2026-10-06):** if imputed array-typed r² ≳ 0.4, scale to 75 pools and run a real eQTL analysis.
- It passes at the sites souporcell itself covers: imputed 0.83–0.93.
- It fails at untyped sites, scored per variant across donors: 0.27–0.28.
- The per-donor untyped number (0.58–0.61) is inflated by allele-frequency variance, so it is not the one to use.
- **Decision pending (user).**

| result | number | where |
|---|---|---|
| A07g: per-donor r² at covered sites, naive → imputed (5 pools) | 0.64–0.73 → **0.83–0.93** | JOURNAL 10-09 (later) |
| A07g: per-variant untyped r² | **0.27–0.28** (simulations: 0.34) | JOURNAL 10-08 (late night), 10-09 |
| A08c: per-variant r² at OneK1K eQTL leads (54 donors) | **0.79**; 0.67 < 10 kb; 0.04 at 100 kb–1 Mb | JOURNAL 10-09 (later) |
| A08c: r² by distance to the nearest souporcell site | 0.89 at the site; 0.43 < 5 kb; 0.04 at 5–20 kb | JOURNAL 10-09 (later) |
| A07j spacing test: held-out typed SNVs, per variant | perfect + random sites 0.49; weak + random 0.31; perfect + RNA sites 0.07; real run 0.037 | JOURNAL 10-09 (A07j result) |
| A07i cis-eQTL (chr6 + chr22, 62 donors): array eGenes recovered, all cells | imputed **0.75**; raw GT / GL 0.57 | JOURNAL 10-09 (night) |
| A07i: same lead SNP as the array / median lead–lead LD | 0.16 / 0.96 | JOURNAL 10-09 (night) |
| A08d SuSiE-RSS (923 loci): PIP r vs in-sample LD | latent 0.67; same donors with array genotypes 0.87; 1000G EUR 0.99 | JOURNAL 10-09 (night) |

**Reading:**
- The limit is the spacing of the RNA sites, not their depth (A07j). Reads cluster in expressed 3′ ends, so
  imputation is good at and near expressed genes, and close to the prior beyond ~5–20 kb.
- cis-eQTL leads sit where the reads are, so lead SNPs are imputed well (0.79).
- Deeper 3′ GEX would barely help; spread-out reads (ATAC, full-length RNA, low-pass WGS) would.
- This matches the simulation prior: A06 Level B per-region GEX 0.37 vs ATAC 0.72 / multiome 0.74.

## Open

- [x] Go/no-go on scaling OneK1K: no-go for GEX only (user, 2026-10-09)
- [ ] **A09 real-data ATAC / multiome arm** (above): souporcell + donor mapping done; imputation / scoring running
- [ ] Optional, never run: the read-based arms A07h (GLIMPSE2 `--bam-list`, QUILT2 stages; scripts staged)
- [ ] Monopogen: analysis plan TBD
- [ ] Low-pass WGS extension (later)
- ASE: an idea only (IDEAS 9); the user does not want it implemented.
- Re-demux with recovered genotypes: stopped at the baselines (2026-10-01). Souporcell singlets are already ~perfect;
  the only gap is n16 GEX doublets.

## Done (details in JOURNAL)

- 2026-09-29 – 10-01: simulations (12 ambisim pools). GLIMPSE2 ambient_b2m5 typed r² 0.92 ATAC / 0.94 GEX, beating
  Minimac4 (A04l). TOPMed server leakage confirmed (donor-private singletons 0.976 vs 0.35 left out), so the
  local 1000G leave-pool-out panel is the main arm. The chr1 truncation is handled: score against the truncated
  `vcf/`.
- 2026-10-02 – 10-03: A06 fine-mapping simulation. Level A (noise injection) and Level B (GLIMPSE2 on emulated
  GLs); realized r² GEX 0.37 / ATAC 0.72 / multiome 0.74; GEX quality is set by the expressed sites near a locus.
- 2026-10-06 – 10-09: A07 OneK1K pilot (A07a–j) and A08 eQTL-anchored evaluation (A08a–d); results above.
- Kockelbergh 2026 (scTAPAS) half-scoops the eQTL use case with known genotypes for demux, 5′, and QUILT2
  (JOURNAL 10-06 late).

## Gotchas

- A09 BAMs: for every pool except ys3i the BAM is newer than its `.bai` (rewritten 2026-01), and samtools region
  reads fail. A09a reads each BAM through a symlink in `cluo_scratch/.../latent_genos/A09/bam/` with a fresh index.
- `ambisim/` (the chr22-design n=8 grid) has had its BAMs deleted, so the n=8 comparison pools come from `ambisim_final/`.
- `pool_design/tmp.sh` deletes the cellranger BAMs of pools whose demuxlet run is complete. Keep the BAMs of any pool on this project's task list.
- Both `/u/project/cluo` and `cluo_scratch` are about 99% full.
- OneK1K genotype VCF (Zenodo 7619796) is already GRCh38 despite `1`..`22` contigs: rename only, never lift over.
  The OneK1K eQTL tables are GRCh37.
- `bcftools query -s` returns samples in list order, not VCF order.
- Pools 70 / 11 / 19: use the `_v2` outputs (`pool_tag()` picks them automatically).
- The `tensorqtl` env was rebuilt with pip pins (JOURNAL 10-09 later); the mamba 0.15 build fails silently.
