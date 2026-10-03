# Progress

> **NEXT STEP (decided 2026-10-02): A06, the eQTL fine-mapping / colocalization go/no-go test.**
> Plan: `docs/A06_FINEMAP_PLAN.md`.
> - Does latent-genotype error (GEX / ATAC / multiome + GLIMPSE2) damage susieR fine-mapping and coloc, compared
>   with true genotypes?
> - 300 1000G EUR donors, 100 gene regions, a pseudobulk eQTL simulator.
> - Estimate: ~75–225 CPU-h, ~3–4 working days.
> - A Level A noise-injection smoke test comes first and gives a first answer on day 1.
>
> **Update 2026-10-02 ~20:30: Level A done** (JOURNAL; `results/A06_finemap/level_a/`).
> - eGene discovery is unaffected by imputed latent genotypes.
> - Fine-mapping credible sets become small and miss the causal SNP, and colocalisation sensitivity drops
>   (0.77 → 0.53 GEX / 0.67 multiome).
> - All of this is likely inflated by Level A's independent per-SNP noise, which breaks LD.
> - **Next: Level B** (GLIMPSE2 on emulated GLs), the same cohort and seeds.

## Status (2026-10-01)

**Update 2026-10-01** (JOURNAL 10-01):
- TOPMed leakage is **confirmed**: donor-private singletons r² 0.976 vs 0.35 left out. The current 1000G-based pools
  cannot be used for a TOPMed comparison. The main arm (local 1000G leave-pool-out) is unaffected and leak-free.
- GLIMPSE2 ambient_b2m5 full run: 215 / 264 done. The ligate-segment and AVX2 fixes are resubmitted (A04g 14988527).
- **Decision pending (user):** whether to build a clean TOPMed arm. Plan, if yes:
  1. leakage-test HGDP donors first (true chr20 genotypes thinned, as in A04d; no simulation needed);
  2. only if clean, re-simulate a few EUR-like pools from HGDP donors with ambisim.
  Watch-outs: the gnomAD HGDP+1KG joint callset also contains 1KG (use HGDP samples only), and HGDP itself may be in
  TOPMed r3 (that is what the test checks).

## Status (2026-09-29)

- [x] Project scaffold, README, analysis plan (`docs/ANALYSIS_PLAN.md`)
- [x] Lit review (`md/LITERATURE.md`), first pass 2026-09-29; still need the full text of Hartoularos 2023 and the other bioRxiv preprints (abstract-only so far)
- [ ] **A03a souporcell**: ATAC tasks from job 14965122; the GEX tasks failed on an empty-array bug under `set -u` (fixed) and were resubmitted as 14966956 (tasks 1-11:2)
- [ ] Scoring notebook: cluster→donor matching, dosage r², NRC, site yield, depth curve
- [ ] **A04a panel msav build**: job 14965520 (chr1-22, full 3,202-sample 1000G 30x panel → `reference/topmed/local_1000G_30x/`), submitted 2026-09-29
- [ ] A04b prep targets (12 tasks), to run after A03a + A04a finish; the leave-pool-out lists were tested (n16 random rep1: 16 donors + 8 relatives excluded, 3,178 kept)
- [ ] A04c Eagle + Minimac4 v4.1.6 (264 tasks = 12 runs × 22 chr; QC sets `server` and `all`)
- [ ] **A04e local leakage calibration**: job 14966649 (held on A04a). Same chr20 input imputed locally with the donors in the panel (`with_donors`) and left out (`loo`); scored with `lib/score_imputation.py` (r² by MAF bin + variants private to the donors). The TOPMed output gets the same scoring
- [ ] If TOPMed does leak: we don't need a new genotype source for the main analysis (local leave-pool-out is leak-free). For a clean TOPMed arm, re-simulate a few pools from HGDP donors (public 30x, not NHLBI cohorts) and repeat the leakage test on them. All of Us genotypes are controlled-tier only and can't leave the Workbench
- [ ] GLIMPSE2 on the GL field (A04b already writes it)
- [ ] TOPMed server comparison through imputationbot (installed, token works, panel id `topmed-r3`). Leakage check first (see ANALYSIS_PLAN)
- [ ] **TOPMed leakage test** submitted 2026-09-29: job-20260929-163357-077 (69 unique donors, true chr20 GTs, every 10th AF>5% SNV = 15,880 sites; hg38, eagle, r2Filter 0, meta). Password in `latent_genos/results/topmed_leakage/password.txt` (not in git). Next: download, then compare rare-variant r² against a local run with and without the donors in the panel
- [ ] Re-demultiplex with latent / imputed genotypes (demuxlet, Alvarez 2025 method)
- [ ] **Monopogen: to run; analysis plan TBD**
- [ ] Low-pass WGS extension (later)

## Gotchas

- `ambisim/` (the chr22-design n=8 grid) has had its BAMs deleted, so the n=8 comparison pools come from `ambisim_final/`.
- `pool_design/tmp.sh` deletes the cellranger BAMs of pools whose demuxlet run is complete. Keep the BAMs of any pool on this project's task list.
- Both `/u/project/cluo` and `cluo_scratch` are about 99% full.
