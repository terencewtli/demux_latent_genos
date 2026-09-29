# Progress

## Status (2026-09-29)

- [x] Project scaffold, README, analysis plan (`docs/ANALYSIS_PLAN.md`)
- [x] Lit review (`md/LITERATURE.md`), first pass 2026-09-29; still need the full text of Hartoularos 2023 and the other bioRxiv preprints (abstract-only so far)
- [ ] **A03a souporcell**: job 14965122 (12 tasks: 4 n16 + 2 n8 pools × GEX/ATAC), submitted 2026-09-29
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
