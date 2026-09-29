# Progress

## Status (2026-09-29)

- [x] Project scaffold, README, analysis plan (`docs/ANALYSIS_PLAN.md`)
- [x] Lit review (`md/LITERATURE.md`), first pass 2026-09-29; still need the full text of Hartoularos 2023 and the other bioRxiv preprints (abstract-only so far)
- [ ] **A03a souporcell**: job 14965122 (12 tasks: 4 n16 + 2 n8 pools × GEX/ATAC), submitted 2026-09-29
- [ ] Scoring notebook: cluster→donor matching, dosage r², NRC, site yield, depth curve
- [ ] Imputation: choose TOPMed server vs local Minimac4 (1000G 30x, leave-pool-out); also GLIMPSE2 on allele counts
- [ ] Re-demultiplex with latent / imputed genotypes (demuxlet, Alvarez 2025 method)
- [ ] **Monopogen: to run; analysis plan TBD**
- [ ] Low-pass WGS extension (later)

## Gotchas

- `ambisim/` (the chr22-design n=8 grid) has had its BAMs deleted, so the n=8 comparison pools come from `ambisim_final/`.
- `pool_design/tmp.sh` deletes the cellranger BAMs of pools whose demuxlet run is complete. Keep the BAMs of any pool on this project's task list.
- Both `/u/project/cluo` and `cluo_scratch` are about 99% full.
