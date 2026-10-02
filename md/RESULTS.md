# Results

## GLIMPSE2 (ambient-aware GLs, burnin 2 / main 5) recovers souporcell latent genotypes (2026-10-01; JOURNAL 10-01 19:00)

- 12 simulated pools (n8 and n16, EUR, GEX and ATAC), genome-wide, scored against truth with a 1000G panel that
  leaves out the pool donors and their relatives.
- Dosage r² at typed sites: **0.68 naive → 0.92 (ATAC) / 0.94 (GEX)**.
- Untyped common (MAF ≥ 5%): **0.74 (ATAC) / 0.56 (GEX)**.
- Tables: `results/A04g_glimpse_impute/`.
- Caveats:
  - a merged souporcell cluster stays bad (0.58);
  - rare variants are not scored yet;
  - Minimac4 comparison pending.
