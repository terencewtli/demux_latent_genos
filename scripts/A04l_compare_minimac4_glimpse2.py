"""A04l: head-to-head of the two imputation arms on the 12 simulated runs, genome-wide.

Arms (both scored by lib/score_latent_imputation.py against the same truth, panel and leave-pool-out):
  glimpse2   A04g ambient-aware GLs, --burnin 2 --main 5 (results/A04g_glimpse_impute/ambient_b2m5)
  mm4_all    A04c Eagle + Minimac4 on souporcell hard calls, QC set 'all'    (results/A04k_minimac4_score/all)
  mm4_server same, QC set 'server' (TOPMed-server-like site filters)         (results/A04k_minimac4_score/server)
Typed-site sets differ slightly between arms (A04b's monomorphic filter for Minimac4), so each arm's naive r2 is
reported on its own typed set. Untyped sites are also compared on the exact shared site ids (mean per-site r2).
Usage (allcools env): python A04l_compare_minimac4_glimpse2.py
Output: results/A04l_compare/{run,depth,untyped_matched}.tsv
"""
from __future__ import annotations

import glob
import os

import pandas as pd

R = '/u/project/cluo/terencew/claude/project_ideas/latent_genos/results'
ARMS = {'glimpse2': f'{R}/A04g_glimpse_impute/ambient_b2m5', 'mm4_all': f'{R}/A04k_minimac4_score/all',
        'mm4_server': f'{R}/A04k_minimac4_score/server'}
OUT = f'{R}/A04l_compare'


def wmean(d: pd.DataFrame, col: str, w: str) -> float:
    d = d.dropna(subset=[col])
    return (d[col] * d[w]).sum() / d[w].sum()


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    runs = sorted(os.path.basename(p) for p in glob.glob(f'{ARMS["glimpse2"]}/*'))
    rows, depth, matched = [], [], []
    for run in runs:
        mod = run.split('__')[-1]
        for arm, d in ARMS.items():
            x = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(f'{d}/{run}/chr*.score.donor.tsv')])
            rows.append({'run': run, 'mod': mod, 'arm': arm, 'n_chrom': len(glob.glob(f'{d}/{run}/chr*.score.donor.tsv')),
                         'typed_per_donor': x.groupby('donor')['n_typed'].sum().mean(),
                         'r2_naive': wmean(x, 'r2_naive_gt', 'n_typed'), 'r2_typed': wmean(x, 'r2_imputed_all_typed', 'n_typed'),
                         'r2_untyped_maf5': wmean(x, 'r2_untyped_maf5', 'n_untyped_maf5'),
                         'conc_imputed_called': wmean(x, 'conc_imputed_called', 'n_called')})
            y = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(f'{d}/{run}/chr*.score.depth.tsv')])
            y = y.dropna(subset=['r2_imputed_all'])
            for b, g in y.groupby('depth'):
                depth.append({'run': run, 'mod': mod, 'arm': arm, 'depth': b, 'n_typed': g['n_typed'].sum(),
                              'r2_imputed': wmean(g, 'r2_imputed_all', 'n_typed')})
        # untyped on exactly the same sites: glimpse2 vs mm4_all
        for f in glob.glob(f'{ARMS["glimpse2"]}/{run}/chr*.score.untyped.sites.tsv.gz'):
            m = f.replace(ARMS['glimpse2'], ARMS['mm4_all'])
            if not os.path.exists(m):
                continue
            a = pd.read_csv(f, sep='\t', usecols=['id', 'af', 'r2'])
            b = pd.read_csv(m, sep='\t', usecols=['id', 'r2'])
            j = a.merge(b, on='id', suffixes=('_glimpse2', '_mm4')).dropna()
            matched.append(j.assign(run=run, mod=mod))
    run_t = pd.DataFrame(rows)
    run_t.to_csv(f'{OUT}/run.tsv', sep='\t', index=False)
    dep = pd.DataFrame(depth).groupby(['mod', 'arm', 'depth']).apply(lambda g: pd.Series(
        {'n_typed': g['n_typed'].sum(), 'r2_imputed': wmean(g, 'r2_imputed', 'n_typed')})).reset_index()
    dep.to_csv(f'{OUT}/depth.tsv', sep='\t', index=False)
    mt = pd.concat(matched)
    mt['af_bin'] = pd.cut(mt['af'].clip(upper=0.5), [0, 0.1, 0.2, 0.3, 0.5], labels=['5-10%', '10-20%', '20-30%', '30-50%'])
    ms = mt.groupby(['mod', 'af_bin']).agg(n_sites=('id', 'size'), r2_glimpse2=('r2_glimpse2', 'mean'), r2_mm4=('r2_mm4', 'mean'),
                                           frac_glimpse2_better=('r2_glimpse2', lambda s: (s > mt.loc[s.index, 'r2_mm4']).mean())).reset_index()
    ms.to_csv(f'{OUT}/untyped_matched.tsv', sep='\t', index=False)
    print(run_t.groupby(['mod', 'arm'])[['typed_per_donor', 'r2_naive', 'r2_typed', 'r2_untyped_maf5', 'conc_imputed_called']].mean().round(3).to_string())
    print(dep.pivot_table(index='depth', columns=['mod', 'arm'], values='r2_imputed').round(3).to_string())
    print(ms.round(3).to_string(index=False))
    print(run_t.pivot_table(index='run', columns='arm', values='r2_typed').round(3).to_string())


if __name__ == '__main__':
    main()
