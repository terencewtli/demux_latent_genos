"""A08d summary: compare SuSiE-RSS results across LD arms, per locus, against the in-sample OneK1K LD ("gold").

For each locus and arm vs onek1k_all:
  pip_r          Pearson r of PIP vectors
  cs_recall      fraction of gold credible sets that share >= 1 variant with some CS of the arm
  cs_precision   fraction of the arm's credible sets that share >= 1 variant with some gold CS
  n_cs, median CS size, lead_in_cs, lead_pip, s_rss (LD-z mismatch; larger = worse)
The key contrast is latent_pilot vs onek1k_pilot (same donors, so the difference is imputation); onek1k_pilot vs
onek1k_all is the small-n cost, kg_eur the usual out-of-sample reference.
Also split by souporcell target density at the lead (n_target_100kb from z.tsv): 0-5 / 6-20 / > 20.

Outputs: results/A08_eqtl_anchor/finemap/summary_locus.tsv, summary_arm.tsv
usage: python scripts/A08d_summarize.py
"""
import glob
import os

import numpy as np
import pandas as pd

OUT = '/u/project/cluo/terencew/claude/project_ideas/latent_genos/results/A08_eqtl_anchor/finemap'
GOLD = 'onek1k_all'


def cs_sets(p: pd.DataFrame) -> list:
    return [set(g.variant) for k, g in p[p.cs > 0].groupby('cs')]


def overlap_frac(a: list, b: list) -> float:
    if not a:
        return np.nan
    return np.mean([any(len(x & y) > 0 for y in b) for x in a])


def main() -> None:
    rows = []
    for d in sorted(glob.glob(f'{OUT}/*/*/')):
        sp, pp, zp = f'{d}susie_summary.tsv', f'{d}susie_pip.tsv.gz', f'{d}z.tsv'
        if not (os.path.exists(sp) and os.path.exists(pp)):
            continue
        s = pd.read_csv(sp, sep='\t')
        p = pd.read_csv(pp, sep='\t')
        z = pd.read_csv(zp, sep='\t')
        dens = int(z.loc[z.is_lead == 1, 'n_target_100kb'].iat[0])
        g = p[p.arm == GOLD]
        if g.empty:
            continue
        gold_cs = cs_sets(g)
        for arm, a in p.groupby('arm'):
            m = g[['variant', 'pip']].merge(a[['variant', 'pip']], on='variant', suffixes=('_gold', ''))
            arm_cs = cs_sets(a)
            r = s[s.arm == arm].iloc[0]
            sizes = [len(x) for x in arm_cs]
            rows.append({'ct': r.ct, 'gene': r.gene, 'arm': arm, 'n_target_100kb_lead': dens, 's_rss': r.s_rss,
                         'n_cs': len(arm_cs), 'n_cs_gold': len(gold_cs),
                         'cs_size_median': np.median(sizes) if sizes else np.nan,
                         'lead_in_cs': r.lead_in_cs, 'lead_pip': r.lead_pip,
                         'pip_r': np.corrcoef(m.pip_gold, m.pip)[0, 1] if m.pip.std() > 0 and m.pip_gold.std() > 0 else np.nan,
                         'cs_recall': overlap_frac(gold_cs, arm_cs), 'cs_precision': overlap_frac(arm_cs, gold_cs)})
    R = pd.DataFrame(rows)
    R.to_csv(f'{OUT}/summary_locus.tsv', sep='\t', index=False)
    R['density'] = pd.cut(R.n_target_100kb_lead, [-1, 5, 20, 1e9], labels=['0-5', '6-20', '>20'])
    cols = ['pip_r', 'cs_recall', 'cs_precision', 'n_cs', 'cs_size_median', 'lead_in_cs', 'lead_pip', 's_rss']
    A = R.groupby('arm')[cols].median().join(R.groupby('arm').size().rename('loci'))
    A.to_csv(f'{OUT}/summary_arm.tsv', sep='\t')
    print(A.round(3).to_string())
    print(R.groupby(['density', 'arm'])[['pip_r', 'cs_recall', 'cs_precision']].median().round(3).to_string())


if __name__ == '__main__':
    main()
