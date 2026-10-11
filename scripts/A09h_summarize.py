#!/usr/bin/env python
"""A09h summary: combine per-chromosome A09h_consistency.py sums into metrics for one or more modalities.

Outputs (results/A09_multiome/consistency/):
  pairs_summary.tsv      same-donor cross-pool r2 / discordance by mod x cls x bin (+ all), with xdonor and wgs baselines
  pairs_by_pool.tsv      same-donor r2 per mod x donor x pool pair (all classes, all bins), with day labels
  kwrong_summary.tsv     share of sites miscalled by k of n pools, by mod x bin x typed_any x WGS non-ref
  consensus_summary.tsv  r2 / concordance / non-ref recall vs WGS for mean-of-pools, single pools, reads-pooled all_<mod>
usage: python scripts/A09h_summarize.py --mods multi atac gex
"""
from __future__ import annotations

import argparse
import glob

import numpy as np
import pandas as pd

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
OUT = f'{PROJ}/results/A09_multiome/consistency'
TIMES = '/u/project/cluo/terencew/igvf/2023_YR2/multiome/txt/demux_time_map.txt'
S = ['n', 'sx', 'sy', 'sxx', 'syy', 'sxy']


def r2(g: pd.DataFrame) -> pd.Series:
    cov = g.sxy - g.sx * g.sy / g.n
    with np.errstate(invalid='ignore', divide='ignore'):
        return cov ** 2 / ((g.sxx - g.sx ** 2 / g.n) * (g.syy - g.sy ** 2 / g.n))


def load(mods: list[str], kind: str) -> pd.DataFrame:
    parts = []
    for m in mods:
        fs = sorted(glob.glob(f'{OUT}/{m}/chr*.{kind}.tsv'))
        print(f'{m} {kind}: {len(fs)} chromosomes', flush=True)
        parts += [pd.read_csv(f, sep='\t', dtype={'p': str, 'q': str}).assign(mod=m) for f in fs]
    return pd.concat(parts)


def agg(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    g = df.groupby(keys, dropna=False).sum(numeric_only=True).reset_index()
    g['r2'] = r2(g)
    for c, num, den in [('disc', 'n_disc', 'n'), ('nr_disc', 'n_nr_disc', 'n_nr'),
                        ('disc_p_or_q_right', None, None), ('conc', 'n_call_conc', 'n'),
                        ('nonref_recall', 'n_y_nr_called', 'n_y_nr')]:
        if c == 'disc_p_or_q_right' and 'n_disc_p_right' in g:
            g[c] = np.where(g.n_disc_p_right.notna() & (g.get('pair', 'same') == 'same'),
                            (g.n_disc_p_right + g.n_disc_q_right) / g.n_disc, np.nan)
        elif num in g:
            g[c] = g[num] / g[den]
    return g


def main() -> None:
    a = argparse.ArgumentParser()
    a.add_argument('--mods', nargs='+', default=['multi', 'atac', 'gex'])
    a = a.parse_args()

    P = load(a.mods, 'pairs')
    same = P[P.pair == 'same']
    base = P[P.pair != 'same']
    keep = ['mod', 'pair', 'cls', 'bin', 'n', 'r2', 'disc', 'nr_disc', 'disc_p_or_q_right']
    s = pd.concat([agg(same, ['mod', 'pair', 'cls', 'bin']),
                   agg(same.assign(bin='all'), ['mod', 'pair', 'cls', 'bin']),
                   agg(same.assign(cls='all', bin='all'), ['mod', 'pair', 'cls', 'bin']),
                   agg(base, ['mod', 'pair', 'cls', 'bin'])])
    s[[c for c in keep if c in s]].to_csv(f'{OUT}/pairs_summary.tsv', sep='\t', index=False)
    print(s[[c for c in keep if c in s]].to_string(index=False, float_format='%.4f'), flush=True)

    t = pd.read_csv(TIMES, sep='\t')
    day = dict(zip(t.Sample, t.TimeSample2.str.split('_').str[1].str[1:].astype(int)))
    bp = agg(same, ['mod', 'donor', 'p', 'q'])
    bp['day_p'] = [day[f'{p}_{d}'] for p, d in zip(bp.p, bp.donor)]
    bp['day_q'] = [day[f'{q}_{d}'] for q, d in zip(bp.q, bp.donor)]
    bp[['mod', 'donor', 'p', 'q', 'day_p', 'day_q', 'n', 'r2', 'disc', 'nr_disc']].to_csv(
        f'{OUT}/pairs_by_pool.tsv', sep='\t', index=False)

    K = load(a.mods, 'kwrong')
    k = K.groupby(['mod', 'bin', 'typed_any', 'y_nr', 'k']).n.sum().reset_index()
    k['frac'] = k.n / k.groupby(['mod', 'bin', 'typed_any', 'y_nr']).n.transform('sum')
    k.to_csv(f'{OUT}/kwrong_summary.tsv', sep='\t', index=False)

    C = load(a.mods, 'consensus')
    C['arm_type'] = np.where(C.arm.str.startswith('single_'), 'single_pool_pooled', C.arm)
    c = pd.concat([agg(C, ['mod', 'arm', 'bin']), agg(C, ['mod', 'arm_type', 'bin']).rename(columns={'arm_type': 'arm'})
                   .query('arm == "single_pool_pooled"')])
    c[['mod', 'arm', 'bin', 'n', 'r2', 'conc', 'nonref_recall']].to_csv(f'{OUT}/consensus_summary.tsv', sep='\t', index=False)
    print(c[c.bin == 'all'][['mod', 'arm', 'n', 'r2', 'conc', 'nonref_recall']].to_string(index=False, float_format='%.4f'))


if __name__ == '__main__':
    main()
