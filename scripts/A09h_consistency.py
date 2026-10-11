#!/usr/bin/env python
"""A09h: consistency of a donor's imputed genotypes across pools (= time points) for one modality and chromosome.

Each pool ys3a..ys3j holds all 4 donors, each at a different time point (multiome/txt/demux_time_map.txt), so a donor
has up to 10 independent imputations (one per pool target ys3?_<mod>) of the same germline genome.
Sites: WGS biallelic SNVs with FT = PASS and GQ >= 20 for that donor (as A09g) that are in the imputed panel output.
Per donor, per pool pair (p, q), per class (both_untyped / one_typed / both_typed; typed = donor DP > 0 in that
pool's target) x panel-MAF bin, sums for:
  r2 of DS_p vs DS_q, hard-call discordance (all sites and where either call is non-ref), and for discordant calls
  which side matches WGS.
Baselines in the same table:
  pair 'xdonor'   DS of donor a in pool p vs donor b in pool q (a != b, p != q): what unrelated genomes give
  pair 'wgs'      WGS of donor a vs donor b (truth-level unrelated baseline)
Per donor x bin also:
  k_wrong   histogram of how many pools miscall the site (systematic vs random error)
  consensus r2 vs WGS of the mean DS over pools, of the single pools (mean), and of all_<mod> (reads pooled, then imputed)
Outputs: results/A09_multiome/consistency/<mod>/<chrom>.{pairs,kwrong,consensus}.tsv
usage: python scripts/A09h_consistency.py --mod multi --chrom chr22
"""
from __future__ import annotations

import argparse
import itertools
import os
import subprocess
from io import StringIO

import numpy as np
import pandas as pd

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
WGS = '/u/project/cluo/terencew/igvf/2023_YR2/multiome/vcf/wgs/pgp_filt1.rm_missing.reheader.vcf.gz'
BCFTOOLS = '/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools'
DONORS = ['C29', 'C37', 'C38', 'C39']
POOLS = [f'ys3{c}' for c in 'abcdefghij']
MIN_GQ = 20
BINS = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
LABELS = ['<1%', '1-5%', '5-10%', '10-20%', '20-50%']
KEY = ['chrom', 'pos', 'ref', 'alt']


def query(args: list[str], names: list[str]) -> pd.DataFrame:
    txt = subprocess.run([BCFTOOLS, 'query', *args], check=True, capture_output=True, text=True).stdout
    return pd.read_csv(StringIO(txt), sep='\t', header=None, names=names, dtype={'chrom': str})


def load_truth(chrom: str) -> tuple[pd.DataFrame, np.ndarray]:
    W = query(['-r', chrom, '-i', 'TYPE="snp" && N_ALT=1', '-s', ','.join(DONORS),
               '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%GQ:%FT]\n', WGS], [*KEY, *DONORS])
    Y = np.full((len(W), len(DONORS)), np.nan)
    for j, d in enumerate(DONORS):
        p = W[d].str.split(':', expand=True)
        gt = p[0].str.replace('|', '/', regex=False)
        ok = (p[2] == 'PASS') & (pd.to_numeric(p[1], errors='coerce') >= MIN_GQ) & gt.isin(['0/0', '0/1', '1/0', '1/1'])
        Y[ok.values, j] = gt[ok].str.count('1').values
    return W[KEY], Y


def load_target(target: str, chrom: str, sites: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    imp = f'{PROJ}/results/A09_multiome/impute/{target}'
    if not os.path.exists(f'{imp}/{chrom}.imputed.bcf.csi'):
        return None
    I = query(['-i', 'TYPE="snp"', '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF[\t%DS]\n', f'{imp}/{chrom}.imputed.bcf'],
              [*KEY, 'raf', *DONORS])
    T = query(['-f', '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n', f'{imp}/{chrom}.target.bcf'], [*KEY, *DONORS])
    m = sites.merge(I.drop_duplicates(KEY), on=KEY, how='left')
    t = sites.merge(T.drop_duplicates(KEY), on=KEY, how='left')
    ds = m[DONORS].apply(pd.to_numeric, errors='coerce').values
    typed = (t[DONORS].apply(pd.to_numeric, errors='coerce').fillna(0).values > 0)
    return ds, typed, m.raf.values


def r2_sums(x: np.ndarray, y: np.ndarray) -> dict:
    return {'n': len(x), 'sx': x.sum(), 'sy': y.sum(), 'sxx': (x * x).sum(), 'syy': (y * y).sum(), 'sxy': (x * y).sum()}


def main() -> None:
    a = argparse.ArgumentParser()
    a.add_argument('--mod', required=True, choices=['multi', 'atac', 'gex'])
    a.add_argument('--chrom', required=True)
    a = a.parse_args()
    out = f'{PROJ}/results/A09_multiome/consistency/{a.mod}'
    os.makedirs(out, exist_ok=True)

    sites, Y = load_truth(a.chrom)
    print(f'{a.chrom}: {len(sites)} WGS biallelic SNVs', flush=True)
    DS, TY, raf = {}, {}, None
    for p in POOLS:
        r = load_target(f'{p}_{a.mod}', a.chrom, sites)
        if r is None:
            print(f'  {p}_{a.mod}: no imputation, skipped', flush=True)
            continue
        DS[p], TY[p] = r[0], r[1]
        raf = r[2] if raf is None else np.where(np.isnan(raf), r[2], raf)
        print(f'  {p}_{a.mod}: loaded', flush=True)
    allp = load_target(f'all_{a.mod}', a.chrom, sites)
    pools = list(DS)
    maf = np.minimum(raf, 1 - raf)
    b = pd.cut(maf, BINS, labels=LABELS, include_lowest=True).astype(str)
    in_panel = ~np.isnan(raf)

    pairs, kw, cons = [], [], []
    for j, d in enumerate(DONORS):
        ok = in_panel & ~np.isnan(Y[:, j]) & np.all([~np.isnan(DS[p][:, j]) for p in pools], 0)
        y = Y[ok, j]
        bb = b[ok]
        for p, q in itertools.combinations(pools, 2):
            x1, x2 = DS[p][ok, j], DS[q][ok, j]
            c1, c2 = np.rint(x1), np.rint(x2)
            nt = TY[p][ok, j].astype(int) + TY[q][ok, j].astype(int)
            cls = np.array(['both_untyped', 'one_typed', 'both_typed'])[nt]
            df = pd.DataFrame({'cls': cls, 'bin': bb, 'x1': x1, 'x2': x2, 'c1': c1, 'c2': c2, 'y': y})
            for (c, bn), g in df.groupby(['cls', 'bin']):
                disc = g.c1 != g.c2
                nr = (g.c1 > 0) | (g.c2 > 0)
                pairs.append({'chrom': a.chrom, 'donor': d, 'pair': 'same', 'p': p, 'q': q, 'cls': c, 'bin': bn,
                              **r2_sums(g.x1.values, g.x2.values), 'n_disc': int(disc.sum()), 'n_nr': int(nr.sum()),
                              'n_nr_disc': int((nr & disc).sum()),
                              'n_disc_p_right': int((disc & (g.c1 == g.y)).sum()),
                              'n_disc_q_right': int((disc & (g.c2 == g.y)).sum())})
        # systematic vs random error: how many pools miscall each site
        calls = np.stack([np.rint(DS[p][ok, j]) for p in pools], 1)
        typed_any = np.any(np.stack([TY[p][ok, j] for p in pools], 1), 1)
        k = (calls != y[:, None]).sum(1)
        df = pd.DataFrame({'bin': bb, 'typed_any': typed_any, 'k': k, 'y_nr': y > 0})
        kw.append(df.groupby(['bin', 'typed_any', 'y_nr', 'k']).size().rename('n').reset_index()
                  .assign(chrom=a.chrom, donor=d, n_pools=len(pools)))
        # consensus mean DS vs single pools vs reads-pooled imputation
        mean_ds = np.stack([DS[p][ok, j] for p in pools], 1).mean(1)
        rows = {'consensus_mean': mean_ds}
        for p in pools:
            rows[f'single_{p}'] = DS[p][ok, j]
        if allp is not None:
            xa = allp[0][ok, j]
            if not np.isnan(xa).any():
                rows[f'all_{a.mod}'] = xa
        for name, x in rows.items():
            for bn in LABELS + ['all']:
                s = bb == bn if bn != 'all' else np.ones(len(y), bool)
                cons.append({'chrom': a.chrom, 'donor': d, 'arm': name, 'bin': bn, **r2_sums(x[s], y[s]),
                             'n_call_conc': int((np.rint(x[s]) == y[s]).sum()),
                             'n_y_nr': int((y[s] > 0).sum()), 'n_y_nr_called': int(((y[s] > 0) & (np.rint(x[s]) > 0)).sum())})

    # unrelated baselines: other donor in another pool (imputed), and WGS vs WGS
    for (i, d1), (j, d2) in itertools.combinations(enumerate(DONORS), 2):
        ok = in_panel & ~np.isnan(Y[:, i]) & ~np.isnan(Y[:, j])
        pairs.append({'chrom': a.chrom, 'donor': f'{d1}-{d2}', 'pair': 'wgs', 'p': '', 'q': '', 'cls': 'all',
                      'bin': 'all', **r2_sums(Y[ok, i], Y[ok, j])})
        for p, q in itertools.permutations(pools, 2):
            ok2 = ok & ~np.isnan(DS[p][:, i]) & ~np.isnan(DS[q][:, j])
            x1, x2 = DS[p][ok2, i], DS[q][ok2, j]
            pairs.append({'chrom': a.chrom, 'donor': f'{d1}-{d2}', 'pair': 'xdonor', 'p': p, 'q': q, 'cls': 'all',
                          'bin': 'all', **r2_sums(x1, x2), 'n_disc': int((np.rint(x1) != np.rint(x2)).sum())})

    pd.DataFrame(pairs).to_csv(f'{out}/{a.chrom}.pairs.tsv', sep='\t', index=False)
    pd.concat(kw).to_csv(f'{out}/{a.chrom}.kwrong.tsv', sep='\t', index=False)
    pd.DataFrame(cons).to_csv(f'{out}/{a.chrom}.consensus.tsv', sep='\t', index=False)
    print(f'{a.chrom}: wrote {out}/{a.chrom}.{{pairs,kwrong,consensus}}.tsv', flush=True)


if __name__ == '__main__':
    main()
