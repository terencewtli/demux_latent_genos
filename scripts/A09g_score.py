#!/usr/bin/env python
"""A09g: score one A09e imputation target against the WGS genotypes of the 4 PGP donors.

Truth: pgp_filt1 WGS, biallelic SNVs, per donor FT = PASS and GQ >= 20 (other donor x site pairs dropped).
Sites: imputed biallelic SNVs (GLIMPSE2 output, 1000G 30x panel sites) that are in the WGS VCF.
Classes per donor x site:
  typed    the donor has >= 1 read at the site in the GLIMPSE2 target (A09e target.bcf, donor DP > 0)
  untyped  everything else (no read for that donor, or site not in the target)
Metrics, per donor x class x panel-MAF bin (RAF from GLIMPSE2), accumulated over chromosomes as sums:
  r2        squared Pearson r of imputed DS vs WGS dosage over donor x site pairs (aggregate r2, the GLIMPSE
            convention; with only 4 donors a per-variant r2 is not meaningful)
  conc      hard-call (round DS) concordance; nonref_conc = concordance where the WGS or the call is non-ref
  naive_r2  typed only: dosage from the target PL argmax (no imputation) at the same pairs
Outputs: results/A09_multiome/score/<target>/{sums.tsv, metrics.tsv}
usage: python scripts/A09g_score.py --target all_multi
"""
from __future__ import annotations

import argparse
import os
import subprocess
from io import StringIO

import numpy as np
import pandas as pd

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
WGS = '/u/project/cluo/terencew/igvf/2023_YR2/multiome/vcf/wgs/pgp_filt1.rm_missing.reheader.vcf.gz'
BCFTOOLS = '/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools'
DONORS = ['C29', 'C37', 'C38', 'C39']
MIN_GQ = 20
BINS = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
LABELS = ['<1%', '1-5%', '5-10%', '10-20%', '20-50%']


def query(args: list[str], names: list[str]) -> pd.DataFrame:
    txt = subprocess.run([BCFTOOLS, 'query', *args], check=True, capture_output=True, text=True).stdout
    return pd.read_csv(StringIO(txt), sep='\t', header=None, names=names, dtype={'chrom': str})


def chrom_sums(target: str, chrom: str) -> pd.DataFrame:
    imp = f'{PROJ}/results/A09_multiome/impute/{target}'
    I = query(['-i', 'TYPE="snp"', '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF[\t%DS]\n', f'{imp}/{chrom}.imputed.bcf'],
              ['chrom', 'pos', 'ref', 'alt', 'raf', *[f'ds_{d}' for d in DONORS]])
    T = query(['-f', '%CHROM\t%POS\t%REF\t%ALT[\t%DP:%PL]\n', f'{imp}/{chrom}.target.bcf'],
              ['chrom', 'pos', 'ref', 'alt', *[f'tg_{d}' for d in DONORS]])
    W = query(['-r', chrom, '-i', 'TYPE="snp" && N_ALT=1', '-s', ','.join(DONORS),
               '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%GQ:%FT]\n', WGS],
              ['chrom', 'pos', 'ref', 'alt', *[f'w_{d}' for d in DONORS]])
    m = I.merge(W, on=['chrom', 'pos', 'ref', 'alt']).merge(T, on=['chrom', 'pos', 'ref', 'alt'], how='left')
    m['maf'] = np.minimum(m.raf, 1 - m.raf)
    m['bin'] = pd.cut(m.maf, BINS, labels=LABELS, include_lowest=True)
    out = []
    for d in DONORS:
        p = m[f'w_{d}'].str.split(':', expand=True)
        gt = p[0].str.replace('|', '/', regex=False)
        ok = (p[2] == 'PASS') & (pd.to_numeric(p[1], errors='coerce') >= MIN_GQ) & gt.isin(['0/0', '0/1', '1/0', '1/1'])
        x = m.loc[ok, ['bin', f'ds_{d}', f'tg_{d}']].copy()
        x['y'] = gt[ok].str.count('1').astype(float)
        x['ds'] = x[f'ds_{d}'].astype(float)
        tg = x[f'tg_{d}'].fillna('0:0,0,0').str.split(':', expand=True)
        x['typed'] = pd.to_numeric(tg[0], errors='coerce').fillna(0) > 0
        pl = tg[1].str.split(',', expand=True).apply(pd.to_numeric, errors='coerce').values
        x['naive'] = np.where(x.typed, np.argmin(np.nan_to_num(pl, nan=1e9), 1), np.nan)
        x['call'] = np.rint(x.ds)
        x['cls'] = np.where(x.typed, 'typed', 'untyped')
        for (cls, b), g in x.groupby(['cls', 'bin'], observed=True):
            nr = (g.y > 0) | (g.call > 0)
            r = {'chrom': chrom, 'donor': d, 'cls': cls, 'bin': b, 'n': len(g), 'sx': g.ds.sum(), 'sy': g.y.sum(),
                 'sxx': (g.ds ** 2).sum(), 'syy': (g.y ** 2).sum(), 'sxy': (g.ds * g.y).sum(),
                 'n_conc': int((g.call == g.y).sum()), 'n_nr': int(nr.sum()), 'n_nr_conc': int((nr & (g.call == g.y)).sum())}
            if cls == 'typed':
                r.update({'nsx': g.naive.sum(), 'nsxx': (g.naive ** 2).sum(), 'nsxy': (g.naive * g.y).sum(),
                      'n_naive_conc': int((g.naive == g.y).sum())})
            out.append(r)
    return pd.DataFrame(out)


def r2(n: pd.Series, sx: pd.Series, sy: pd.Series, sxx: pd.Series, syy: pd.Series, sxy: pd.Series) -> pd.Series:
    cov = sxy - sx * sy / n
    with np.errstate(invalid='ignore', divide='ignore'):
        return cov ** 2 / ((sxx - sx ** 2 / n) * (syy - sy ** 2 / n))


def metrics(S: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    g = S.groupby(keys, observed=True).sum(numeric_only=True).reset_index()
    g['r2'] = r2(g.n, g.sx, g.sy, g.sxx, g.syy, g.sxy)
    g['conc'] = g.n_conc / g.n
    g['nonref_conc'] = g.n_nr_conc / g.n_nr
    if 'nsx' in g:
        g['naive_r2'] = np.where(g.cls == 'typed', r2(g.n, g.nsx, g.sy, g.nsxx, g.syy, g.nsxy), np.nan)
        g['naive_conc'] = np.where(g.cls == 'typed', g.n_naive_conc / g.n, np.nan)
    return g[[*keys, 'n', 'r2', 'conc', 'nonref_conc', *[c for c in ['naive_r2', 'naive_conc'] if c in g]]]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--target', required=True)
    a = p.parse_args()
    out = f'{PROJ}/results/A09_multiome/score/{a.target}'
    os.makedirs(out, exist_ok=True)
    parts = []
    for c in range(1, 23):
        chrom = f'chr{c}'
        if not os.path.exists(f'{PROJ}/results/A09_multiome/impute/{a.target}/{chrom}.imputed.bcf.csi'):
            print(f'{chrom}: not imputed, skipped', flush=True)
            continue
        parts.append(chrom_sums(a.target, chrom))
        print(f'{chrom}: {parts[-1].n.sum()} donor x site pairs', flush=True)
    S = pd.concat(parts)
    S.insert(0, 'target', a.target)
    S.to_csv(f'{out}/sums.tsv', sep='\t', index=False)
    M = pd.concat([metrics(S, ['target', 'cls', 'bin']).assign(donor='all'),
                   metrics(S, ['target', 'cls']).assign(donor='all', bin='all'),
                   metrics(S, ['target', 'donor', 'cls']).assign(bin='all')])
    M['n_chrom'] = S.chrom.nunique()
    M.to_csv(f'{out}/metrics.tsv', sep='\t', index=False, float_format='%.4f')
    print(M.to_string(index=False, float_format=lambda x: f'{x:.3f}'))


if __name__ == '__main__':
    main()
