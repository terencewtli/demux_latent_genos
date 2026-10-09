"""A07j scoring: per-variant r2 across the 54 donors for the four simulated arms (A07j_impute) and the real latent
run (A07f), on the same held-out variants.

Held-out set: OneK1K array-TYPED biallelic SNVs with MAF >= 1%, allele-matched to the 1000G panel, that were never an
input to any arm: not in any pool's A07f target (all_target_sites) and not an A07j input site (input_sites).
Truth = OneK1K DS at typed sites (array genotypes). r2 = squared Pearson correlation across donors of imputed DS
vs truth, per variant (variants monomorphic among the donors are dropped). Chance floor with 54 donors ~ 0.02.
Bins: MAF; distance to the nearest souporcell site (union over pools, the axis of A08c); distance to the nearest
input site of the arm itself (soup arms and latent: souporcell sites; random arms: their random sites).
Output: results/A07j_spacing/{per_variant.tsv.gz, summary.tsv, summary.txt}
usage: python scripts/A07j_score.py [--chroms chr20,chr22]
"""
import argparse
import os
import subprocess
from typing import Dict, List

import numpy as np
import pandas as pd

from A07j_targets import ONEK1K, PANEL, POOLS, PROJ, SCR, pool_clusters, run, sample_order

ARMS = ['perfect_soup', 'perfect_random', 'weak_soup', 'weak_random']
OUT = f'{PROJ}/results/A07j_spacing'
DIST_BINS = [-1, 0, 2000, 5000, 20000, 100000, 1e9]
DIST_LABELS = ['0', '<2kb', '2-5kb', '5-20kb', '20-100kb', '>100kb']
MAF_BINS = [0.01, 0.05, 0.1, 0.2, 0.5]


def heldout(chrom: str, donors: List[str]) -> pd.DataFrame:
    expr = 'TYPE="snp" && N_ALT=1 && INFO/TYPED=1 && INFO/MAF>=0.01'
    order = sample_order(ONEK1K, donors)
    txt = run(['bcftools', 'query', '-r', chrom.replace('chr', ''), '-s', ','.join(donors), '-i', expr,
               '-f', '%POS\t%REF\t%ALT\t%INFO/MAF[\t%DS]\n', ONEK1K])
    d = pd.DataFrame([l.split('\t') for l in txt.splitlines()], columns=['pos', 'ref', 'alt', 'maf'] + order)
    d['pos'] = d.pos.astype(int)
    for c in ['maf'] + order:
        d[c] = d[c].astype(float)
    ps = run(['bcftools', 'query', '-f', '%POS\t%REF\t%ALT\n', f'{PANEL}/1000G_30x.{chrom}.sites.vcf.gz'])
    p = pd.DataFrame([l.split('\t') for l in ps.splitlines()], columns=['pos', 'ref', 'alt'])
    p['pos'] = p.pos.astype(int)
    d = d.merge(p.drop_duplicates(), on=['pos', 'ref', 'alt']).drop_duplicates('pos', keep=False)
    inp = pd.read_csv(f'{OUT}/inputs/{chrom}.input_sites.tsv.gz', sep='\t')
    tgt = pd.read_csv(f'{OUT}/inputs/{chrom}.all_target_sites.tsv.gz', sep='\t')
    used = set(inp.pos) | set(tgt.pos)
    return d[~d.pos.isin(used)].reset_index(drop=True)


def read_ds(bcf: str, chrom: str, pos: pd.DataFrame, rename: Dict[str, str]) -> pd.DataFrame:
    """DS at the held-out sites (pos, ref, alt), columns renamed to donors."""
    keep = [s for s in run(['bcftools', 'query', '-l', bcf]).split() if s in rename]
    reg = f'{OUT}/inputs/.{os.path.basename(bcf)}.regions.tsv'
    pos[['pos']].assign(chrom=chrom)[['chrom', 'pos']].to_csv(reg, sep='\t', header=False, index=False)
    txt = run(['bcftools', 'query', '-T', reg, '-s', ','.join(keep), '-f', '%POS\t%REF\t%ALT[\t%DS]\n', bcf])
    order = sample_order(bcf, keep)
    d = pd.DataFrame([l.split('\t') for l in txt.splitlines()], columns=['pos', 'ref', 'alt'] + [rename[s] for s in order])
    d['pos'] = d.pos.astype(int)
    for s in d.columns[3:]:
        d[s] = d[s].astype(float)
    return pos[['pos', 'ref', 'alt']].merge(d, on=['pos', 'ref', 'alt'], how='left')


def r2_rows(truth: np.ndarray, est: np.ndarray) -> np.ndarray:
    t = truth - truth.mean(1, keepdims=True)
    e = est - np.nanmean(est, 1, keepdims=True)
    num = np.nansum(t * e, 1)
    den = np.sqrt(np.nansum(t ** 2, 1) * np.nansum(e ** 2, 1))
    with np.errstate(invalid='ignore', divide='ignore'):
        return (num / den) ** 2


def nearest(pos: np.ndarray, sites: np.ndarray) -> np.ndarray:
    s = np.sort(np.unique(sites))
    i = np.searchsorted(s, pos)
    lo = np.abs(pos - s[np.clip(i - 1, 0, len(s) - 1)])
    hi = np.abs(s[np.clip(i, 0, len(s) - 1)] - pos)
    return np.minimum(lo, hi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--chroms', default='chr20,chr22')
    a = ap.parse_args()
    rows = []
    for chrom in a.chroms.split(','):
        cmaps = {p: pool_clusters(p, s) for p, s in POOLS.items()}
        donors = sorted({d for m in cmaps.values() for d in m.values()})
        H = heldout(chrom, donors)
        T = H[donors].values
        inp = pd.read_csv(f'{OUT}/inputs/{chrom}.input_sites.tsv.gz', sep='\t')
        tgt = pd.read_csv(f'{OUT}/inputs/{chrom}.all_target_sites.tsv.gz', sep='\t')
        d_soup = nearest(H.pos.values, tgt.pos.values)
        print(f'{chrom}: {len(H):,} held-out typed variants, {len(donors)} donors', flush=True)
        est = {}
        for arm in ARMS:
            est[arm] = read_ds(f'{SCR}/A07j/impute/{arm}.{chrom}.imputed.bcf', chrom, H, {d: d for d in donors})
        parts = []
        for p, s in POOLS.items():
            bcf = f'{PROJ}/results/A07_onek1k/impute/pool{p}{s}/{chrom}.imputed.bcf'
            parts.append(read_ds(bcf, chrom, H, cmaps[p]).set_index(['pos', 'ref', 'alt']))
        est['latent'] = pd.concat(parts, axis=1).reset_index()
        for arm, E in est.items():
            r2 = r2_rows(T, E[donors].values)
            own = d_soup if arm in ('perfect_soup', 'weak_soup', 'latent') else \
                nearest(H.pos.values, inp.loc[inp.spacing == 'random', 'pos'].values)
            rows.append(pd.DataFrame({'chrom': chrom, 'pos': H.pos, 'maf': H.maf, 'arm': arm, 'r2': r2,
                                      'dist_soup': d_soup, 'dist_own_input': own}))
    R = pd.concat(rows, ignore_index=True)
    R = R[R.r2.notna()]
    R['maf_bin'] = pd.cut(R.maf, MAF_BINS, include_lowest=True).astype(str)
    R['dist_soup_bin'] = pd.cut(R.dist_soup, DIST_BINS, labels=DIST_LABELS).astype(str)
    R['dist_own_bin'] = pd.cut(R.dist_own_input, DIST_BINS, labels=DIST_LABELS).astype(str)
    R.to_csv(f'{OUT}/per_variant.tsv.gz', sep='\t', index=False, float_format='%.4f')
    order = ['perfect_random', 'perfect_soup', 'weak_random', 'weak_soup', 'latent']
    lines = ['per-variant r2 across donors at held-out array-typed SNVs (median / mean):']
    S = []
    for key in [None, 'maf_bin', 'dist_soup_bin', 'dist_own_bin']:
        g = R.groupby(['arm'] + ([key] if key else [])).r2.agg(['size', 'median', 'mean']).reset_index()
        g.insert(0, 'by', key or 'all')
        S.append(g)
        piv = g.pivot_table(index=key or 'by', columns='arm', values='median')[order]
        if key and key.startswith('dist'):
            piv = piv.reindex([x for x in DIST_LABELS if x in piv.index])
        n = g.groupby(key or 'by')['size'].first()
        piv.insert(0, 'n', n.reindex(piv.index))
        lines.append(f'\nmedian r2 by {key or "all"}:')
        lines += ['  ' + l for l in piv.round(3).to_string().split('\n')]
    pd.concat(S).to_csv(f'{OUT}/summary.tsv', sep='\t', index=False, float_format='%.4f')
    open(f'{OUT}/summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
