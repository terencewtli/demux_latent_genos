"""Score imputed latent genotypes (souporcell clusters) against true 1000G 30x genotypes.

Works on any imputed file whose samples are souporcell clusters (c0..c{k-1}) and
which has FORMAT/DS: GLIMPSE2 (A04g) or Minimac4 (A04c) output. Clusters are
mapped to donors with A03b's assign.tsv (cluster -> donor, Hungarian on -r);
unassigned clusters are dropped.

Two parts:
Truth = the pool VCF (<tree>/vcf/<pool>.vcf.gz): ambisim simulated reads from
those genotypes, so they are the truth for these pools (user decision 2026-09-29;
do not use vcf_fixed/). They hold MAF >= 5% biallelic SNVs only, and chr1 stops at
190.67 Mb, so untyped rare-variant bins are empty by construction. Panel AF/AC
(3,202 samples) come from the A04a sites file (--panel-sites).

  untyped: sites NOT in the imputation input, same definitions and output format
           as lib/score_imputation.py (MAF bins from the 3,202-sample panel AF,
           aggregate r2 over site x sample, mean per-site r2, private variants).
           score_imputation.py itself calls `bcftools query -m2 -M2 -v snps`,
           which bcftools 1.11 (the module) does not support (-v = --vcf-list),
           and silently scores 0 sites; this module uses -i 'TYPE="snp"' and
           checks the exit status. Its bins and r2() are imported unchanged.
  typed:   the input (souporcell-covered) sites, per donor, A03b style (Pearson r2
           over sites within a donor): souporcell GT hard call and GN expected
           dosage vs imputed DS, on the same sites (those where souporcell calls
           a GT, = A03b's r2_gt site set), and by cluster read depth.

Writes:
  <out>.untyped.sites.tsv.gz  <out>.untyped.bins.tsv
  <out>.donor.tsv   per donor: typed r2 (naive GT / GN / imputed DS) + untyped r2
  <out>.depth.tsv   per donor x depth bin, typed sites
"""
import argparse
import math
import os
import subprocess
import sys
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_imputation import LABELS, MAF_BINS, r2  # noqa: E402

GT_DOSAGE = {'0|0': 0.0, '0|1': 1.0, '1|0': 1.0, '1|1': 2.0, '0/0': 0.0, '0/1': 1.0, '1/0': 1.0, '1/1': 2.0}
DEPTH_BINS = [0, 1, 3, 6, 11, 10 ** 9]
DEPTH_LABELS = ['depth0', 'depth1-2', 'depth3-5', 'depth6-10', 'depth>10']
KEY = ['chrom', 'pos', 'ref', 'alt']


def query_df(cmd: List[str], names: List[str]) -> pd.DataFrame:
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    df = pd.read_csv(p.stdout, sep='\t', header=None, names=names, dtype={'chrom': str}, na_values=['.'])
    if p.wait() != 0:
        raise RuntimeError(f'failed: {" ".join(cmd)}')
    return df


def biallelic_snv(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df.ref.str.len() == 1) & (df.alt.str.len() == 1)]


def read_souporcell(path: str, chrom: str, clusters: List[str]) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.startswith(chrom + '\t'):
                continue
            f = line.rstrip('\n').split('\t')
            if 'BACKGROUND' in f[6] or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            keys = f[8].split(':')
            i_gt, i_gn = keys.index('GT'), keys.index('GN')
            gt, gn = [], []
            for s in f[9:]:
                v = s.split(':')
                gt.append(np.nan if '.' in v[i_gt] else float(v[i_gt].count('1')))
                x = [float(y) for y in v[i_gn].split(',')]
                if len(x) != 3 or any(math.isnan(y) for y in x):
                    gn.append(np.nan)
                else:
                    pr = np.exp(np.array(x) - max(x))
                    pr /= pr.sum()
                    gn.append(2 * pr[1] + pr[2])  # order 0/0, 1/1, 0/1
            rows.append([f[0], int(f[1]), f[3], f[4]] + gt + gn)
    k = (len(rows[0]) - 4) // 2 if rows else 0
    cols = KEY + [f'gt_c{i}' for i in range(k)] + [f'gn_c{i}' for i in range(k)]
    return pd.DataFrame(rows, columns=cols)


def rowwise_r2(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    num = (xc * yc).sum(axis=1) ** 2
    den = (xc ** 2).sum(axis=1) * (yc ** 2).sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0, num / den, np.nan)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--dose', required=True, help='imputed VCF/BCF with FORMAT/DS, samples = clusters')
    p.add_argument('--typed', required=True, help='imputation input (sites = typed); FORMAT/DP or AD for depth')
    p.add_argument('--soup', required=True, help='souporcell cluster_genotypes.vcf (naive GT / GN)')
    p.add_argument('--assign', required=True, help='A03b assign.tsv (cluster, donor)')
    p.add_argument('--truth', required=True, help='pool truth VCF (<tree>/vcf/<pool>.vcf.gz), samples = donors')
    p.add_argument('--panel-sites', required=True, help='A04a 1000G_30x.chrN.sites.vcf.gz (INFO/AF, AC over 3,202)')
    p.add_argument('--region', required=True)
    p.add_argument('--out', required=True, help='output prefix')
    a = p.parse_args()

    A = pd.read_csv(a.assign, sep='\t')
    c2d: Dict[str, str] = dict(zip(A.cluster, A.donor))
    dose_samples = subprocess.run(['bcftools', 'query', '-l', a.dose], capture_output=True, text=True,
                                  check=True).stdout.split()
    clusters = [c for c in dose_samples if c in c2d]
    donors = [c2d[c] for c in clusters]
    print(f'clusters scored: {len(clusters)} of {len(dose_samples)}; unassigned: '
          f'{[c for c in dose_samples if c not in c2d]}', flush=True)

    fmt = '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n'
    dose = biallelic_snv(query_df(['bcftools', 'query', '-r', a.region, '-s', ','.join(clusters),
                                   '-i', 'TYPE="snp"', '-f', fmt, a.dose], KEY + clusters))
    typed = query_df(['bcftools', 'query', '-r', a.region, '-s', ','.join(clusters), '-f',
                      '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n', a.typed], KEY + clusters)
    for c in clusters:
        typed[c] = typed[c].astype(str).str.split(',').apply(lambda v: sum(int(x) for x in v if x != '.'))
    truth = biallelic_snv(query_df(['bcftools', 'query', '-r', a.region, '-s', ','.join(donors), '-i', 'TYPE="snp"',
                                    '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n', a.truth], KEY + donors))
    for d in donors:
        truth[d] = truth[d].map(GT_DOSAGE)
    panel = biallelic_snv(query_df(['bcftools', 'query', '-r', a.region, '-i', 'TYPE="snp"', '-f',
                                    '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n', a.panel_sites],
                                   KEY + ['af', 'ac_panel'])).drop_duplicates(KEY)
    truth = truth.merge(panel, on=KEY)
    print(f'dose SNVs: {len(dose)}  typed sites: {len(typed)}  truth SNVs (with panel AF): {len(truth)}', flush=True)

    m = dose.merge(truth, on=KEY)  # cluster (c*) and donor columns never collide
    ids = m.chrom + ':' + m.pos.astype(str) + ':' + m.ref + ':' + m.alt
    typed_ids = set(typed.chrom + ':' + typed.pos.astype(str) + ':' + typed.ref + ':' + typed.alt)
    is_typed = ids.isin(typed_ids).values
    X = m[donors].values.astype(float)                          # truth dosage
    Y = m[clusters].values.astype(float)                        # imputed DS

    ### untyped: score_imputation.py format
    u = ~is_typed
    Xu, Yu = X[u], Y[u]
    n_dose_untyped = int((~(dose.chrom + ':' + dose.pos.astype(str) + ':' + dose.ref + ':' + dose.alt).isin(typed_ids)).sum())
    print(f'untyped SNVs imputed: {n_dose_untyped}; with truth (scored): {int(u.sum())}', flush=True)
    ac_s = np.nansum(Xu, axis=1).astype(int)
    df = pd.DataFrame({'id': ids.values[u], 'af': m.af.values[u], 'ac_panel': m.ac_panel.values[u].astype(int),
                       'ac_samples': ac_s})
    df['private'] = ((df.ac_samples > 0) & (df.ac_samples == df.ac_panel)).astype(int)
    df['r2'] = rowwise_r2(Xu, Yu)
    df['bin'] = pd.cut(np.minimum(df.af, 1 - df.af), MAF_BINS, labels=LABELS, right=False)
    df.to_csv(f'{a.out}.untyped.sites.tsv.gz', sep='\t', index=False)
    out = []
    groups = [('maf_' + str(b), (df.bin == b).values) for b in LABELS]
    groups += [('carried_by_samples', (df.ac_samples > 0).values), ('private_to_samples', (df.private == 1).values),
               ('private_singleton', ((df.private == 1) & (df.ac_panel == 1)).values),
               ('all_untyped', np.ones(len(df), bool))]
    for name, sel in groups:
        if sel.sum() == 0:
            out.append((name, 0, np.nan, np.nan))
            continue
        out.append((name, int(sel.sum()), r2(Xu[sel].ravel(), Yu[sel].ravel()), df.r2[sel].mean()))
    bins = pd.DataFrame(out, columns=['group', 'n_sites', 'aggregate_r2', 'mean_site_r2'])
    bins.to_csv(f'{a.out}.untyped.bins.tsv', sep='\t', index=False)
    print(bins.to_string(index=False), flush=True)

    ### typed: per donor, A03b style
    soup = read_souporcell(a.soup, a.region.split(':')[0], clusters)
    t = m[is_typed][KEY + donors].copy()
    t[[f'ds_{c}' for c in clusters]] = Y[is_typed]
    t = t.merge(soup, on=KEY).merge(typed.rename(columns={c: f'dp_{c}' for c in clusters}), on=KEY)
    maf_u = np.minimum(df.af.values, 1 - df.af.values)
    rows, drows = [], []
    for j, (c, d) in enumerate(zip(clusters, donors)):
        tr, ds = t[d].values, t[f'ds_{c}'].values
        gt, gn, dp = t[f'gt_{c}'].values, t[f'gn_{c}'].values, t[f'dp_{c}'].values
        called = ~np.isnan(gt)
        ucom = maf_u >= 0.05
        rows.append((c, d, len(t), int(called.sum()), r2(gt, tr), r2(gn[called], tr[called]), r2(ds[called], tr[called]),
                     r2(ds, tr), float(np.mean(gt[called] == tr[called])), float(np.mean(np.round(ds[called]) == tr[called])),
                     int(u.sum()), r2(Yu[:, j], Xu[:, j]), int(ucom.sum()), r2(Yu[ucom, j], Xu[ucom, j])))
        dbin = np.array(pd.cut(dp, DEPTH_BINS, labels=DEPTH_LABELS, right=False).astype(str))
        for b in DEPTH_LABELS:
            s = dbin == b
            sc = s & called
            drows.append((c, d, b, int(s.sum()), int(sc.sum()), r2(gt[sc], tr[sc]), r2(ds[sc], tr[sc]), r2(ds[s], tr[s])))
    D = pd.DataFrame(rows, columns=['cluster', 'donor', 'n_typed', 'n_called', 'r2_naive_gt', 'r2_naive_gn',
                                    'r2_imputed_called', 'r2_imputed_all_typed', 'conc_naive_gt', 'conc_imputed_called',
                                    'n_untyped', 'r2_untyped', 'n_untyped_maf5', 'r2_untyped_maf5'])
    S = pd.DataFrame(drows, columns=['cluster', 'donor', 'depth', 'n_typed', 'n_called', 'r2_naive_gt',
                                     'r2_imputed_called', 'r2_imputed_all'])
    D.to_csv(f'{a.out}.donor.tsv', sep='\t', index=False)
    S.to_csv(f'{a.out}.depth.tsv', sep='\t', index=False)
    print(D.to_string(index=False))
    print(S.groupby('depth')[['n_called', 'r2_naive_gt', 'r2_imputed_called', 'r2_imputed_all']].mean().to_string())
    print('mean over donors: ' + ' '.join(f'{k}={D[k].mean():.4f}' for k in
          ['r2_naive_gt', 'r2_naive_gn', 'r2_imputed_called', 'r2_imputed_all_typed', 'r2_untyped', 'r2_untyped_maf5']))


if __name__ == '__main__':
    main()
