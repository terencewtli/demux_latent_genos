"""Naive dosage r2 of souporcell latent genotypes against the true pool genotypes.

For every finished A03a row of txt/latent_genos_souporcell_tasks.txt:
  1. read cluster_genotypes.vcf (BACKGROUND loci and non-SNVs dropped)
       gt dosage: GT hard call (0/1/2, ./. -> NaN)
       gl dosage: expected dosage from GN (natural-log posteriors, int-truncated,
                  order [0/0, 1/1, 0/1]); NaN where the cluster has 0 reads
  2. match to the pool truth VCF (<tree>/vcf/<pool>.vcf.gz, 1000G 30x common
     biallelic SNVs, AC=0 sites included) on CHROM/POS; REF/ALT exact or swapped
     (swapped -> 2 - dosage); anything else dropped
  3. cluster -> donor: Pearson r of gt dosage vs truth dosage over each cluster's
     called sites, Hungarian assignment on -r; margin = r(assigned) - best r to
     any other donor
  4. per matched pair: r2 of gt and gl dosage; stratified by EUR MAF (1000G
     AF_EUR), pool MAF (from the truth genotypes), and cluster read depth (AO+RO)

Outputs, in <outdir>/<tree>__<pool>__<mod>/ (skipped if assign.tsv exists):
  assign.tsv  r_matrix.tsv  r2_donor.tsv  r2_strata.tsv  site_counts.tsv
Then summary tables over all finished rows in <outdir>/:
  summary_donor.n<rows>.tsv  summary_pool.n<rows>.tsv  summary_strata.n<rows>.tsv

Usage (allcools env, bcftools on PATH):
  python A03b_souporcell_dosage_r2.py [--tasks ...] [--outdir ...]
"""
import argparse
import math
import os
import subprocess
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

SAMPLE = '20220928-IGVF-D0'
PROJDIR = '/u/project/cluo/terencew/claude/project_ideas/pool_design'
OUTDIR = '/u/project/cluo/terencew/claude/project_ideas/latent_genos/results/A03b_souporcell_dosage_r2'
EUR_BINS = [0, 0.05, 0.2, 0.5000001]
EUR_LABELS = ['eurMAF<5%', 'eurMAF5-20%', 'eurMAF>20%']
POOL_LABELS = ['poolMAF=0', 'poolMAF<10%', 'poolMAF10-25%', 'poolMAF>25%']
DEPTH_BINS = [1, 3, 6, 11, 10 ** 9]
DEPTH_LABELS = ['depth1-2', 'depth3-5', 'depth6-10', 'depth>10']


def read_souporcell(path: str) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    keys, gt, gl, dp = [], [], [], []
    n_bg = n_nonsnv = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if 'BACKGROUND' in f[6]:
                n_bg += 1
                continue
            if len(f[3]) != 1 or len(f[4]) != 1 or f[4] not in 'ACGT':
                n_nonsnv += 1
                continue
            fmt = f[8].split(':')
            i_gt, i_ao, i_ro, i_gn = (fmt.index(x) for x in ('GT', 'AO', 'RO', 'GN'))
            g_row, l_row, d_row = [], [], []
            for s in f[9:]:
                v = s.split(':')
                g = v[i_gt]
                g_row.append(np.nan if '.' in g else float(g.count('1')))
                d = int(v[i_ao]) + int(v[i_ro])
                d_row.append(d)
                gn = [float(x) for x in v[i_gn].split(',')]
                if d == 0 or len(gn) != 3 or any(math.isnan(x) for x in gn):
                    l_row.append(np.nan)
                else:
                    p = np.exp(np.array(gn) - max(gn))
                    p /= p.sum()
                    l_row.append(2 * p[1] + p[2])
            keys.append((f[0].replace('chr', ''), int(f[1]), f[3], f[4]))
            gt.append(g_row)
            gl.append(l_row)
            dp.append(d_row)
    print(f'  souporcell: kept={len(keys)} background={n_bg} non_snv={n_nonsnv}', flush=True)
    sites = pd.DataFrame(keys, columns=['chrom', 'pos', 'ref', 'alt'])
    return sites, np.array(gt), np.array(gl), np.array(dp)


def read_truth(vcf: str, sites: pd.DataFrame, tmp: str) -> Tuple[List[str], pd.DataFrame, np.ndarray]:
    donors = subprocess.run(['bcftools', 'query', '-l', vcf], capture_output=True, text=True,
                            check=True).stdout.split()
    sites[['chrom', 'pos']].assign(chrom='chr' + sites.chrom).drop_duplicates().to_csv(
        tmp, sep='\t', header=False, index=False)
    # truth VCF is already biallelic SNVs only (bcftools 1.11 query has no -m/-v)
    cmd = ['bcftools', 'query', '-T', tmp,
           '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_EUR[\t%GT]\n', vcf]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert p.stdout is not None
    rows, dos = [], []
    for line in p.stdout:
        r = line.rstrip('\n').split('\t')
        if len(r[2]) != 1 or len(r[3]) != 1:
            continue
        rows.append((r[0].replace('chr', ''), int(r[1]), r[2], r[3], float(r[4])))
        dos.append([np.nan if '.' in g else float(g.count('1')) for g in r[5:]])
    if p.wait() != 0:
        raise RuntimeError(f'bcftools query failed on {vcf}')
    os.remove(tmp)
    return donors, pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt', 'af_eur']), np.array(dos)


def pearson(x: np.ndarray, y: np.ndarray) -> Tuple[float, int]:
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if len(x) < 3 or x.std() == 0 or y.std() == 0:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(x, y)[0, 1]), int(ok.sum())


def score_task(tree: str, pool: str, mod: str, outdir: str) -> None:
    soup = f'{PROJDIR}/{tree}/{pool}/demux/souporcell/{mod}/{SAMPLE}/cluster_genotypes.vcf'
    truth_vcf = f'{PROJDIR}/{tree}/vcf/{pool}.vcf.gz'
    os.makedirs(outdir, exist_ok=True)

    sites, gt, gl, dp = read_souporcell(soup)
    donors, truth, tdos = read_truth(truth_vcf, sites, f'{outdir}/.sites.tmp')
    print(f'  truth records at souporcell positions: {len(truth)} donors={len(donors)}', flush=True)

    t = truth.reset_index().rename(columns={'index': 'ti'})
    s = sites.reset_index().rename(columns={'index': 'si'})
    m = s.merge(t, on=['chrom', 'pos'], suffixes=('', '_t'))
    same = (m.ref == m.ref_t) & (m.alt == m.alt_t)
    swap = (m.ref == m.alt_t) & (m.alt == m.ref_t)
    m = m[same | swap].copy()
    m['swap'] = swap[same | swap].values
    m = m.drop_duplicates('si')
    counts = {'souporcell_snv': len(sites), 'matched_exact': int((~m.swap).sum()),
              'matched_swapped': int(m.swap.sum()), 'unmatched': len(sites) - len(m)}

    si, ti = m.si.values, m.ti.values
    sw = m.swap.values[:, None]
    G = np.where(sw, 2 - gt[si], gt[si])
    L = np.where(sw, 2 - gl[si], gl[si])
    D = dp[si]
    T = tdos[ti]
    af_eur = truth.af_eur.values[ti]
    k = G.shape[1]

    R = np.full((k, len(donors)), np.nan)
    for c in range(k):
        for j in range(len(donors)):
            R[c, j] = pearson(G[:, c], T[:, j])[0]
    # clusters with no called GTs (souporcell leaves a cluster empty when two
    # donors collapse into one) cannot be assigned; their donors stay unmatched
    valid = np.where(~np.all(np.isnan(R), axis=1))[0]
    counts['empty_clusters'] = ','.join(f'c{c}' for c in range(k) if c not in valid) or 'none'
    vr, cols = linear_sum_assignment(-np.nan_to_num(R[valid], nan=-1))
    rows = valid[vr]
    assign = []
    for c, j in zip(rows, cols):
        other = np.delete(R[c], j)
        other_c = np.delete(R[:, j], c)
        assign.append((f'c{c}', donors[j], R[c, j], np.nanmax(other), R[c, j] - np.nanmax(other),
                       R[c, j] - np.nanmax(other_c), int(np.nanargmax(R[c]) == j),
                       donors[int(np.nanargmax(R[c]))]))
    A = pd.DataFrame(assign, columns=['cluster', 'donor', 'r', 'r_next_donor', 'margin_donor',
                                      'margin_cluster', 'argmax_agrees', 'argmax_donor'])
    A['ambiguous'] = ((A[['margin_donor', 'margin_cluster']].min(axis=1) < 0.1) | (A.argmax_agrees == 0)).astype(int)
    counts['unmatched_donors'] = ','.join(d for i, d in enumerate(donors) if i not in set(cols)) or 'none'

    pool_af = np.nanmean(T, axis=1) / 2
    pool_maf = np.minimum(pool_af, 1 - pool_af)
    pool_bin = np.select([pool_maf == 0, pool_maf < 0.1, pool_maf < 0.25],
                         POOL_LABELS[:3], POOL_LABELS[3])
    eur_bin = pd.cut(np.minimum(af_eur, 1 - af_eur), EUR_BINS, labels=EUR_LABELS, right=False).astype(str)

    per_donor, strata = [], []
    for c, j in zip(rows, cols):
        g, l, d, tr = G[:, c], L[:, c], D[:, c], T[:, j]
        r_gt, n_gt = pearson(g, tr)
        r_gl, n_gl = pearson(l, tr)
        ok = ~np.isnan(g) & ~np.isnan(tr)
        conc = float((g[ok] == tr[ok]).mean())
        per_donor.append((f'c{c}', donors[j], n_gt, r_gt ** 2, n_gl, r_gl ** 2, conc, float(np.median(d[d > 0]))))
        depth_bin = np.array(pd.cut(d, DEPTH_BINS, labels=DEPTH_LABELS, right=False).astype(str))
        for kind, lab, labels in (('eur_maf', eur_bin, EUR_LABELS), ('pool_maf', pool_bin, POOL_LABELS),
                                  ('depth', depth_bin, DEPTH_LABELS)):
            for b in labels:
                sel = lab == b
                rg, ng = pearson(g[sel], tr[sel])
                rl, nl = pearson(l[sel], tr[sel])
                okb = sel & ok
                cb = float((g[okb] == tr[okb]).mean()) if okb.any() else np.nan
                strata.append((f'c{c}', donors[j], kind, b, ng, rg ** 2, nl, rl ** 2, cb))
    Dn = pd.DataFrame(per_donor, columns=['cluster', 'donor', 'n_sites_gt', 'r2_gt', 'n_sites_gl', 'r2_gl',
                                          'gt_concordance', 'median_depth'])
    S = pd.DataFrame(strata, columns=['cluster', 'donor', 'stratum', 'bin', 'n_sites_gt', 'r2_gt',
                                      'n_sites_gl', 'r2_gl', 'gt_concordance'])

    pd.DataFrame(R, index=[f'c{c}' for c in range(k)], columns=donors).to_csv(f'{outdir}/r_matrix.tsv', sep='\t')
    Dn.to_csv(f'{outdir}/r2_donor.tsv', sep='\t', index=False)
    S.to_csv(f'{outdir}/r2_strata.tsv', sep='\t', index=False)
    pd.Series(counts).to_csv(f'{outdir}/site_counts.tsv', sep='\t', header=False)
    A.to_csv(f'{outdir}/assign.tsv', sep='\t', index=False)  # written last = done marker
    print(A.to_string(index=False))
    print(Dn.to_string(index=False))
    print(f'  {counts}', flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--tasks', default=f'{PROJDIR}/txt/latent_genos_souporcell_tasks.txt')
    p.add_argument('--outdir', default=OUTDIR)
    a = p.parse_args()

    done: List[Tuple[str, str, str, str]] = []
    for line in open(a.tasks):
        if not line.strip():
            continue
        tree, pool, _, mod = line.split()
        sdir = f'{PROJDIR}/{tree}/{pool}/demux/souporcell/{mod}/{SAMPLE}'
        tag = f'{tree}__{pool}__{mod}'
        if not (os.path.exists(f'{sdir}/consensus.done') and os.path.getsize(f'{sdir}/cluster_genotypes.vcf') > 0
                if os.path.exists(f'{sdir}/cluster_genotypes.vcf') else False):
            print(f'{tag}: souporcell not finished, skipping', flush=True)
            continue
        odir = f'{a.outdir}/{tag}'
        if os.path.exists(f'{odir}/assign.tsv'):
            print(f'{tag}: already scored, reusing', flush=True)
        else:
            print(f'{tag}: scoring', flush=True)
            score_task(tree, pool, mod, odir)
        done.append((tree, pool, mod, odir))

    if not done:
        return
    donor, strata = [], []
    for tree, pool, mod, odir in done:
        meta = {'tree': tree, 'pool': pool, 'mod': mod}
        donor.append(pd.read_csv(f'{odir}/r2_donor.tsv', sep='\t').merge(
            pd.read_csv(f'{odir}/assign.tsv', sep='\t')[['cluster', 'r', 'margin_donor', 'margin_cluster', 'ambiguous']],
            on='cluster').assign(**meta))
        strata.append(pd.read_csv(f'{odir}/r2_strata.tsv', sep='\t').assign(**meta))
    donor_df = pd.concat(donor)
    strata_df = pd.concat(strata)
    n = len(done)
    pool_df = donor_df.groupby(['tree', 'pool', 'mod']).agg(
        n_donors=('donor', 'size'), n_ambiguous=('ambiguous', 'sum'), min_margin=('margin_donor', 'min'),
        median_sites_gt=('n_sites_gt', 'median'), mean_r2_gt=('r2_gt', 'mean'), mean_r2_gl=('r2_gl', 'mean'),
        mean_concordance=('gt_concordance', 'mean')).reset_index()
    strata_sum = strata_df.groupby(['mod', 'stratum', 'bin']).agg(
        median_sites=('n_sites_gt', 'median'), mean_r2_gt=('r2_gt', 'mean'), mean_r2_gl=('r2_gl', 'mean'),
        mean_concordance=('gt_concordance', 'mean')).reset_index()
    donor_df.to_csv(f'{a.outdir}/summary_donor.n{n}.tsv', sep='\t', index=False)
    pool_df.to_csv(f'{a.outdir}/summary_pool.n{n}.tsv', sep='\t', index=False)
    strata_sum.to_csv(f'{a.outdir}/summary_strata.n{n}.tsv', sep='\t', index=False)
    print(pool_df.to_string(index=False))
    print(strata_sum.to_string(index=False))
    print('overall mean r2 across donors: gt={:.4f} gl={:.4f}'.format(donor_df.r2_gt.mean(), donor_df.r2_gl.mean()))


if __name__ == '__main__':
    main()
