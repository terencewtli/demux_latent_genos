"""OneK1K pilot scoring for one pool: souporcell clusters -> OneK1K donors, then latent / imputed genotype
accuracy against the array-typed truth (A07c), using the simulation scorer (truth = array-typed GRCh38 sites) per chromosome.

Steps:
  1. assign: Pearson r between each cluster's souporcell hard-call dosage and every one of the 1,098
     OneK1K donors at shared typed sites (genome-wide); Hungarian assignment on -r (scipy), plus
     argmax / runner-up / margin columns in the A03b assign.tsv format.
  2. GEO check: the pool's GEO donor list uses IDs "<a>_<b>"; report whether assigned donors are
     OneK1K_<a> or OneK1K_<b> (or neither = assignment problem).
  3. cell check: souporcell singlet cluster -> donor vs GEO per-cell donor labels (Individual_Barcodes).
  4. truth.pool.vcf.gz = truth restricted to the assigned donors; score each chromosome with
     pool_design/.../lib/score_latent_imputation.py (same as the simulations, A04g step 4).
     Its "typed" = souporcell-covered sites that are on the array; "untyped" = other array sites.
  5. genome-wide per-donor summary: means of the per-chromosome r2 weighted by site counts
     (n_called for typed, n_untyped for untyped).

usage: python scripts/A07g_onek1k_score.py --pool 1
"""
import argparse
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
LIB = '/u/project/cluo/terencew/claude/project_ideas/pool_design/scripts/latent_genos/lib'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A07'
PANEL = '/u/project/cluo/terencew/reference/topmed/local_1000G_30x'
TRUTH = f'{SCR}/truth/onek1k.typed.b38.vcf.gz'
GEO = f'{PROJ}/reference/onek1k_geo'
GT_DOSAGE = {'0|0': 0.0, '0|1': 1.0, '1|0': 1.0, '1|1': 2.0, '0/0': 0.0, '0/1': 1.0, '1/0': 1.0, '1/1': 2.0}


def log(msg: str) -> None:
    print(msg, flush=True)


def run(cmd: list) -> None:
    log('$ ' + ' '.join(cmd))
    subprocess.run(cmd, check=True)


def read_soup_gt(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if 'BACKGROUND' in f[6] or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            i_gt = f[8].split(':').index('GT')
            gts = [s.split(':')[i_gt] for s in f[9:]]
            rows.append([f[0], int(f[1]), f[3], f[4]] + [np.nan if '.' in g else float(g.count('1')) for g in gts])
    k = len(rows[0]) - 4
    return pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt'] + [f'c{i}' for i in range(k)])


def truth_at(sites: pd.DataFrame, outdir: str) -> pd.DataFrame:
    reg = f'{outdir}/soup_sites.tsv'
    sites[['chrom', 'pos']].drop_duplicates().to_csv(reg, sep='\t', header=False, index=False)
    samples = subprocess.run(['bcftools', 'query', '-l', TRUTH], check=True, capture_output=True, text=True).stdout.split()
    p = subprocess.Popen(['bcftools', 'query', '-R', reg, '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n', TRUTH],
                         stdout=subprocess.PIPE, text=True)
    df = pd.read_csv(p.stdout, sep='\t', header=None, names=['chrom', 'pos', 'ref', 'alt'] + samples, dtype={'chrom': str})
    if p.wait() != 0:
        raise RuntimeError('bcftools query failed')
    for s in samples:
        df[s] = df[s].map(GT_DOSAGE)
    return df


def assign(soup: pd.DataFrame, truth: pd.DataFrame, outdir: str) -> pd.DataFrame:
    clusters = [c for c in soup.columns if c.startswith('c') and c[1:].isdigit()]
    donors = [c for c in truth.columns if c.startswith('OneK1K_')]
    m = soup.merge(truth, on=['chrom', 'pos', 'ref', 'alt'])
    log(f'assign: {len(soup)} souporcell SNVs, {len(m)} on the array (allele-matched)')
    D = m[donors].values.astype(float)
    R = np.full((len(clusters), len(donors)), np.nan)
    for i, c in enumerate(clusters):
        x = m[c].values
        ok = ~np.isnan(x)
        X = x[ok] - x[ok].mean()
        Y = D[ok] - np.nanmean(D[ok], axis=0)
        Y = np.nan_to_num(Y)
        den = np.sqrt((X ** 2).sum() * (Y ** 2).sum(axis=0))
        with np.errstate(invalid='ignore', divide='ignore'):
            R[i] = (X @ Y) / den
    Rz = np.nan_to_num(R, nan=-1.0)
    ri, ci = linear_sum_assignment(-Rz)
    out = []
    for i, j in zip(ri, ci):
        row = np.sort(Rz[i])[::-1]
        col = np.sort(Rz[:, j])[::-1]
        am = int(np.argmax(Rz[i]))
        out.append({'cluster': clusters[i], 'donor': donors[j], 'r': Rz[i, j],
                    'r_next_donor': row[1], 'margin_donor': Rz[i, j] - row[1],
                    'margin_cluster': Rz[i, j] - (col[1] if len(col) > 1 else np.nan),
                    'argmax_agrees': int(am == j), 'argmax_donor': donors[am],
                    'n_sites': int((~np.isnan(m[clusters[i]].values)).sum())})
    A = pd.DataFrame(out)
    A['ambiguous'] = ((A.margin_donor < 0.1) | (A.argmax_agrees == 0)).astype(int)
    A.to_csv(f'{outdir}/assign.tsv', sep='\t', index=False)
    log(A.to_string(index=False))
    return A


def geo_checks(A: pd.DataFrame, gsm: str, soup_dir: str, outdir: str) -> None:
    listed = pd.read_csv(glob.glob(f'{GEO}/{gsm}_*GenotypeSamples.txt.gz')[0], header=None)[0].astype(str).tolist()
    first = {f'OneK1K_{x.split("_")[0]}': x for x in listed}
    second = {f'OneK1K_{x.split("_")[1]}': x for x in listed}
    n1, n2 = A.donor.isin(first).sum(), A.donor.isin(second).sum()
    log(f'GEO donor list ({len(listed)}): assigned donors matching OneK1K_<a>: {n1}, OneK1K_<b>: {n2}, of {len(A)} clusters')
    conv = first if n1 >= n2 else second
    A = A.assign(geo_id=A.donor.map(conv))
    cells = pd.read_csv(f'{soup_dir}/clusters.tsv', sep='\t')
    cells = cells[cells.status == 'singlet'].copy()
    cells['geo_id'] = cells.assignment.astype(str).map(dict(zip(A.cluster.str[1:], A.geo_id)))
    geo = pd.read_csv(glob.glob(f'{GEO}/{gsm}_*Individual_Barcodes.csv.gz')[0])
    geo['barcode'] = geo.Barcode.str.replace(r'-1$', '', regex=True)
    cells['barcode'] = cells.barcode.str.replace(r'-1$', '', regex=True)
    j = cells.merge(geo, on='barcode')
    conc = (j.geo_id == j['Individual ID'].astype(str)).mean() if len(j) else np.nan
    log(f'cells: souporcell singlets {len(cells)}, GEO-labelled cells {len(geo)}, shared {len(j)}, donor concordance {conc:.4f}')
    pd.DataFrame([{'geo_listed': len(listed), 'match_first_id': n1, 'match_second_id': n2,
                   'soup_singlets': len(cells), 'geo_cells': len(geo), 'shared_cells': len(j),
                   'cell_concordance': conc}]).to_csv(f'{outdir}/geo_checks.tsv', sep='\t', index=False)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--pool', required=True)
    a = p.parse_args()
    pools = pd.read_csv(f'{PROJ}/txt/onek1k_pilot_pools.txt', sep='\t', header=None, names=['pool', 'gsm', 'k'], dtype=str)
    gsm = pools.loc[pools.pool == a.pool, 'gsm'].item()
    soup_dir = f'{SCR}/souporcell/pool{a.pool}'
    imp = f'{PROJ}/results/A07_onek1k/impute/pool{a.pool}'
    outdir = f'{PROJ}/results/A07_onek1k/score/pool{a.pool}'
    os.makedirs(outdir, exist_ok=True)

    soup = read_soup_gt(f'{soup_dir}/cluster_genotypes.vcf')
    A = assign(soup, truth_at(soup, outdir), outdir)
    geo_checks(A, gsm, soup_dir, outdir)

    pool_truth = f'{outdir}/truth.pool.vcf.gz'
    run(['bcftools', 'view', '-s', ','.join(A.donor), '-Oz', '-o', pool_truth, TRUTH])
    run(['bcftools', 'index', '-t', '-f', pool_truth])

    donor_tabs = []
    for c in range(1, 23):
        chrom = f'chr{c}'
        dose = f'{imp}/{chrom}.imputed.bcf'
        if not os.path.exists(dose + '.csi'):
            log(f'{chrom}: no imputed output, skipped')
            continue
        prefix = f'{outdir}/{chrom}.score'
        if not os.path.exists(prefix + '.donor.tsv'):
            run([sys.executable, f'{LIB}/score_latent_imputation.py', '--dose', dose, '--typed', f'{imp}/{chrom}.target.bcf',
                 '--soup', f'{soup_dir}/cluster_genotypes.vcf', '--assign', f'{outdir}/assign.tsv', '--truth', pool_truth,
                 '--panel-sites', f'{PANEL}/1000G_30x.{chrom}.sites.vcf.gz', '--region', chrom, '--out', prefix])
        donor_tabs.append(pd.read_csv(prefix + '.donor.tsv', sep='\t').assign(chrom=chrom))
    D = pd.concat(donor_tabs)

    def wmean(g: pd.DataFrame, col: str, w: str) -> float:
        ok = g[col].notna() & (g[w] > 0)
        return np.average(g.loc[ok, col], weights=g.loc[ok, w]) if ok.any() else np.nan

    S = D.groupby(['cluster', 'donor']).apply(lambda g: pd.Series({
        'n_called': g.n_called.sum(), 'r2_naive_gt': wmean(g, 'r2_naive_gt', 'n_called'),
        'r2_imputed_called': wmean(g, 'r2_imputed_called', 'n_called'),
        'n_untyped': g.n_untyped.sum(), 'r2_untyped': wmean(g, 'r2_untyped', 'n_untyped'),
        'n_chrom': len(g)})).reset_index()
    S.insert(0, 'pool', a.pool)
    S.to_csv(f'{outdir}/summary_donor.tsv', sep='\t', index=False)
    log(S.to_string(index=False))
    log(f'pool {a.pool} mean over donors: naive {S.r2_naive_gt.mean():.3f}  imputed typed {S.r2_imputed_called.mean():.3f}  '
        f'imputed untyped(array) {S.r2_untyped.mean():.3f}')


if __name__ == '__main__':
    main()
