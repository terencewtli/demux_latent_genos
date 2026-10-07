"""A07h scoring: accuracy of each imputation arm against the OneK1K array genotypes (typed sites, GRCh38), pooled
over the pilot pools so that per-variant r2 is computed across ~70 donors (the scTAPAS metric), not only per donor.

Arms (each: results dir with pool<P>/chrN.imputed.bcf, sample naming):
  soupsites_glimpse   A07f: GLIMPSE2 from souporcell cluster GLs at souporcell's sites   (samples = clusters)
  bam_glimpse_soup    A07h: GLIMPSE2 --bam-list on souporcell-singlet BAMs                (samples = c<k>)
  bam_glimpse_oracle  A07h: GLIMPSE2 --bam-list on GEO-labelled donor BAMs                (samples = OneK1K_<id>)
  bam_quilt_soup      A07h: QUILT2 on souporcell-singlet BAMs
  bam_quilt_oracle    A07h: QUILT2 on GEO-labelled donor BAMs
Cluster -> donor for the soup arms: A07g assign.tsv (Hungarian match on souporcell genotypes).
Strata per variant: cohort MAF (OneK1K array, all 1,098 donors) < 5% / 5-20% / >= 20%, and whether the variant lies
in an expressed gene (gene body of a GENCODE gene with >= MIN_UMI UMIs per 1,000 cells, pooled pilot pseudobulk) -
the sites where scRNA reads exist.
Outputs (results/A07_onek1k/A07h_score/): per_variant.tsv.gz (arm, variant, n, r2, maf, expressed),
  per_donor.tsv, summary.tsv (median / mean per-variant r2 by arm x stratum, mean per-donor r2)
usage: python scripts/A07h_score.py --chroms chr6,chr22
"""
import argparse
import glob
import gzip
import os
import re
import subprocess

import numpy as np
import pandas as pd
import scipy.io

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A07'
RES = f'{PROJ}/results/A07_onek1k'
TRUTH = f'{SCR}/truth/onek1k.typed.b38.vcf.gz'
GTF = '/u/project/cluo/terencew/reference/hg38_igvf/gencode.v43.chr_patch_hapl_scaff.annotation.gtf'
ARMS = {'soupsites_glimpse': f'{RES}/impute', 'bam_glimpse_soup': f'{RES}/impute_bam/soup',
        'bam_glimpse_oracle': f'{RES}/impute_bam/oracle', 'bam_quilt_soup': f'{RES}/impute_quilt/soup',
        'bam_quilt_oracle': f'{RES}/impute_quilt/oracle'}
MIN_UMI = 50
GT = {'0|0': 0, '0|1': 1, '1|0': 1, '1|1': 2, '0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2}


def query(vcf: str, fmt: str, region: str, samples: list = None) -> pd.DataFrame:
    cmd = ['bcftools', 'query', '-r', region, '-f', f'%CHROM:%POS:%REF:%ALT[\\t%{fmt}]\\n']
    if samples:
        cmd[2:2] = ['-s', ','.join(samples)]
    names = samples or subprocess.run(['bcftools', 'query', '-l', vcf], capture_output=True, text=True, check=True).stdout.split()
    out = subprocess.run(cmd + [vcf], capture_output=True, text=True, check=True).stdout
    rows = [l.split('\t') for l in out.splitlines()]
    d = pd.DataFrame(rows, columns=['id'] + names).set_index('id')
    return d


def pools() -> pd.DataFrame:
    return pd.read_csv(f'{PROJ}/txt/onek1k_pilot_pools.txt', sep='\t', header=None, names=['pool', 'gsm', 'k'], dtype=str)


def expressed_genes(chroms: list) -> pd.DataFrame:
    tot, ncell, feats = None, 0, None
    for p in pools().pool:
        d = f'{SCR}/starsolo/pool{p}/Solo.out/Gene/filtered'
        m = scipy.io.mmread(f'{d}/matrix.mtx').tocsr()
        s = np.asarray(m.sum(axis=1)).ravel()
        tot = s if tot is None else tot + s
        ncell += m.shape[1]
        feats = pd.read_csv(f'{d}/features.tsv', sep='\t', header=None)[0].values
    keep = set(feats[tot / ncell * 1000 >= MIN_UMI])
    rows = []
    with open(GTF) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 9)
            if f[2] == 'gene' and f[0] in chroms:
                gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
                if gid in keep:
                    rows.append((f[0], int(f[3]), int(f[4])))
    g = pd.DataFrame(rows, columns=['chrom', 'start', 'end'])
    print(f'expressed genes (>= {MIN_UMI} UMI / 1k cells) on {chroms}: {len(g)}', flush=True)
    return g


def in_genes(ids: pd.Index, g: pd.DataFrame) -> np.ndarray:
    chrom = ids.str.split(':').str[0]
    pos = ids.str.split(':').str[1].astype(int)
    out = np.zeros(len(ids), dtype=bool)
    for c, gg in g.groupby('chrom'):
        m = (chrom == c)
        p = pos[m].values
        hit = np.zeros(len(p), dtype=bool)
        for s, e in zip(gg.start, gg.end):
            hit |= (p >= s) & (p <= e)
        out[m.values] = hit
    return out


def main() -> None:
    a = argparse.ArgumentParser()
    a.add_argument('--chroms', default='chr6,chr22')
    args = a.parse_args()
    chroms = args.chroms.split(',')
    out = f'{RES}/A07h_score'
    os.makedirs(out, exist_ok=True)
    P = pools()
    assign = {p: pd.read_csv(f'{RES}/score/pool{p}/assign.tsv', sep='\t') for p in P.pool}
    donors = sorted(set(pd.concat(assign.values()).donor))
    allsamp = subprocess.run(['bcftools', 'query', '-l', TRUTH], capture_output=True, text=True, check=True).stdout.split()
    G = expressed_genes(chroms)
    per_var, per_don = [], []
    for c in chroms:
        cohort = query(TRUTH, 'GT', c).replace(GT).apply(pd.to_numeric, errors='coerce')
        cohort.columns = allsamp
        af = cohort.mean(axis=1) / 2
        maf = np.minimum(af, 1 - af)
        truth = cohort[donors]
        expr = pd.Series(in_genes(truth.index, G), index=truth.index)
        for arm, d in ARMS.items():
            mats = []
            for p in P.pool:
                f = f'{d}/pool{p}/{c}.imputed.bcf'
                if not os.path.exists(f + '.csi'):
                    continue
                ds = query(f, 'DS', c).apply(pd.to_numeric, errors='coerce')
                A = assign[p].assign(cl=assign[p].cluster.astype(str))
                ren = {}
                for s in ds.columns:
                    key = s if s.startswith(('c', 'OneK1K_')) else f'c{s}'
                    ren[s] = key if key.startswith('OneK1K_') else dict(zip(A.cl, A.donor)).get(key)
                ds = ds.rename(columns=ren)
                ds = ds[[s for s in ds.columns if s in truth.columns]]
                mats.append(ds)
            if not mats:
                print(f'{arm} {c}: no output yet', flush=True)
                continue
            D = pd.concat(mats, axis=1)
            D = D.loc[:, ~D.columns.duplicated()]
            ids = D.index.intersection(truth.index)
            T, X = truth.loc[ids, D.columns], D.loc[ids]
            for s in D.columns:
                ok = T[s].notna() & X[s].notna()
                per_don.append({'arm': arm, 'chrom': c, 'donor': s, 'n': int(ok.sum()),
                                'r2': np.corrcoef(T.loc[ok, s], X.loc[ok, s])[0, 1] ** 2})
            t, x = T.values.astype(float), X.values.astype(float)
            ok = ~np.isnan(t) & ~np.isnan(x)
            n = ok.sum(1)
            tm = np.where(ok, t, 0).sum(1) / n
            xm = np.where(ok, x, 0).sum(1) / n
            tc, xc = np.where(ok, t - tm[:, None], 0), np.where(ok, x - xm[:, None], 0)
            with np.errstate(invalid='ignore', divide='ignore'):
                r = (tc * xc).sum(1) / np.sqrt((tc ** 2).sum(1) * (xc ** 2).sum(1))
            per_var.append(pd.DataFrame({'arm': arm, 'variant': ids, 'n': n, 'r2': r ** 2,
                                         'maf': maf.loc[ids].values, 'expressed': expr.loc[ids].values}))
            print(f'{arm} {c}: {len(D.columns)} donors, {len(ids)} array sites, median per-variant r2 '
                  f'{np.nanmedian(r ** 2):.3f}', flush=True)
    V = pd.concat(per_var)
    V.to_csv(f'{out}/per_variant.tsv.gz', sep='\t', index=False, float_format='%.4g')
    Dn = pd.DataFrame(per_don)
    Dn.to_csv(f'{out}/per_donor.tsv', sep='\t', index=False, float_format='%.4g')
    V['maf_bin'] = pd.cut(V.maf, [0, 0.05, 0.2, 0.5], labels=['<5%', '5-20%', '>=20%'], include_lowest=True)
    rows = []
    for (arm, e, mb), g in V[V.n >= 20].groupby(['arm', 'expressed', 'maf_bin'], observed=True):
        rows.append({'arm': arm, 'expressed_gene': e, 'maf_bin': mb, 'variants': len(g), 'median_r2': g.r2.median(),
                     'mean_r2': g.r2.mean(), 'share_r2_gt_0.8': (g.r2 > 0.8).mean()})
    S = pd.DataFrame(rows)
    S = S.merge(Dn.groupby('arm').r2.mean().rename('mean_donor_r2_all_sites').reset_index(), on='arm')
    S.to_csv(f'{out}/summary.tsv', sep='\t', index=False, float_format='%.3f')
    pd.set_option('display.width', 200)
    print(S.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
