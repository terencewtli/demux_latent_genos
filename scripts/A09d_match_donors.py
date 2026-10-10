#!/usr/bin/env python
"""A09d: map souporcell clusters to donors (WGS genotype concordance) and score against ambimux singlets.

Per run (pool x modality):
  1. cluster x donor genotype concordance at covered sites with a confident WGS call (FT=PASS, GQ>=20)
     -> one-to-one cluster->donor map (Hungarian on discordance).
  2. souporcell singlet calls (mapped to donors) vs ambimux joint singlets.
Per pool: GEX vs ATAC per-barcode agreement after mapping.

Outputs: results/A09_multiome/match/{cluster_donor_concordance.tsv, cluster_map.tsv, barcode_calls.tsv,
         run_summary.tsv, summary.txt}
"""
import os
import subprocess
from io import StringIO

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A09/souporcell_v2'
MULTI = '/u/project/cluo/terencew/igvf/2023_YR2/multiome'
WGS = f'{MULTI}/vcf/wgs/pgp_filt1.rm_missing.reheader.vcf.gz'
AMBI = f'{MULTI}/csv/demux/ambimux/wgs/ambimux_joint_sings.csv'
TIMEMAP = f'{MULTI}/txt/demux_time_map.txt'
BCFTOOLS = '/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools'
OUT = f'{PROJ}/results/A09_multiome/match'
DONORS = ['C29', 'C37', 'C38', 'C39']
MIN_GQ = 20


def read_cluster_vcf(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            gts = []
            for s in f[9:13]:
                g = s.split(':', 1)[0]
                gts.append(np.nan if '.' in g else sum(int(a) for a in g.replace('|', '/').split('/')))
            rows.append((f[0], int(f[1]), f[3], f[4], *gts))
    return pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt', 'c0', 'c1', 'c2', 'c3'])


def wgs_at(sites: pd.DataFrame, tmp: str) -> pd.DataFrame:
    sites[['chrom', 'pos']].drop_duplicates().sort_values(['chrom', 'pos']).to_csv(
        tmp, sep='\t', header=False, index=False)
    fmt = '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%GQ:%FT]\n'
    txt = subprocess.run([BCFTOOLS, 'query', '-R', tmp, '-s', ','.join(DONORS), '-f', fmt, WGS],
                         check=True, capture_output=True, text=True).stdout
    w = pd.read_csv(StringIO(txt), sep='\t', header=None, names=['chrom', 'pos', 'ref', 'alt', *DONORS])
    w = w[~w.alt.str.contains(',')]
    for d in DONORS:
        p = w[d].str.split(':', expand=True)
        gt = p[0].str.replace('|', '/', regex=False)
        ok = (p[2] == 'PASS') & (pd.to_numeric(p[1], errors='coerce') >= MIN_GQ) & ~gt.str.contains(r'\.')
        w[d] = np.where(ok, gt.str.count('1'), np.nan)
    return w.drop_duplicates(['chrom', 'pos', 'ref', 'alt'])


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    runs = pd.read_csv(f'{PROJ}/txt/A09_runs.txt', sep=r'\s+', header=None, names=['pool', 'mod'])
    ambi = pd.read_csv(AMBI, sep='\t', usecols=['sample', 'barcode', 'assignment', 'time'])
    tmap = pd.read_csv(TIMEMAP, sep='\t')

    print('loading cluster VCFs', flush=True)
    cv = {f'{p}_{m}': read_cluster_vcf(f'{SCR}/{p}_{m}/cluster_genotypes.vcf') for p, m in runs.values}
    allsites = pd.concat([v[['chrom', 'pos']] for v in cv.values()])
    print(f'querying WGS at {allsites.drop_duplicates().shape[0]} sites', flush=True)
    w = wgs_at(allsites, f'{OUT}/.sites.tsv')
    os.remove(f'{OUT}/.sites.tsv')

    conc_rows, map_rows, bc_rows, run_rows = [], [], [], []
    for (pool, mod) in runs.values:
        run = f'{pool}_{mod}'
        m = cv[run].merge(w, on=['chrom', 'pos', 'ref', 'alt'])
        disc = np.zeros((4, 4))
        for i in range(4):
            for j, d in enumerate(DONORS):
                ok = m[f'c{i}'].notna() & m[d].notna()
                n = int(ok.sum())
                c = float((m.loc[ok, f'c{i}'] == m.loc[ok, d]).mean()) if n else np.nan
                disc[i, j] = 1 - c
                conc_rows.append((pool, mod, i, d, n, c))
        r, c = linear_sum_assignment(disc)
        cmap = {int(i): DONORS[j] for i, j in zip(r, c)}
        for i in range(4):
            srt = np.sort(1 - disc[i])
            map_rows.append((pool, mod, i, cmap[i], 1 - disc[i, DONORS.index(cmap[i])], srt[-1] - srt[-2]))

        cl = pd.read_csv(f'{SCR}/{run}/clusters.tsv', sep='\t', usecols=['barcode', 'status', 'assignment'])
        cl['donor'] = np.where(cl.status == 'singlet',
                               cl.assignment.map(lambda a: cmap.get(int(a)) if str(a).isdigit() else None), None)
        cl.insert(0, 'mod', mod)
        cl.insert(0, 'pool', pool)
        bc_rows.append(cl)

        a = ambi[ambi['sample'] == pool].merge(cl, on='barcode', how='left', suffixes=('_ambi', ''))
        sing = a[a.status == 'singlet']
        run_rows.append(dict(pool=pool, mod=mod, n_soup_singlet=int((cl.status == 'singlet').sum()),
                             n_ambi_singlet=len(a), frac_ambi_soup_singlet=len(sing) / len(a),
                             concord_vs_ambi=float((sing.donor == sing.assignment_ambi).mean()),
                             min_cluster_wgs_conc=min(x[4] for x in map_rows[-4:]),
                             min_margin=min(x[5] for x in map_rows[-4:]),
                             map_one_to_one_matches_ambi_majority=bool(
                                 (sing.groupby('assignment')['assignment_ambi'].agg(lambda s: s.mode()[0])
                                  .rename(index=int).to_dict() == cmap) if len(sing) else False)))
        print(run, run_rows[-1], flush=True)

    conc = pd.DataFrame(conc_rows, columns=['pool', 'mod', 'cluster', 'donor', 'n_sites', 'gt_concordance'])
    cmapdf = pd.DataFrame(map_rows, columns=['pool', 'mod', 'cluster', 'donor', 'gt_concordance', 'margin_vs_2nd'])
    bc = pd.concat(bc_rows)
    rs = pd.DataFrame(run_rows)

    # GEX vs ATAC per barcode
    piv = bc.pivot_table(index=['pool', 'barcode'], columns='mod', values='donor', aggfunc='first')
    both = piv.dropna()
    xa = both.assign(agree=both.gex == both.atac).groupby('pool').agg(n_both_singlet=('agree', 'size'),
                                                                     gex_atac_agree=('agree', 'mean'))
    rs = rs.merge(xa, left_on='pool', right_index=True, how='left')
    # time label per (pool, donor) from the user's map
    tmap[['pool', 'donor']] = tmap['Sample'].str.split('_', expand=True)
    tmap['time'] = tmap['TimeSample2'].str.split('_').str[1]
    cmapdf = cmapdf.merge(tmap[['pool', 'donor', 'time']], on=['pool', 'donor'], how='left')

    conc.to_csv(f'{OUT}/cluster_donor_concordance.tsv', sep='\t', index=False)
    cmapdf.to_csv(f'{OUT}/cluster_map.tsv', sep='\t', index=False)
    bc.to_csv(f'{OUT}/barcode_calls.tsv.gz', sep='\t', index=False)
    rs.to_csv(f'{OUT}/run_summary.tsv', sep='\t', index=False)
    with open(f'{OUT}/summary.txt', 'w') as fh:
        fh.write('A09d souporcell (souporcell_gpu, k=4) cluster->donor mapping vs WGS and ambimux\n\n')
        fh.write(rs.to_string(index=False, float_format=lambda x: f'{x:.3f}') + '\n\n')
        fh.write(cmapdf.to_string(index=False, float_format=lambda x: f'{x:.3f}') + '\n')
    print(open(f'{OUT}/summary.txt').read())


if __name__ == '__main__':
    main()
