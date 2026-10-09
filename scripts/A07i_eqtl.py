"""A07i step 2: cis-eQTL mapping (tensorQTL) with each genotype arm, compared against OneK1K's own array-imputed
genotypes on the same donors and the same phenotypes - the scTAPAS comparison (Kockelbergh 2026, Fig. 3), here on
3' data and, for the soup arms, without known genotypes.

Genotype arms (chroms of --chroms only):
  array          OneK1K.noGP.vcf.gz (array + Minimac4), DS, imputation R2 >= 0.8: the reference
  <A07h arms>    DS from results/A07_onek1k/{impute,impute_bam/*,impute_quilt/*}/pool<P>/chrN.imputed.bcf
                 (cluster names -> donors via A07g assign.tsv; see A07h_score.ARMS)
  raw_gt, raw_gl no-imputation baselines at souporcell's own sites (A07i_raw_arms.py): hard calls, and the
                 single-site posterior mean (same PLs as GLIMPSE2's input + EUR AF prior, no LD)
Pools: pool<P>_v2 (A07d rerun) is used when its A07g score/assign.tsv exists, else pool<P>.
Missing dosages (raw arms; pools with different site sets): variants need call rate >= MIN_CALL among the analysed
donors, and the rest are filled with the variant mean (standard for missing genotypes; a no-op for complete arms).
Phenotypes: A07i_pseudobulk counts (--labels oracle by default, so every arm shares one phenotype matrix and only
genotypes differ; --labels soup gives the fully genotype-free pipeline). Per cell type: genes with CPM >= 1 in
>= 50% of donors; log2(CPM + 1); rank-inverse-normal per gene. Covariates: N_PC expression PCs + pool indicators.
Variants: biallelic SNVs with MAF >= 5% among the analysed donors and present in the arm. TSS from GENCODE v43.
tensorQTL cis.map_cis (nperm 1000, 1 Mb window); eGene = BH FDR < 0.05 on pval_beta within arm x cell type.
Comparison with the array arm: recovery of array eGenes, extra eGenes, slope correlation of shared eGenes, and LD
(array-dosage r2 among the analysed donors) between the two lead variants.
Output: results/A07_onek1k/eqtl/<labels>/{<arm>.<celltype>.cis.tsv.gz, compare.tsv}
usage: (tensorqtl env) python scripts/A07i_eqtl.py --chroms chr6,chr22 [--labels oracle] [--arms all]
"""
import argparse
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A07'
RES = f'{PROJ}/results/A07_onek1k'
ARRAY = f'{SCR}/truth/OneK1K.noGP.vcf.gz'
GTF = '/u/project/cluo/terencew/reference/hg38_igvf/gencode.v43.chr_patch_hapl_scaff.annotation.gtf'
ARMS = {'soupsites_glimpse': f'{RES}/impute', 'bam_glimpse_soup': f'{RES}/impute_bam/soup',
        'bam_glimpse_oracle': f'{RES}/impute_bam/oracle', 'bam_quilt_soup': f'{RES}/impute_quilt/soup',
        'bam_quilt_oracle': f'{RES}/impute_quilt/oracle', 'raw_gt': f'{RES}/raw/gt', 'raw_gl': f'{RES}/raw/gl'}
CELLTYPES = ['all', 'CD4T', 'CD8T', 'NK', 'B', 'Mono']
N_PC = 5
MIN_MAF = 0.05
MIN_CALL = 0.5


def pool_tag(p: str) -> str:
    return f'pool{p}_v2' if os.path.exists(f'{RES}/score/pool{p}_v2/assign.tsv') else f'pool{p}'


def bcf_ds(vcf: str, region: str, samples: list = None, extra: list = ()) -> pd.DataFrame:
    names = samples or subprocess.run(['bcftools', 'query', '-l', vcf], capture_output=True, text=True, check=True).stdout.split()
    cmd = ['bcftools', 'view', '-r', region, '-m2', '-M2', '-v', 'snps', *extra]
    if samples:
        cmd += ['-s', ','.join(samples)]
    v = subprocess.Popen(cmd + ['-Ou', vcf], stdout=subprocess.PIPE)
    q = subprocess.run(['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n'], stdin=v.stdout,
                       capture_output=True, text=True, check=True).stdout
    rows = [l.split('\t') for l in q.splitlines()]
    d = pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt'] + names)
    c = d.chrom.str.replace('^chr', '', regex=True)
    d.index = 'chr' + c + ':' + d.pos + ':' + d.ref + ':' + d.alt
    d = d[~d.index.duplicated()]
    return d[names].apply(pd.to_numeric, errors='coerce')


def arm_dosage(arm: str, chroms: list, donors: list) -> pd.DataFrame:
    if arm == 'array':
        parts = [bcf_ds(ARRAY, c.replace('chr', ''), donors, ['-i', 'INFO/R2>=0.8 || INFO/TYPED_ONLY=1']) for c in chroms]
        return pd.concat(parts)
    pools = pd.read_csv(f'{PROJ}/txt/onek1k_pilot_pools.txt', sep='\t', header=None, names=['pool', 'gsm', 'k'], dtype=str)
    out = []
    for c in chroms:
        mats = []
        for p in pools.pool:
            tag = pool_tag(p)
            f = f'{ARMS[arm]}/{tag}/{c}.imputed.bcf'
            if not os.path.exists(f + '.csi') or not os.path.exists(f'{RES}/score/{tag}/assign.tsv'):
                continue
            A = pd.read_csv(f'{RES}/score/{tag}/assign.tsv', sep='\t')
            c2d = dict(zip(A.cluster.astype(str), A.donor))
            ds = bcf_ds(f, c)
            ds = ds.rename(columns={s: (s if s.startswith('OneK1K_') else c2d.get(s if s.startswith('c') else f'c{s}'))
                                    for s in ds.columns})
            mats.append(ds[[s for s in ds.columns if s in donors]])
        if mats:
            out.append(pd.concat(mats, axis=1))
    if not out:
        return pd.DataFrame()
    D = pd.concat(out)
    return D.loc[:, ~D.columns.duplicated()]


def tss(chroms: list) -> pd.DataFrame:
    rows = []
    with open(GTF) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 9)
            if f[2] != 'gene' or f[0] not in chroms:
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            t = int(f[3]) if f[6] == '+' else int(f[4])
            rows.append((gid, f[0], t, t))
    return pd.DataFrame(rows, columns=['gene_id', 'chr', 'start', 'end']).set_index('gene_id')


def phenotypes(labels: str, ct: str, pos: pd.DataFrame) -> pd.DataFrame:
    m = pd.read_csv(f'{RES}/eqtl/pseudobulk/{labels}/{ct}.counts.tsv.gz', sep='\t', index_col=0)
    cpm = m / m.sum() * 1e6
    keep = (cpm >= 1).mean(axis=1) >= 0.5
    x = np.log2(cpm[keep] + 1)
    x = x.loc[x.index.intersection(pos.index)]
    r = x.apply(lambda v: pd.Series(norm.ppf((rankdata(v) - 0.5) / len(v)), index=v.index), axis=1)
    return r


def covariates(ph: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    z = ph.sub(ph.mean(axis=1), axis=0)
    u, s, vt = np.linalg.svd(z.values, full_matrices=False)
    pcs = pd.DataFrame(vt[:N_PC].T, index=ph.columns, columns=[f'PC{i + 1}' for i in range(N_PC)])
    pool = cells.dropna(subset=['donor']).groupby('donor').pool.agg(lambda s: s.value_counts().index[0])
    dummies = pd.get_dummies(pool.reindex(ph.columns), prefix='pool', drop_first=True).astype(float)
    return pd.concat([pcs, dummies], axis=1)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def main() -> None:
    from tensorqtl import cis
    ap = argparse.ArgumentParser()
    ap.add_argument('--chroms', default='chr6,chr22')
    ap.add_argument('--labels', default='oracle', choices=['oracle', 'soup'])
    ap.add_argument('--arms', default='all')
    a = ap.parse_args()
    chroms = a.chroms.split(',')
    out = f'{RES}/eqtl/{a.labels}'
    os.makedirs(out, exist_ok=True)
    arms = ['array'] + (list(ARMS) if a.arms == 'all' else a.arms.split(','))
    pos = tss(chroms)
    cells = pd.read_csv(f'{RES}/eqtl/pseudobulk/{a.labels}/cells.tsv.gz', sep='\t', dtype={'pool': str})
    donors = sorted(cells.donor.dropna().unique())
    G = {}
    for arm in arms:
        d = arm_dosage(arm, chroms, donors)
        if d.empty:
            print(f'{arm}: no dosages yet, skipped', flush=True)
            continue
        G[arm] = d
        print(f'{arm}: {d.shape[0]:,} variants x {d.shape[1]} donors', flush=True)
    res = {}
    for ct in CELLTYPES:
        ph_all = phenotypes(a.labels, ct, pos)
        for arm, g in G.items():
            ds = [s for s in ph_all.columns if s in g.columns]
            if len(ds) < 20:
                continue
            ph = ph_all[ds]
            gg = g[ds]
            gg = gg[gg.notna().mean(axis=1) >= MIN_CALL]
            gg = gg.apply(lambda v: v.fillna(v.mean()), axis=1)
            af = gg.mean(axis=1) / 2
            gg = gg[(np.minimum(af, 1 - af) >= MIN_MAF)]
            vdf = pd.DataFrame({'chrom': gg.index.str.split(':').str[0], 'pos': gg.index.str.split(':').str[1].astype(int)},
                               index=gg.index)
            cov = covariates(ph, cells)
            r = cis.map_cis(gg, vdf, ph, pos.loc[ph.index], covariates_df=cov, window=1_000_000, nperm=1000, seed=1,
                            verbose=False)
            r['qval_bh'] = bh(r.pval_beta.values)
            r.to_csv(f'{out}/{arm}.{ct}.cis.tsv.gz', sep='\t')
            res[(arm, ct)] = r
            print(f'{ct} {arm}: {len(ds)} donors, {len(ph)} genes, {len(gg):,} variants, eGenes '
                  f'{int((r.qval_bh < 0.05).sum())}', flush=True)
    rows = []
    for (arm, ct), r in res.items():
        if arm == 'array' or ('array', ct) not in res:
            continue
        ref = res[('array', ct)]
        e_ref = set(ref.index[ref.qval_bh < 0.05])
        e_arm = set(r.index[r.qval_bh < 0.05])
        shared = sorted(e_ref & e_arm)
        ld = []
        for gid in shared:
            v1, v2 = ref.loc[gid, 'variant_id'], r.loc[gid, 'variant_id']
            if v1 in G['array'].index and v2 in G['array'].index:
                x, y = G['array'].loc[v1], G['array'].loc[v2]
                ld.append(np.corrcoef(x, y)[0, 1] ** 2)
        rows.append({'labels': a.labels, 'celltype': ct, 'arm': arm, 'array_eGenes': len(e_ref), 'arm_eGenes': len(e_arm),
                     'shared': len(shared), 'recovery': len(shared) / len(e_ref) if e_ref else np.nan,
                     'arm_only': len(e_arm - e_ref),
                     'slope_r2_shared': (np.corrcoef(ref.loc[shared, 'slope'], r.loc[shared, 'slope'])[0, 1] ** 2
                                         if len(shared) > 2 else np.nan),
                     'same_lead': float(np.mean([ref.loc[s, 'variant_id'] == r.loc[s, 'variant_id'] for s in shared]))
                     if shared else np.nan,
                     'lead_ld_median': float(np.median(ld)) if ld else np.nan})
    C = pd.DataFrame(rows)
    C.to_csv(f'{out}/compare.tsv', sep='\t', index=False, float_format='%.3f')
    pd.set_option('display.width', 200)
    print(C.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
