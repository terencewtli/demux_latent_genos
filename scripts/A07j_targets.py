"""A07j: spacing vs per-site information. Why does latent-genotype imputation decay within ~5-20 kb of a read-covered
site? Builds GLIMPSE2 GL targets for a 2 x 2 of simulated inputs on the same 54 OneK1K pilot donors:

                       spacing = the donor's own souporcell sites      spacing = random sites, same count, MAF-matched
  perfect calls        perfect_soup                                    perfect_random
  weak calls           weak_soup                                       weak_random

- Donors: clean clusters of the 5 pilot pools (A07g preimpute assign.tsv, as A08b), renamed to OneK1K donors.
- Soup sites of a donor: sites of its pool's A07f target.bcf where its cluster has DP >= 1 (the sites that carry any
  information in the real run), restricted to candidates (below).
- Candidates: OneK1K biallelic SNVs with R2 >= MIN_R2 and MAF >= 1%, allele-matched to the 1000G panel sites.
  Input genotype = round(OneK1K DS). Random sites are drawn from a fixed random subset R of NON-typed candidates
  (so the typed sites stay available for scoring), MAF-bin matched to the donor's soup sites.
- perfect: PL 0 at the true genotype, PERFECT_PL elsewhere. weak: alt reads ~ Binomial(DP, p_g), p = (e, 0.5, 1 - e),
  with DP = the donor's real cluster depth at that site (soup) or its depths permuted onto the random sites
  (random); PL from the same binomial likelihoods. No ambient term, so weak_soup vs the real latent run is the cost
  of ambient RNA / cluster errors, and weak vs perfect is the cost of low depth.
- Sites absent for a donor get DP 0 and PL 0,0,0 (no information), as in the real targets.
Scoring is A07j_score.py (array-typed sites never used as input by any arm).
Output: <scratch>/A07j/targets/<arm>.<chrom>.bcf (+ .csi); results/A07j_spacing/inputs/<chrom>.input_sites.tsv.gz
usage: python scripts/A07j_targets.py --chrom chr22
"""
import argparse
import gzip
import os
import subprocess
from typing import Dict, List

import numpy as np
import pandas as pd

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos'
ONEK1K = f'{SCR}/A07/truth/OneK1K.noGP.vcf.gz'
PANEL = '/u/project/cluo/terencew/reference/topmed/local_1000G_30x'
POOLS = {'70': '_v2', '1': '', '55': '', '11': '_v2', '19': '_v2'}   # as A08b
MIN_R2 = 0.8
MIN_MAF = 0.01
ERR = 0.01
PERFECT_PL = 60
R_FACTOR = 3          # random subset R = R_FACTOR x the largest per-donor soup-site count
MAF_BINS = [0, 0.02, 0.05, 0.1, 0.2, 0.5]
ARMS = ['perfect_soup', 'perfect_random', 'weak_soup', 'weak_random']


def run(cmd: List[str]) -> str:
    return subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout


def sample_order(path: str, samples: List[str]) -> List[str]:
    return run(['bcftools', 'view', '-h', '-s', ','.join(samples), path]).splitlines()[-1].split('\t')[9:]


def pool_clusters(pool: str, sfx: str) -> Dict[str, str]:
    a = pd.read_csv(f'{PROJ}/results/A07_onek1k/preimpute/pool{pool}{sfx}/assign.tsv', sep='\t')
    a = a[(a.ambiguous == 0) & (a.r >= 0.6)]
    return dict(zip(a.cluster, a.donor))


def target_depth(pool: str, sfx: str, chrom: str, cmap: Dict[str, str]) -> pd.DataFrame:
    """Long table: donor, pos, ref, alt, dp (DP >= 1 only) from the A07f GL target of one pool."""
    bcf = f'{PROJ}/results/A07_onek1k/impute/pool{pool}{sfx}/{chrom}.target.bcf'
    keep = [c for c in run(['bcftools', 'query', '-l', bcf]).split() if c in cmap]
    order = sample_order(bcf, keep)
    txt = run(['bcftools', 'query', '-s', ','.join(keep), '-f', '%POS\t%REF\t%ALT[\t%DP]\n', bcf])
    d = pd.DataFrame([l.split('\t') for l in txt.splitlines()], columns=['pos', 'ref', 'alt'] + [cmap[c] for c in order])
    d = d.melt(id_vars=['pos', 'ref', 'alt'], var_name='donor', value_name='dp')
    d['pos'] = d.pos.astype(int)
    d['dp'] = pd.to_numeric(d.dp, errors='coerce').fillna(0).astype(int)
    return d[d.dp >= 1]


def candidates(chrom: str, donors: List[str]) -> pd.DataFrame:
    expr = f'TYPE="snp" && N_ALT=1 && INFO/R2>={MIN_R2} && INFO/MAF>={MIN_MAF}'
    order = sample_order(ONEK1K, donors)
    txt = run(['bcftools', 'query', '-r', chrom.replace('chr', ''), '-s', ','.join(donors), '-i', expr,
               '-f', '%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/TYPED[\t%DS]\n', ONEK1K])
    d = pd.DataFrame([l.split('\t') for l in txt.splitlines()], columns=['pos', 'ref', 'alt', 'maf', 'typed'] + order)
    d['pos'] = d.pos.astype(int)
    d['maf'] = d.maf.astype(float)
    d['typed'] = d.typed == '1'
    for s in order:
        d[s] = np.clip(np.rint(d[s].astype(float)), 0, 2).astype(np.int8)
    ps = run(['bcftools', 'query', '-f', '%POS\t%REF\t%ALT\n', f'{PANEL}/1000G_30x.{chrom}.sites.vcf.gz'])
    p = pd.DataFrame([l.split('\t') for l in ps.splitlines()], columns=['pos', 'ref', 'alt'])
    p['pos'] = p.pos.astype(int)
    d = d.merge(p.drop_duplicates(), on=['pos', 'ref', 'alt'])
    return d.drop_duplicates('pos', keep=False).reset_index(drop=True)


def pl(g: np.ndarray, dp: np.ndarray, weak: bool, rng: np.random.Generator) -> np.ndarray:
    """[n x 3] Phred-scaled likelihoods for true genotypes g at depths dp."""
    if not weak:
        out = np.full((len(g), 3), PERFECT_PL, dtype=int)
        out[np.arange(len(g)), g] = 0
        return out
    p = np.array([ERR, 0.5, 1 - ERR])
    a = rng.binomial(dp, p[g])
    ll = a[:, None] * np.log10(p)[None, :] + (dp - a)[:, None] * np.log10(1 - p)[None, :]
    return np.rint(-10 * (ll - ll.max(axis=1, keepdims=True))).astype(int)


def write_target(path: str, chrom: str, sites: pd.DataFrame, donors: List[str], P: Dict[str, Dict[int, tuple]]) -> None:
    """sites: pos, ref, alt (sorted); P[donor][pos] = (dp, pl0, pl1, pl2)."""
    tmp = path + '.tmp.vcf.gz'
    with gzip.open(tmp, 'wt') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write(f'##contig=<ID={chrom}>\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Unset; GLIMPSE2 reads PL">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Simulated depth (perfect arms: 1)">\n')
        f.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods (A07j simulation)">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(donors) + '\n')
        for pos, ref, alt in zip(sites.pos, sites.ref, sites.alt):
            cols = []
            for d in donors:
                v = P[d].get(pos)
                cols.append('./.:0:0,0,0' if v is None else f'./.:{v[0]}:{v[1]},{v[2]},{v[3]}')
            f.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP:PL\t' + '\t'.join(cols) + '\n')
    run(['bcftools', 'view', tmp, '-Ob', '-o', path])
    run(['bcftools', 'index', '-f', path])
    os.replace(tmp, path.replace('.bcf', '.vcf.gz'))   # keep the text version next to it for inspection


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    tdir = f'{SCR}/A07j/targets'
    idir = f'{PROJ}/results/A07j_spacing/inputs'
    os.makedirs(tdir, exist_ok=True)
    os.makedirs(idir, exist_ok=True)

    dep = []
    for pool, sfx in POOLS.items():
        cm = pool_clusters(pool, sfx)
        dep.append(target_depth(pool, sfx, a.chrom, cm).assign(pool=pool))
    dep = pd.concat(dep, ignore_index=True)
    donors = sorted(dep.donor.unique())
    C = candidates(a.chrom, donors)
    print(f'{a.chrom}: {len(donors)} donors; {len(C):,} candidates ({C.typed.sum():,} typed)', flush=True)
    C['maf_bin'] = pd.cut(C.maf, MAF_BINS, labels=False, include_lowest=True)
    cidx = pd.Series(np.arange(len(C)), index=C.pos)

    soup = dep.merge(C[['pos', 'ref', 'alt']], on=['pos', 'ref', 'alt'])
    n_d = soup.groupby('donor').size()
    pool_R = C[~C.typed & ~C.pos.isin(soup.pos)]
    R = pool_R.sample(n=min(len(pool_R), R_FACTOR * int(n_d.max())), random_state=a.seed)
    print(f'soup sites per donor: median {n_d.median():.0f} (range {n_d.min()}-{n_d.max()}); '
          f'union {soup.pos.nunique():,}; random subset R {len(R):,}', flush=True)

    rand_rows = []
    for d, s in soup.groupby('donor'):
        want = C.loc[cidx[s.pos].values, 'maf_bin'].value_counts()
        dps = rng.permutation(s.dp.values)
        picks = []
        for b, k in want.items():
            pool_b = R[R.maf_bin == b]
            picks.append(pool_b.sample(n=min(k, len(pool_b)), random_state=int(rng.integers(1 << 30))))
        pk = pd.concat(picks)
        rand_rows.append(pd.DataFrame({'donor': d, 'pos': pk.pos.values, 'dp': dps[:len(pk)]}))
    rand = pd.concat(rand_rows, ignore_index=True)

    for arm in ARMS:
        src = soup if arm.endswith('_soup') else rand
        weak = arm.startswith('weak')
        P: Dict[str, Dict[int, tuple]] = {d: {} for d in donors}
        for d, s in src.groupby('donor'):
            ii = cidx[s.pos].values
            g = C[d].values[ii].astype(int)
            L = pl(g, s.dp.values.astype(int), weak, rng)
            dpv = s.dp.values if weak else np.ones(len(s), dtype=int)
            P[d] = {int(p): (int(x), *map(int, l)) for p, x, l in zip(s.pos.values, dpv, L)}
        sites = C[C.pos.isin(src.pos)][['pos', 'ref', 'alt']].sort_values('pos')
        out = f'{tdir}/{arm}.{a.chrom}.bcf'
        write_target(out, a.chrom, sites, donors, P)
        print(f'{arm}: {len(sites):,} sites, {len(src):,} donor x site inputs -> {out}', flush=True)

    inp = pd.concat([soup[['donor', 'pos']].assign(spacing='soup'), rand[['donor', 'pos']].assign(spacing='random')])
    inp.groupby(['spacing', 'pos']).size().rename('n_donors').reset_index().to_csv(
        f'{idir}/{a.chrom}.input_sites.tsv.gz', sep='\t', index=False)
    dep.groupby(['pos']).size().rename('n_donors').reset_index().to_csv(
        f'{idir}/{a.chrom}.all_target_sites.tsv.gz', sep='\t', index=False)


if __name__ == '__main__':
    main()
