"""A06d: Level B targets. Emulate latent-genotype allele counts for the A06 cohort and write GLIMPSE2 GL targets.

Same cohort, regions and donor -> real-cluster assignment as A06a (profile_assign.tsv). For each test donor and
modality, at every site in its assigned cluster's souporcell site list (within the region's TSS +- half window):
  depth d   = that real cluster's AO + RO at the site (0 allowed: listed, no reads, as souporcell writes it)
  alt reads ~ Binomial(d, p_g),  p_g = (1 - rho) [e + (1 - 2e) g/2] + rho * a_s
              g = the donor's TRUE alt-allele count (1000G 30x), rho = the run's souporcell ambient fraction
              (ambient_rna.txt), a_s = cohort alt-allele frequency at the site (the ambient pool), e = 0.01
Genotype likelihoods use exactly the A04g ambient model (lib/souporcell_to_gl.py alt_probs / phred_pl), with a_s
re-estimated from the emulated counts pooled over donors, as the real pipeline does.
  multiome: per donor, GEX + ATAC counts summed per site (union of site lists), rho = mean of the two runs'.
Output (outdir/<region>/): target_<arm>.vcf.gz, one sample column per test donor; arms gex, atac, multiome.
Usage: python A06d_emulate_targets.py <a06a_outdir> [seed=2]   (allcools or scFates env; bcftools 1.11 on PATH)
"""
from __future__ import annotations

import gzip
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/u/project/cluo/terencew/claude/project_ideas/pool_design/scripts/latent_genos/lib')
from A06a_finemap_cohort import BCFTOOLS, POOLS, VCF  # noqa: E402
from souporcell_to_gl import HEADER, alt_probs, phred_pl  # noqa: E402

A03B = '/u/project/cluo/terencew/claude/project_ideas/latent_genos/results/A03b_souporcell_dosage_r2'
ERR = 0.01


def run_sites(run: str, chroms: set[str]) -> tuple[pd.DataFrame, dict[str, np.ndarray], float]:
    # listed SNV sites (chr16-22) with REF/ALT and per-cluster depth; plus the run's ambient fraction
    tree, rest = run.split('__', 1)
    pool, mod = rest.rsplit('__', 1)
    sdir = f'{POOLS}/{tree}/{pool}/demux/souporcell/{mod}'
    sdir = f'{sdir}/{os.listdir(sdir)[0]}'
    rho = float(re.search(r'([0-9.]+)%', open(f'{sdir}/ambient_rna.txt').read()).group(1)) / 100
    pat = '|'.join(sorted(chroms))
    txt = subprocess.run(f"awk -F'\\t' '$1 ~ /^({pat})$/' {sdir}/cluster_genotypes.vcf", shell=True,
                         capture_output=True, text=True, check=True).stdout
    recs = [l.split('\t') for l in txt.splitlines()]
    recs = [r for r in recs if 'BACKGROUND' not in r[6] and len(r[3]) == 1 and len(r[4]) == 1]
    sites = pd.DataFrame({'chrom': [r[0] for r in recs], 'pos': [int(r[1]) for r in recs],
                          'ref': [r[3] for r in recs], 'alt': [r[4] for r in recs]})
    k = len(recs[0]) - 9
    depth = {f'c{c}': np.array([sum(int(x) for x in r[9 + c].split(':')[1:3] if x not in ('', '.')) for r in recs],
                               dtype=np.int32) for c in range(k)}
    return sites, depth, rho


def main(a06a: str, seed: int = 2, half: int = 500_000) -> None:
    rng = np.random.default_rng(seed)
    reg = pd.read_csv(f'{a06a}/regions.tsv', sep='\t')
    asg = pd.read_csv(f'{a06a}/profile_assign.tsv', sep='\t')
    samples = asg['test_donor'].tolist()
    chroms = set(reg['chrom'])
    runs = sorted(set(asg['run_gex']) | set(asg['run_atac']))
    cache = {}
    for run in runs:
        cache[run] = run_sites(run, chroms)
        print(run, len(cache[run][0]), 'listed SNVs; rho', round(cache[run][2], 4), flush=True)

    for r in reg.itertuples():
        rd = f'{a06a}/{r.region}'
        if os.path.exists(f'{rd}/target_multiome.vcf.gz'):
            continue
        lo, hi = r.tss - half, r.tss + half
        # union of listed sites in the region over the runs used, with per-run index
        per_run = {}
        allsites = []
        for run in runs:
            s = cache[run][0]
            m = (s['chrom'] == r.chrom) & (s['pos'] >= lo) & (s['pos'] <= hi)
            per_run[run] = np.where(m.values)[0]
            allsites.append(s[m])
        u = pd.concat(allsites).drop_duplicates(['pos', 'ref', 'alt']).sort_values('pos').reset_index(drop=True)
        key = u['pos'].astype(str) + ':' + u['ref'] + ':' + u['alt']
        # true genotypes of the cohort at these sites (any MAF)
        tf = f'{rd}/emul_sites.tsv'
        u[['chrom', 'pos']].to_csv(tf, sep='\t', header=False, index=False)
        cmd = (f"{BCFTOOLS} view -s {','.join(samples)} -R {tf} -m2 -M2 -v snps {VCF.format(chrom=r.chrom)} | "
               f"{BCFTOOLS} query -f '%POS:%REF:%ALT[\\t%GT]\\n'")
        txt = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
        gt = {l.split('\t', 1)[0]: np.array([g.count('1') for g in l.rstrip('\n').split('\t')[1:]], dtype=np.int8)
              for l in txt.splitlines()}
        keep = key.isin(gt.keys()).values
        if not keep.any():   # e.g. pericentromeric regions: no listed site is a biallelic 1000G SNV
            for arm in ('gex', 'atac', 'multiome'):
                open(f'{rd}/target_{arm}.none', 'w').write('no emulated sites\n')
            open(f'{rd}/target_multiome.vcf.gz', 'w').close()   # completion marker for the skip check
            print(r.region, r.gene, 'no emulated sites; Level B arms = no information (2 x AF)', flush=True)
            continue
        u, key = u[keep].reset_index(drop=True), key[keep].reset_index(drop=True)
        G = np.vstack([gt[k] for k in key]).T            # donors x sites
        af = G.mean(0) / 2
        idx = {k: i for i, k in enumerate(key)}
        n, p = G.shape
        AO = {m: np.zeros((n, p), np.int32) for m in ('gex', 'atac')}
        RO = {m: np.zeros((n, p), np.int32) for m in ('gex', 'atac')}
        listed = {m: np.zeros((n, p), bool) for m in ('gex', 'atac')}
        rhos = {m: np.zeros(n) for m in ('gex', 'atac')}
        for j, a in enumerate(asg.itertuples()):
            for mod, run, cl in (('gex', a.run_gex, a.cluster_gex), ('atac', a.run_atac, a.cluster_atac)):
                s, dep, rho = cache[run]
                ii = per_run[run]
                ks = s['pos'].values[ii].astype(str)
                ks = pd.Series(ks) + ':' + s['ref'].values[ii] + ':' + s['alt'].values[ii]
                cols = np.array([idx.get(k, -1) for k in ks])
                ok = cols >= 0
                d = dep[cl][ii][ok]
                c = cols[ok]
                g = G[j, c]
                pg = (1 - rho) * (ERR + (1 - 2 * ERR) * g / 2) + rho * af[c]
                ao = rng.binomial(d, pg)
                AO[mod][j, c], RO[mod][j, c], listed[mod][j, c] = ao, d - ao, True
                rhos[mod][j] = rho
        arms = {'gex': (AO['gex'], RO['gex'], listed['gex'], rhos['gex'].mean()),
                'atac': (AO['atac'], RO['atac'], listed['atac'], rhos['atac'].mean()),
                'multiome': (AO['gex'] + AO['atac'], RO['gex'] + RO['atac'], listed['gex'] | listed['atac'],
                             (rhos['gex'].mean() + rhos['atac'].mean()) / 2)}
        for arm, (ao, ro, lst, rho) in arms.items():
            cols = np.where(lst.any(0))[0]
            if len(cols) == 0:
                open(f'{rd}/target_{arm}.none', 'w').write('no listed sites\n')
                continue
            out = f'{rd}/target_{arm}.vcf'
            with open(out, 'w') as fh:
                fh.write('\n'.join(HEADER) + '\n')
                fh.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples) + '\n')
                for c in cols:
                    tot = ao[:, c].sum() + ro[:, c].sum()
                    probs = alt_probs(ERR, rho, ao[:, c].sum() / tot if tot else 0.0)
                    calls = [f'./.:{a + b}:{b},{a}:{phred_pl(int(a), int(b), probs)}' for a, b in zip(ao[:, c], ro[:, c])]
                    fh.write('\t'.join([r.chrom, str(u['pos'][c]), '.', u['ref'][c], u['alt'][c], '.', 'PASS', '.',
                                        'GT:DP:AD:PL'] + calls) + '\n')
            subprocess.run(f'{BCFTOOLS} view -Oz -o {out}.gz {out} && rm {out}', shell=True, check=True)
        print(r.region, r.gene, 'emulated sites', p, 'gex listed/donor', int(listed['gex'].sum(1).mean()),
              'atac', int(listed['atac'].sum(1).mean()), 'mean depth gex (listed)',
              round(float((AO['gex'] + RO['gex'])[listed['gex']].mean()), 2), flush=True)


if __name__ == '__main__':
    main(sys.argv[1], *(int(x) for x in sys.argv[2:3]))
