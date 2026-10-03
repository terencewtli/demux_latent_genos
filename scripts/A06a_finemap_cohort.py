"""A06a: cohort, gene regions, true dosages and per-donor latent-genotype site classes for the A06 fine-mapping test.

Plan: docs/A06_FINEMAP_PLAN.md (mirror). Level A smoke test inputs; no reads are simulated.
1. Test cohort: n unrelated EUR from the 1000G 30x set (founders: no parent in the 3,202 set), seeded.
2. Regions: protein-coding genes on chr16-22 (GENCODE v43), TSS +- half_window, TSS >= min_gap apart.
3. Per region: biallelic SNPs, cohort MAF >= min_maf (default 5%: latent / imputed accuracy is only measured at
   MAF >= 5%, A04g), true dosage matrix (donors x SNPs).
4. Per donor, per SNP, a latent-genotype site class from a REAL souporcell cluster (A03a runs; cluster -> donor from
   A03b assign.tsv). Each test donor is assigned one real pool donor at random, and gets that donor's GEX and ATAC
   cluster depth (AO + RO) at each site:
     0 = site not in the cluster run's souporcell site list (untyped; filled by imputation)
     1 = in the list, depth 0; 2 = depth 1-2; 3 = 3-5; 4 = 6-10; 5 = > 10
   (the depth bins of results/A04g_glimpse_impute/summary_depth.ambient_b2m5.tsv).
Outputs (outdir/): samples.txt, regions.tsv, profile_assign.tsv, r<NNN>/{snps.tsv,truth.tsv.gz,class_gex.tsv.gz,class_atac.tsv.gz}
Usage: python A06a_finemap_cohort.py <outdir> [n=300] [n_regions=100] [seed=1]
"""
from __future__ import annotations

import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd

PD = '/u/project/cluo/terencew/demux_benchmark/pool_design'
VCF = PD + '/vcf/1000G/by_chrom/1000G.{chrom}.vcf.gz'
PED = PD + '/csv/1000G/meta/20130606_g1k_3202_samples_ped_population.txt'
GTF = '/u/project/cluo/terencew/reference/hg38_igvf/gencode.v43.chr_patch_hapl_scaff.annotation.gtf'
POOLS = '/u/project/cluo/terencew/claude/project_ideas/pool_design'
A03B = '/u/project/cluo/terencew/claude/project_ideas/latent_genos/results/A03b_souporcell_dosage_r2'
BCFTOOLS = '/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools'
CHROMS = [f'chr{i}' for i in range(16, 23)]
DEPTH_BINS = [0, 1, 3, 6, 11]  # class = 1 + number of bin edges <= depth, for listed sites


def cohort(n: int, rng: np.random.Generator) -> list[str]:
    ped = pd.read_csv(PED, sep=r'\s+')
    eur = ped[(ped['Superpopulation'] == 'EUR') & (ped['FatherID'].astype(str) == '0') & (ped['MotherID'].astype(str) == '0')]
    return sorted(rng.choice(eur['SampleID'].values, n, replace=False).tolist())


def regions(n: int, half: int, gap: int, rng: np.random.Generator) -> pd.DataFrame:
    cmd = (f"awk -F'\\t' '$3==\"gene\" && $1 ~ /^chr(1[6-9]|2[0-2])$/ && $9 ~ /gene_type \"protein_coding\"/' {GTF}")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split('\t')
        name = f[8].split('gene_name "')[1].split('"')[0]
        tss = int(f[3]) if f[6] == '+' else int(f[4])
        rows.append((f[0], tss, name))
    g = pd.DataFrame(rows, columns=['chrom', 'tss', 'gene']).drop_duplicates('gene')
    g = g.iloc[rng.permutation(len(g))]
    keep = []
    for r in g.itertuples():
        if r.tss - half < 1:
            continue
        if any(k.chrom == r.chrom and abs(k.tss - r.tss) < gap for k in keep):
            continue
        keep.append(r)
        if len(keep) == n:
            break
    df = pd.DataFrame(keep)[['chrom', 'tss', 'gene']]
    df['start'], df['end'] = df['tss'] - half, df['tss'] + half
    return df.sort_values(['chrom', 'start']).reset_index(drop=True)


def cluster_profiles() -> tuple[dict, pd.DataFrame]:
    # depth of every listed site (chr16-22) per (run, cluster); profiles = pool donors with both a GEX and an ATAC cluster
    depth, rows = {}, []
    for run in sorted(os.listdir(A03B)):
        if not run.endswith(('__gex', '__atac')):
            continue
        tree, pool, mod = run.split('__', 1)[0], run.split('__', 1)[1].rsplit('__', 1)[0], run.rsplit('__', 1)[1]
        a = pd.read_csv(f'{A03B}/{run}/assign.tsv', sep='\t')
        a = a[a['ambiguous'] == 0]
        sdir = f'{POOLS}/{tree}/{pool}/demux/souporcell/{mod}'
        vcf = f'{sdir}/{os.listdir(sdir)[0]}/cluster_genotypes.vcf'
        cmd = f"awk -F'\\t' '$1 ~ /^chr(1[6-9]|2[0-2])$/' {vcf}"
        txt = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
        recs = [l.split('\t') for l in txt.splitlines()]
        key = pd.Index([f'{r[0]}:{r[1]}' for r in recs])
        for c in a['cluster']:
            ci = 9 + int(c[1:])
            d = np.array([sum(int(x) if x not in ('', '.') else 0 for x in r[ci].split(':')[1:3]) for r in recs], dtype=np.int32)
            depth[(run, c)] = pd.Series(d, index=key)
        for r in a.itertuples():
            rows.append((tree + '__' + pool, r.donor, mod, run, r.cluster))
        print(run, len(recs), 'sites', len(a), 'clusters', flush=True)
    prof = pd.DataFrame(rows, columns=['pool', 'donor', 'mod', 'run', 'cluster'])
    prof = prof.pivot_table(index=['pool', 'donor'], columns='mod', values=['run', 'cluster'], aggfunc='first').dropna()
    prof.columns = [f'{a}_{b}' for a, b in prof.columns]
    return depth, prof.reset_index()


def site_class(dep: pd.Series, keys: pd.Index) -> np.ndarray:
    d = dep.reindex(keys)
    listed = d.notna().values
    cls = np.zeros(len(keys), dtype=np.int8)
    cls[listed] = 1 + np.searchsorted(DEPTH_BINS[1:], d.values[listed], side='right')
    return cls


def main(outdir: str, n: int = 300, n_regions: int = 100, seed: int = 1, half: int = 500_000, gap: int = 1_200_000,
         min_maf: float = 0.05) -> None:
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)
    samples = cohort(n, rng)
    open(f'{outdir}/samples.txt', 'w').write('\n'.join(samples) + '\n')
    reg = regions(n_regions, half, gap, rng)
    reg['region'] = [f'r{i:03d}' for i in range(len(reg))]
    reg.to_csv(f'{outdir}/regions.tsv', sep='\t', index=False)
    print('regions', len(reg), reg['chrom'].value_counts().to_dict(), flush=True)

    depth, prof = cluster_profiles()
    pick = rng.integers(0, len(prof), n)
    asg = prof.iloc[pick].reset_index(drop=True)
    asg.insert(0, 'test_donor', samples)
    asg.to_csv(f'{outdir}/profile_assign.tsv', sep='\t', index=False)
    print('profiles', len(prof), flush=True)

    for r in reg.itertuples():
        rd = f'{outdir}/{r.region}'
        if os.path.exists(f'{rd}/class_atac.tsv.gz'):
            continue
        os.makedirs(rd, exist_ok=True)
        cmd = (f"{BCFTOOLS} view -s {','.join(samples)} -r {r.chrom}:{r.start}-{r.end} -m2 -M2 -v snps "
               f"{VCF.format(chrom=r.chrom)} | {BCFTOOLS} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
        txt = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
        recs = [l.split('\t') for l in txt.splitlines()]
        pos = np.array([int(x[0]) for x in recs])
        G = np.array([[g.count('1') for g in x[3:]] for x in recs], dtype=np.int8).T  # donors x SNPs
        af = G.mean(0) / 2
        ok = np.minimum(af, 1 - af) >= min_maf
        snps = pd.DataFrame({'pos': pos[ok], 'ref': [x[1] for x, k in zip(recs, ok) if k],
                             'alt': [x[2] for x, k in zip(recs, ok) if k], 'af': af[ok]})
        snps.to_csv(f'{rd}/snps.tsv', sep='\t', index=False)
        G = G[:, ok]
        keys = pd.Index([f'{r.chrom}:{p}' for p in snps['pos']])
        cg = np.vstack([site_class(depth[(a.run_gex, a.cluster_gex)], keys) for a in asg.itertuples()])
        ca = np.vstack([site_class(depth[(a.run_atac, a.cluster_atac)], keys) for a in asg.itertuples()])
        for name, m in [('truth', G), ('class_gex', cg), ('class_atac', ca)]:
            with gzip.open(f'{rd}/{name}.tsv.gz', 'wt') as fh:
                np.savetxt(fh, m, fmt='%d', delimiter='\t')
        print(r.region, r.gene, r.chrom, len(snps), 'SNPs; GEX observed (depth>=1)', round(float((cg >= 2).mean()), 4),
              'ATAC', round(float((ca >= 2).mean()), 4), flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], *(int(x) for x in a[1:4]))
