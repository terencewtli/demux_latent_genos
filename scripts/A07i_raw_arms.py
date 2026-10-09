"""A07i no-imputation baselines: per-cluster dosages at souporcell's own sites, written in the layout A07i_eqtl.py reads
for imputed arms (results/A07_onek1k/raw/<arm>/pool<P><sfx>/chrN.imputed.bcf, FORMAT/DS, samples c0..c{k-1}).

Arms:
  gt  souporcell cluster_genotypes.vcf hard calls (0/1/2; ./. = missing). Non-BACKGROUND biallelic SNVs. Nothing
      beyond souporcell itself: the "skip imputation" baseline.
  gl  posterior mean dosage at the GLIMPSE2 input sites (A07f chrN.target.bcf: the same ambient-model PLs GLIMPSE2
      sees) under an HWE prior from the panel's AF_EUR_unrel. Same per-site information + allele-frequency prior,
      but no haplotypes / LD: imputed minus gl isolates what LD-based imputation adds at the typed sites.

usage: python scripts/A07i_raw_arms.py --pool 1 [--suffix _v2] [--chroms chr22]
"""
import argparse
import os
import subprocess

import numpy as np
import pandas as pd

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A07'
RES = f'{PROJ}/results/A07_onek1k'
PANEL = '/u/project/cluo/terencew/reference/topmed/local_1000G_30x'
GT_DOSAGE = {'0/0': 0.0, '0/1': 1.0, '1/0': 1.0, '1/1': 2.0, '0|0': 0.0, '0|1': 1.0, '1|0': 1.0, '1|1': 2.0}


def log(msg: str) -> None:
    print(msg, flush=True)


def write_ds(rows: list, samples: list, out: str) -> None:
    # rows: (chrom, pos, id, ref, alt, [dosage or nan per sample]); bgzipped BCF + csi
    tmp = out + '.tmp.vcf'
    with open(tmp, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n')
        for c in [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']:
            fh.write(f'##contig=<ID={c}>\n')
        fh.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Alt allele dosage (no imputation)">\n')
        fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
        for chrom, pos, vid, ref, alt, ds in sorted(rows, key=lambda r: r[1]):
            vals = ['.' if np.isnan(x) else f'{x:.4f}' for x in ds]
            fh.write(f'{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\tDS\t' + '\t'.join(vals) + '\n')
    subprocess.run(['bcftools', 'view', '-Ob', '-o', out, tmp], check=True)
    subprocess.run(['bcftools', 'index', '-f', out], check=True)
    os.remove(tmp)


def gt_arm(soup_vcf: str, outdir: str, chroms: list) -> None:
    by_chrom = {}
    samples = []
    with open(soup_vcf) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'):
                samples = [s if s.startswith('c') else f'c{s}' for s in f[9:]]
                continue
            if 'BACKGROUND' in f[6] or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            i_gt = f[8].split(':').index('GT')
            ds = [GT_DOSAGE.get(s.split(':')[i_gt], np.nan) for s in f[9:]]
            by_chrom.setdefault(f[0], []).append((f[0], int(f[1]), f[2], f[3], f[4], ds))
    for chrom, rows in by_chrom.items():
        if chrom not in chroms:
            continue
        write_ds(rows, samples, f'{outdir}/{chrom}.imputed.bcf')
        log(f'gt {chrom}: {len(rows):,} sites, called {np.mean([not np.isnan(x) for r in rows for x in r[5]]):.3f}')


def panel_af(chrom: str) -> dict:
    q = subprocess.run(['bcftools', 'query', '-f', '%POS\t%REF\t%ALT\t%INFO/AF_EUR_unrel\n',
                        f'{PANEL}/1000G_30x.{chrom}.sites.vcf.gz'], capture_output=True, text=True, check=True).stdout
    af = {}
    for l in q.splitlines():
        p, r, a, x = l.split('\t')
        af[(int(p), r, a)] = float(x)
    return af


def gl_arm(imp: str, outdir: str, chroms: list) -> None:
    for chrom in chroms:
        target = f'{imp}/{chrom}.target.bcf'
        if not os.path.exists(target):
            log(f'gl {chrom}: no target.bcf, skipped')
            continue
        samples = subprocess.run(['bcftools', 'query', '-l', target], capture_output=True, text=True, check=True).stdout.split()
        q = subprocess.run(['bcftools', 'query', '-f', '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%PL]\n', target],
                           capture_output=True, text=True, check=True).stdout
        af = panel_af(chrom)
        rows, n_noaf = [], 0
        for l in q.splitlines():
            f = l.split('\t')
            p = af.get((int(f[1]), f[3], f[4]))
            if p is None:
                n_noaf += 1
                continue
            p = min(max(p, 1e-4), 1 - 1e-4)
            prior = np.log10([(1 - p) ** 2, 2 * p * (1 - p), p ** 2])
            pl = np.array([[float(x) for x in s.split(',')] for s in f[5:]])  # samples x 3
            lp = -pl / 10 + prior
            w = 10 ** (lp - lp.max(axis=1, keepdims=True))
            ds = (w @ np.array([0.0, 1.0, 2.0])) / w.sum(axis=1)
            rows.append((f[0], int(f[1]), f[2], f[3], f[4], list(ds)))
        write_ds(rows, samples, f'{outdir}/{chrom}.imputed.bcf')
        log(f'gl {chrom}: {len(rows):,} sites ({n_noaf} without panel AF dropped)')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--pool', required=True)
    ap.add_argument('--suffix', default='')
    ap.add_argument('--arms', default='gt,gl')
    ap.add_argument('--chroms', default=','.join(f'chr{i}' for i in range(1, 23)))
    a = ap.parse_args()
    tag = f'pool{a.pool}{a.suffix}'
    for arm in a.arms.split(','):
        outdir = f'{RES}/raw/{arm}/{tag}'
        os.makedirs(outdir, exist_ok=True)
        if arm == 'gt':
            gt_arm(f'{SCR}/souporcell/{tag}/cluster_genotypes.vcf', outdir, a.chroms.split(','))
        else:
            gl_arm(f'{RES}/impute/{tag}', outdir, a.chroms.split(','))


if __name__ == '__main__':
    main()
