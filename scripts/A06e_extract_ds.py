"""A06e helper: GLIMPSE2 dosages at the A06 cohort SNPs -> donors x SNPs matrix (same layout as truth.tsv.gz).

Rows follow <a06a_outdir>/samples.txt (matched by sample name, not position); columns follow <region>/snps.tsv
(matched by POS:REF:ALT). SNPs missing from the imputed output (absent from the leave-out panel) get 2 x cohort AF,
i.e. no information; the count is printed.
Usage: python A06e_extract_ds.py <imputed.bcf | NONE> <a06a_outdir> <region> <out.tsv.gz>
(NONE: the arm had no target sites, so every SNP gets 2 x AF.)
"""
from __future__ import annotations

import gzip
import subprocess
import sys

import numpy as np
import pandas as pd

BCFTOOLS = '/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools'


def main(bcf: str, a06a: str, region: str, out: str) -> None:
    samples = open(f'{a06a}/samples.txt').read().split()
    snps = pd.read_csv(f'{a06a}/{region}/snps.tsv', sep='\t')
    keys = snps['pos'].astype(str) + ':' + snps['ref'] + ':' + snps['alt']
    ds = {}
    txt = ''
    if bcf != 'NONE':
        hdr = subprocess.run([BCFTOOLS, 'query', '-l', bcf], capture_output=True, text=True, check=True).stdout.split()
        order = [hdr.index(s) for s in samples]
        txt = subprocess.run([BCFTOOLS, 'query', '-f', '%POS:%REF:%ALT[\t%DS]\n', bcf], capture_output=True,
                             text=True, check=True).stdout
    for line in txt.splitlines():
        f = line.split('\t')
        ds[f[0]] = np.array(f[1:], dtype=float)[order]
    n, p = len(samples), len(snps)
    X = np.empty((n, p))
    miss = 0
    for j, k in enumerate(keys):
        if k in ds:
            X[:, j] = ds[k]
        else:
            X[:, j] = 2 * snps['af'].iloc[j]
            miss += 1
    with gzip.open(out, 'wt') as fh:
        np.savetxt(fh, X, fmt='%.3f', delimiter='\t')
    print(region, out.rsplit('/', 1)[-1], 'SNPs', p, 'missing from imputed output', miss, flush=True)


if __name__ == '__main__':
    main(*sys.argv[1:5])
