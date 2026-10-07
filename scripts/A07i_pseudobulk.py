"""A07i step 1: donor x cell-type pseudobulk counts from the pilot STARsolo matrices, for cis-eQTL mapping (A07i_eqtl.py).

Cells -> donors, two ways (same as A07h):
  oracle  GEO per-cell donor labels (convention resolved from A07g assign.tsv)
  soup    souporcell singlets, cluster -> donor via A07g assign.tsv (genotype-free labels)
Cell types: a marker-score lineage call per cell (pilot-scale stand-in for OneK1K's annotation): log1p CP10k, mean
over each marker set, argmax if > 0.5 else 'other'. T cells split CD8 / CD4 by CD8A + CD8B.
  T: CD3D CD3E    NK: NKG7 GNLY KLRD1 (and not T)    B: MS4A1 CD79A CD79B    Mono: CD14 LYZ FCGR3A S100A8
Pseudobulk = summed UMI counts per donor x cell type (types: all, CD4T, CD8T, NK, B, Mono), donors with >= MIN_CELLS.
Output: results/A07_onek1k/eqtl/pseudobulk/<labels>/<celltype>.counts.tsv.gz (genes x donors),
        cells.tsv.gz (pool, barcode, donor, lineage), n_cells.tsv
usage: python scripts/A07i_pseudobulk.py --labels oracle
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A07h_split_bams import oracle_labels, soup_labels  # noqa: E402

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A07'
RES = f'{PROJ}/results/A07_onek1k'
MIN_CELLS = 10
MARKERS = {'T': ['CD3D', 'CD3E'], 'NK': ['NKG7', 'GNLY', 'KLRD1'], 'B': ['MS4A1', 'CD79A', 'CD79B'],
           'Mono': ['CD14', 'LYZ', 'FCGR3A', 'S100A8']}


def lineage(m: sp.csr_matrix, symbols: np.ndarray) -> np.ndarray:
    lib = np.asarray(m.sum(axis=0)).ravel()
    x = m.multiply(1e4 / np.maximum(lib, 1)).tocsr()
    x.data = np.log1p(x.data)
    score = {}
    for k, genes in MARKERS.items():
        idx = [i for i, s in enumerate(symbols) if s in genes]
        score[k] = np.asarray(x[idx].mean(axis=0)).ravel()
    S = pd.DataFrame(score)
    S.loc[S['T'] > 0.5, 'NK'] = 0
    lab = np.where(S.max(axis=1) > 0.5, S.idxmax(axis=1), 'other')
    cd8 = [i for i, s in enumerate(symbols) if s in ('CD8A', 'CD8B')]
    cd8s = np.asarray(x[cd8].mean(axis=0)).ravel()
    return np.where(lab == 'T', np.where(cd8s > 0.5, 'CD8T', 'CD4T'), lab)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--labels', choices=['oracle', 'soup'], required=True)
    a = p.parse_args()
    out = f'{RES}/eqtl/pseudobulk/{a.labels}'
    os.makedirs(out, exist_ok=True)
    pools = pd.read_csv(f'{PROJ}/txt/onek1k_pilot_pools.txt', sep='\t', header=None, names=['pool', 'gsm', 'k'], dtype=str)
    sums, cells, genes = {}, [], None
    for pool in pools.pool:
        d = f'{SCR}/starsolo/pool{pool}/Solo.out/Gene/filtered'
        m = scipy.io.mmread(f'{d}/matrix.mtx').tocsr()
        f = pd.read_csv(f'{d}/features.tsv', sep='\t', header=None)
        bc = pd.read_csv(f'{d}/barcodes.tsv', header=None)[0].str.replace(r'-1$', '', regex=True).values
        genes = f[0].values
        if a.labels == 'oracle':
            lab = oracle_labels(pool)
        else:
            A = pd.read_csv(f'{RES}/score/pool{pool}/assign.tsv', sep='\t')
            c2d = dict(zip(A.cluster.astype(str), A.donor))
            lab = {b: c2d.get(c) for b, c in soup_labels(pool).items()}
        donor = np.array([lab.get(b) for b in bc], dtype=object)
        lin = lineage(m, f[1].values)
        cells.append(pd.DataFrame({'pool': pool, 'barcode': bc, 'donor': donor, 'lineage': lin}))
        for ct in ['all', 'CD4T', 'CD8T', 'NK', 'B', 'Mono']:
            for dn in pd.unique(donor[donor != None]):  # noqa: E711
                sel = (donor == dn) & ((lin == ct) if ct != 'all' else (lin != 'other'))
                if sel.sum() < MIN_CELLS:
                    continue
                sums[(ct, dn)] = (np.asarray(m[:, np.where(sel)[0]].sum(axis=1)).ravel(), int(sel.sum()))
        print(f'pool {pool}: {len(bc)} cells, labelled {np.sum(donor != None)}, lineages '  # noqa: E711
              f'{pd.Series(lin).value_counts().to_dict()}', flush=True)
    C = pd.concat(cells)
    C.to_csv(f'{out}/cells.tsv.gz', sep='\t', index=False)
    nrow = []
    for ct in ['all', 'CD4T', 'CD8T', 'NK', 'B', 'Mono']:
        ks = [k for k in sums if k[0] == ct]
        M = pd.DataFrame({k[1]: sums[k][0] for k in ks}, index=genes)
        M.index.name = 'gene_id'
        M.to_csv(f'{out}/{ct}.counts.tsv.gz', sep='\t')
        nrow += [{'celltype': ct, 'donor': k[1], 'cells': sums[k][1]} for k in ks]
    N = pd.DataFrame(nrow)
    N.to_csv(f'{out}/n_cells.tsv', sep='\t', index=False)
    print(N.groupby('celltype').agg(donors=('donor', 'size'), median_cells=('cells', 'median')).to_string())


if __name__ == '__main__':
    main()
