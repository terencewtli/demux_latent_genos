"""OneK1K significant cis-eQTLs (A08a step 1 output) -> independent leads and fine-mapping loci.

Lead = per (cell type, gene, round) the SNP with the smallest P among localFDR < 0.05 rows, as in the paper's
conditional analysis (round r tests residuals after regressing out leads 1..r-1). Each lead is one of the
paper's "independent cis-eQTLs".

Outputs in results/A08_eqtl_anchor/:
  leads.tsv     every independent lead (cell type, gene, round), with chr/pos (GRCh37, see note) and z
  loci.tsv      fine-mapping loci for A08d: round-1 lead per (cell type, gene) with |z| >= --min-z, in --cts
  summary.tsv   counts per cell type (eGenes, leads, leads per round, typed fraction)

Note on build: the table's CHR/POS are from OneK1K's HRC imputation and are GRCh37; A08b lifts them.
z = sign(rho) * Phi^-1(1 - P/2) (two-sided Spearman P -> normal z).

usage: python scripts/A08a_select_eqtls.py [--min-z 5] [--cts cd4nc,cd8et,nk,bin,monoc]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
OUT = f'{PROJ}/results/A08_eqtl_anchor'
COLS = ['CELL_ID', 'RSID', 'SNPID', 'GENE', 'GENE_ID', 'CHR', 'POS', 'A1', 'A2', 'A2_FREQ_ONEK1K',
        'SPEARMANS_RHO', 'P_VALUE', 'Q_VALUE', 'FDR', 'RSQUARE', 'GENOTYPED', 'ROUND']


def load(path: str) -> pd.DataFrame:
    d = pd.read_csv(path, sep='\t', usecols=COLS)
    d['POS'] = d.POS.astype(int)
    d['z'] = np.sign(d.SPEARMANS_RHO) * norm.isf(d.P_VALUE.clip(lower=1e-300) / 2)
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-z', type=float, default=5.0)
    ap.add_argument('--cts', default='cd4nc,cd8et,nk,bin,monoc')
    a = ap.parse_args()
    files = sorted(glob.glob(f'{OUT}/sig/*.sig.tsv.gz'))
    print(f'{len(files)} cell types with significant tables', flush=True)
    leads, summ = [], []
    for f in files:
        d = load(f)
        L = d.sort_values('P_VALUE').groupby(['CELL_ID', 'GENE_ID', 'ROUND'], as_index=False).first()
        leads.append(L)
        row = {'cell_type': os.path.basename(f).split('.')[0], 'sig_rows': len(d), 'egenes': L.GENE_ID.nunique(),
               'leads': len(L), 'lead_typed_frac': (L.GENOTYPED == 'Genotyped').mean()}
        for r in range(1, 6):
            row[f'leads_round{r}'] = int((L.ROUND == r).sum())
        summ.append(row)
        print(row, flush=True)
    L = pd.concat(leads, ignore_index=True)
    L.to_csv(f'{OUT}/leads.tsv', sep='\t', index=False)
    pd.DataFrame(summ).to_csv(f'{OUT}/summary.tsv', sep='\t', index=False)
    cts = a.cts.split(',')
    loci = L[(L.ROUND == 1) & (L.z.abs() >= a.min_z) & L.CELL_ID.isin(cts)]
    loci.to_csv(f'{OUT}/loci.tsv', sep='\t', index=False)
    print(f'leads {len(L)}; loci for fine-mapping {len(loci)} ({loci.groupby("CELL_ID").size().to_dict()})', flush=True)


if __name__ == '__main__':
    main()
