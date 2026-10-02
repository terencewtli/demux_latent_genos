"""A05b: score demultiplexing with recovered genotypes (A05a) against the ambisim truth.

Per run (row of pool_design/txt/latent_genos_souporcell_tasks.txt) and method:
  true_GT     pool_design demuxlet on the TRUE pool genotypes (same pileup, same sites; the baseline)
  souporcell  souporcell's own clusters.tsv (no genotypes needed; what the latent genotypes came from)
  latent_GT / imputed_GT / imputed_GP   A05a demuxlet arms
Cluster-named calls (souporcell and the A05a arms) are mapped to donors with A03b assign.tsv; a cluster that A03b
flags ambiguous (e.g. a merged pair) keeps its matched donor, so its cells count as errors for the other donor.

Metrics (truth = drop_data_rand.txt, n = donors in the droplet; barcodes restricted to the filtered CR list):
  singlet_recall     true singlets called singlet with the right donor / true singlets
  singlet_precision  correct singlet calls / all singlet calls
  doublet_recall     true doublets called doublet (or ambiguous) / true doublets
  best_guess_acc     true singlets whose best single-donor guess is right, whatever the droplet call: genotype
                     informativeness apart from doublet calling. demuxlet (no ambient model) calls ~20% of true GEX
                     singlets DBL even with true genotypes while its best guess is right (2026-10-01), so singlet_recall
                     mostly measures doublet calling. souporcell has no best guess for its doublet calls (counted wrong).
  plus singlet_recall by reads-per-cell quartile (demuxlet NUM.READS from the true_GT run, so the bins are shared)
Output: results/A05b_score_demux/{summary,by_depth}.tsv
Usage (allcools env): python A05b_score_demux.py
"""
from __future__ import annotations

import os

import pandas as pd

PD = '/u/project/cluo/terencew/claude/project_ideas/pool_design'
LG = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SAMPLE = '20220928-IGVF-D0'
ARMS = ['latent_GT', 'imputed_GT', 'imputed_GP']
OUT = f'{LG}/results/A05b_score_demux'


def demuxlet_calls(path: str, cl2donor: dict[str, str] | None) -> pd.DataFrame:
    b = pd.read_csv(path, sep='\t', usecols=['BARCODE', 'NUM.READS', 'DROPLET.TYPE', 'SNG.BEST.GUESS'])
    call = b['SNG.BEST.GUESS'].astype(str)
    if cl2donor is not None:
        call = call.map(cl2donor)
    return pd.DataFrame({'barcode': b['BARCODE'].str.replace('-1', '', regex=False), 'reads': b['NUM.READS'],
                         'type': b['DROPLET.TYPE'].map({'SNG': 'singlet', 'DBL': 'doublet', 'AMB': 'ambiguous'}),
                         'donor': call})


def souporcell_calls(path: str, cl2donor: dict[str, str]) -> pd.DataFrame:
    s = pd.read_csv(path, sep='\t', usecols=['barcode', 'status', 'assignment'], dtype={'assignment': str})
    return pd.DataFrame({'barcode': s['barcode'].str.replace('-1', '', regex=False), 'reads': float('nan'),
                         'type': s['status'].replace({'unassigned': 'ambiguous'}),
                         'donor': ('c' + s['assignment']).map(cl2donor)})


def score(calls: pd.DataFrame, truth: pd.DataFrame) -> dict[str, float]:
    d = truth.merge(calls, on='barcode', how='left')
    sng, dbl = d['n'] == 1, d['n'] >= 2
    called_sng = d['type'] == 'singlet'
    correct = sng & called_sng & (d['donor'] == d['sam'])
    return {'n_true_singlets': int(sng.sum()), 'n_true_doublets': int(dbl.sum()),
            'singlet_recall': correct.sum() / sng.sum(), 'singlet_precision': correct.sum() / max(called_sng.sum(), 1),
            'best_guess_acc': (sng & (d['donor'] == d['sam'])).sum() / sng.sum(),
            'doublet_recall': (dbl & d['type'].isin(['doublet', 'ambiguous'])).sum() / max(dbl.sum(), 1),
            'missing_calls': int(d['type'].isna().sum())}


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    tasks = pd.read_csv(f'{PD}/txt/latent_genos_souporcell_tasks.txt', sep=' ', header=None, names=['tree', 'pool', 'k', 'mod'])
    rows, depth_rows = [], []
    for t in tasks.itertuples():
        tag = f'{t.tree}__{t.pool}__{t.mod}'
        a5 = f'{LG}/results/A05a_demuxlet_latent/{tag}'
        truth = pd.read_csv(f'{PD}/{t.tree}/{t.pool}/{SAMPLE}/drop_data_rand.txt', sep='\t', usecols=['RNA_BC', 'n', 'sam'],
                            dtype={'sam': str}).rename(columns={'RNA_BC': 'barcode'})
        bc = pd.read_csv(f'{PD}/{t.tree}/{t.pool}/cr_arc/{SAMPLE}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz', header=None)[0]
        truth = truth[truth['barcode'].isin(bc.str.replace('-1', '', regex=False))]
        assign = pd.read_csv(f'{LG}/results/A03b_souporcell_dosage_r2/{tag}/assign.tsv', sep='\t')
        cl2donor = dict(zip(assign['cluster'], assign['donor']))
        true_best = f'{PD}/{t.tree}/{t.pool}/demux/demuxlet/{t.mod}/{SAMPLE}.best'
        if not os.path.exists(true_best):
            print(f'{tag}: no true-genotype demuxlet baseline (no pileup); skipped')
            continue
        methods = {'true_GT': demuxlet_calls(true_best, None),
                   'souporcell': souporcell_calls(f'{PD}/{t.tree}/{t.pool}/demux/souporcell/{t.mod}/{SAMPLE}/clusters.tsv', cl2donor)}
        for arm in ARMS:
            if os.path.exists(f'{a5}/{arm}.best'):
                methods[arm] = demuxlet_calls(f'{a5}/{arm}.best', cl2donor)
        reads = methods['true_GT'].set_index('barcode')['reads']
        q = pd.qcut(reads, 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
        for m, calls in methods.items():
            rows.append({'run': tag, 'mod': t.mod, 'k': t.k, 'method': m, 'n_ambiguous_clusters': int(assign['ambiguous'].sum()),
                         **score(calls, truth)})
            for b in q.cat.categories:
                sub = truth[truth['barcode'].isin(q.index[q == b])]
                depth_rows.append({'run': tag, 'method': m, 'reads_quartile': b, 'median_reads': reads[q == b].median(),
                                   **score(calls, sub)})
        print(f'{tag}: {", ".join(methods)}')
    pd.DataFrame(rows).to_csv(f'{OUT}/summary.tsv', sep='\t', index=False)
    pd.DataFrame(depth_rows).to_csv(f'{OUT}/by_depth.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).pivot_table(index='run', columns='method', values='singlet_recall').round(3).to_string())


if __name__ == '__main__':
    main()
