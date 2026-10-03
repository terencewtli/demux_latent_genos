"""A06c: score the A06 fine-mapping / coloc test (docs/A06_FINEMAP_PLAN.md) against truth and apply the go/no-go rule.

Input: <outdir>/r*.tsv.gz from A06b (one row per region x scenario x h2 x seed x arm [x coloc type]).
Every arm is compared with `truth` on the SAME phenotype (paired by region, scenario, h2, seed), so differences are
due to the genotypes only.
Outputs (results dir):
  summary.tsv     mean (and SE over phenotypes) of each metric by arm x scenario x h2
  paired.tsv      arm - truth differences for the key metrics, by scenario x h2
  coloc.tsv       PP.H4 > 0.8 rate by arm x coloc type (shared = sensitivity, distinct = false colocalisation)
  gonogo.tsv      the plan's criteria evaluated on the realistic cells (untyped causal; h2 0.05 / 0.10)
Usage: python A06c_score_finemap.py <a06b_outdir> <results_dir>
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

KEYS = ['region', 'scenario', 'h2', 'seed']
METRICS = ['egene_bonf', 'n_cs', 'causal_in_cs', 'all_causal_in_cs', 'cs_size_median', 'causal_pip_mean',
           'causal_pip_gt05', 'lead_is_causal', 'lead_ld_causal', 'lead_observed_gex']


def se(x: pd.Series) -> float:
    x = x.dropna()
    return float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else np.nan


def main(indir: str, resdir: str) -> None:
    os.makedirs(resdir, exist_ok=True)
    d = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{indir}/r*.tsv.gz'))], ignore_index=True)
    for c in ['egene_bonf', 'all_causal_in_cs', 'lead_is_causal', 'lead_observed_gex', 'susie_ok']:
        d[c] = d[c].astype(str).str.upper().eq('TRUE').astype(float)
    print('rows', len(d), 'regions', d['region'].nunique(), 'arms', sorted(d['arm'].unique()), flush=True)

    fm = d.drop_duplicates(KEYS + ['arm'])          # coloc rows duplicate the fine-mapping row
    g = fm.groupby(['arm', 'scenario', 'h2'])
    summ = g[METRICS].mean().join(g[METRICS].agg(se).add_suffix('_se')).join(g.size().rename('n_pheno'))
    summ = summ.join(fm.groupby(['arm', 'scenario', 'h2'])['realized_r2'].first())
    summ.to_csv(f'{resdir}/summary.tsv', sep='\t')

    tr = fm[fm['arm'] == 'truth'].set_index(KEYS)[METRICS]
    rows = []
    for arm, a in fm[fm['arm'] != 'truth'].groupby('arm'):
        diff = a.set_index(KEYS)[METRICS] - tr.reindex(a.set_index(KEYS).index)
        diff = diff.reset_index()
        for (sc, h2), x in diff.groupby(['scenario', 'h2']):
            row = {'arm': arm, 'scenario': sc, 'h2': h2, 'n_pheno': len(x)}
            for m in METRICS:
                row[f'd_{m}'] = x[m].mean()
                row[f'd_{m}_se'] = se(x[m])
            rows.append(row)
    paired = pd.DataFrame(rows)
    paired.to_csv(f'{resdir}/paired.tsv', sep='\t', index=False)

    c = d.dropna(subset=['coloc_type'])
    coloc = c.assign(h4=(c['pp_h4'] > 0.8).astype(float)).groupby(['arm', 'coloc_type', 'scenario'])['h4'].agg(['mean', se, 'size'])
    coloc.columns = ['pp_h4_gt08', 'se', 'n']
    coloc.to_csv(f'{resdir}/coloc.tsv', sep='\t')

    # go / no-go on realistic cells: untyped causal, h2 0.05 and 0.10, per imputed arm
    real = paired[(paired['scenario'] == 'untyped') & (paired['h2'].isin([0.05, 0.10]))]
    tr_obs = summ.loc['truth'].reset_index()
    tr_obs = tr_obs[(tr_obs['scenario'] == 'untyped') & (tr_obs['h2'].isin([0.05, 0.10]))]['lead_observed_gex'].mean()
    go = []
    for arm, x in real.groupby('arm'):
        arm_obs = summ.loc[arm].reset_index()
        arm_obs = arm_obs[(arm_obs['scenario'] == 'untyped') & (arm_obs['h2'].isin([0.05, 0.10]))]['lead_observed_gex'].mean()
        cs_drop = -x['d_causal_in_cs'].mean()
        cl = coloc.reset_index()
        cl_arm = cl[(cl['arm'] == arm) & (cl['scenario'] == 'untyped')].set_index('coloc_type')['pp_h4_gt08']
        cl_tr = cl[(cl['arm'] == 'truth') & (cl['scenario'] == 'untyped')].set_index('coloc_type')['pp_h4_gt08']
        coloc_shift = float(np.nanmax(np.abs((cl_arm - cl_tr).values))) if len(cl_arm) else np.nan
        go.append({'arm': arm, 'cs_coverage_drop': cs_drop,
                   'lead_observed_rate': arm_obs, 'lead_observed_rate_truth': tr_obs,
                   'lead_shift_ratio': arm_obs / tr_obs if tr_obs > 0 else np.nan,
                   'coloc_max_abs_change': coloc_shift,
                   # lead shift: >= 2x the truth rate AND >= 10 points higher (handles a truth rate of 0)
                   'GO': (cs_drop >= 0.15) or (arm_obs >= 2 * tr_obs and arm_obs - tr_obs >= 0.10) or (coloc_shift >= 0.10)})
    go = pd.DataFrame(go)
    go.to_csv(f'{resdir}/gonogo.tsv', sep='\t', index=False)
    pd.set_option('display.width', 200)
    print(summ[['causal_in_cs', 'cs_size_median', 'causal_pip_mean', 'lead_is_causal', 'egene_bonf', 'realized_r2', 'n_pheno']].round(3).to_string())
    print(coloc.round(3).to_string())
    print(go.round(3).to_string())


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
