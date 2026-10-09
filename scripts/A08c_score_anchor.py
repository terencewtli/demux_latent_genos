"""A08c: eQTL-anchored accuracy and LD fidelity of latent (souporcell -> GLIMPSE2) genotypes vs OneK1K genotypes.

Inputs: A08b dosage matrices (results/A08_eqtl_anchor/dosage/chrN.*) and A08a leads.tsv (OneK1K published
independent cis-eQTL leads, GRCh37; lifted here to GRCh38 with pyliftover and matched on position + alleles).

Per variant (all chromosomes present): r2 between latent and OneK1K dosage across the pilot donors, annotated with
  MAF (1000G), typed, OneK1K R2, souporcell target density / distance (A08b), lead status and distance to the nearest
  lead. Scored variants: OneK1K R2 >= --min-ok-r2 (OneK1K's own imputation error otherwise dominates).
Summaries:
  summary_dist_target.tsv   r2 by distance to the nearest souporcell GL site (the A06 Level B mechanism)
  summary_lead.tsv          r2 at lead SNPs vs by distance-to-lead bin, x distance-to-target bin
  summary_lead_ct.tsv       r2 at lead SNPs per cell type and conditional round
  donor_r2.tsv              per-donor r2 over typed sites and over lead SNPs
  lead_ld.tsv               LD fidelity per unique lead: within +-ld_win, variants with MAF >= 5%;
                            r(lead, v) from latent vs from OneK1K dosages on the same donors:
                            corr of r2 vectors, mean |delta r2|, Jaccard of the r2 >= 0.8 proxy sets
  variant_r2.tsv.gz         the per-variant table
Note: with ~60 donors, per-variant r2 has a chance floor of ~1/(n-1) and a wide spread; compare groups, not values.

usage: python scripts/A08c_score_anchor.py [--min-ok-r2 0.8] [--ld-win 250000]
"""
import argparse
import glob
import os
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from pyliftover import LiftOver

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
OUT = f'{PROJ}/results/A08_eqtl_anchor'
CHAIN = '/u/project/cluo/terencew/reference/liftover/hg19ToHg38.over.chain.gz'
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def log(msg: str) -> None:
    print(msg, flush=True)


def rowr(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = a - a.mean(1, keepdims=True)
    b = b - b.mean(1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))


def lift_leads(path: str) -> pd.DataFrame:
    L = pd.read_csv(path, sep='\t')
    lo = LiftOver(CHAIN)
    pos38, chr38, strand = [], [], []
    for c, p in zip(L.CHR.astype(str), L.POS.astype(int)):
        r = lo.convert_coordinate(f'chr{c}', p - 1)
        if r and len(r) == 1:
            chr38.append(r[0][0]); pos38.append(r[0][1] + 1); strand.append(r[0][2])
        else:
            chr38.append(None); pos38.append(-1); strand.append(None)
    L['chrom'], L['pos'], L['strand'] = chr38, pos38, strand
    L = L[L.chrom.notna() & (L.chrom == 'chr' + L.CHR.astype(str))].copy()
    flip = L.strand == '-'
    L.loc[flip, 'A1'] = L.loc[flip, 'A1'].map(COMP)
    L.loc[flip, 'A2'] = L.loc[flip, 'A2'].map(COMP)
    log(f'leads: {len(L)} lifted to GRCh38 (same chromosome), {flip.sum()} on the minus strand')
    return L


def match_leads(L: pd.DataFrame, v: pd.DataFrame) -> pd.DataFrame:
    k = v[['chrom', 'pos', 'ref', 'alt']].reset_index().rename(columns={'index': 'vidx'})
    m = L.merge(k, on=['chrom', 'pos'])
    ok = ((m.A1 == m.ref) & (m.A2 == m.alt)) | ((m.A1 == m.alt) & (m.A2 == m.ref))
    return m[ok]


def nearest(pos: np.ndarray, anchors: np.ndarray) -> np.ndarray:
    if len(anchors) == 0:
        return np.full(len(pos), 10 ** 9)
    anchors = np.sort(np.unique(anchors))
    i = np.clip(np.searchsorted(anchors, pos), 1, max(len(anchors) - 1, 1))
    lo = anchors[np.maximum(i - 1, 0)]
    hi = anchors[np.minimum(i, len(anchors) - 1)]
    return np.minimum(np.abs(pos - lo), np.abs(hi - pos))


def ld_fidelity(v: pd.DataFrame, lat: np.ndarray, ok: np.ndarray, lead_idx: np.ndarray, win: int) -> pd.DataFrame:
    maf = np.minimum(v.raf.values, 1 - v.raf.values)
    pos = v.pos.values
    rows = []
    for i in np.unique(lead_idx):
        lo, hi = np.searchsorted(pos, pos[i] - win), np.searchsorted(pos, pos[i] + win, side='right')
        j = np.arange(lo, hi)
        j = j[(maf[j] >= 0.05) & (j != i)]
        if len(j) < 10:
            continue
        r_ok = rowr(np.repeat(ok[[i]], len(j), 0), ok[j]) ** 2
        r_lat = rowr(np.repeat(lat[[i]], len(j), 0), lat[j]) ** 2
        good = np.isfinite(r_ok) & np.isfinite(r_lat)
        if good.sum() < 10:
            continue
        r_ok, r_lat = r_ok[good], r_lat[good]
        p_ok, p_lat = r_ok >= 0.8, r_lat >= 0.8
        union = (p_ok | p_lat).sum()
        rows.append({'chrom': v.chrom.iat[i], 'pos': pos[i], 'n_vars': int(good.sum()),
                     'corr_r2': np.corrcoef(r_ok, r_lat)[0, 1], 'mean_abs_dr2': np.abs(r_ok - r_lat).mean(),
                     'n_proxy_ok': int(p_ok.sum()), 'n_proxy_lat': int(p_lat.sum()),
                     'proxy_jaccard': (p_ok & p_lat).sum() / union if union else np.nan,
                     'n_target_100kb': int(v.n_target_100kb.iat[i])})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-ok-r2', type=float, default=0.8)
    ap.add_argument('--ld-win', type=int, default=250000)
    a = ap.parse_args()
    out = f'{OUT}/score'
    os.makedirs(out, exist_ok=True)
    leads = lift_leads(f'{OUT}/leads.tsv')

    V, LM, donor_rows, LD = [], [], [], []
    chroms = sorted(glob.glob(f'{OUT}/dosage/chr*.variants.tsv.gz'),
                    key=lambda p: int(os.path.basename(p).split('.')[0][3:]))
    for p in chroms:
        c = os.path.basename(p).split('.')[0]
        v = pd.read_csv(p, sep='\t')
        lat = np.load(f'{OUT}/dosage/{c}.latent.npy')
        ok = np.load(f'{OUT}/dosage/{c}.onek1k.npy')
        donors = open(f'{OUT}/dosage/{c}.donors.txt').read().split()
        v['r2'] = rowr(lat, ok) ** 2
        lm = match_leads(leads[leads.chrom == c], v)
        v['is_lead'] = v.index.isin(lm.vidx).astype(np.int8)
        v['dist_lead'] = nearest(v.pos.values, v.pos.values[lm.vidx.unique()])
        lm = lm.assign(r2=v.r2.values[lm.vidx], dist_target=v.dist_target.values[lm.vidx],
                       r2_onek1k=v.r2_onek1k.values[lm.vidx])
        LM.append(lm.drop(columns=['vidx']))
        t = v.typed.values == 1
        li = lm.vidx.unique()
        for j, d in enumerate(donors):
            donor_rows.append({'chrom': c, 'donor': d, 'n_typed': int(t.sum()), 'n_lead': len(li),
                               'ss_typed': np.corrcoef(lat[t, j], ok[t, j])[0, 1] ** 2 * t.sum() if t.sum() > 2 else np.nan,
                               'ss_lead': np.corrcoef(lat[li, j], ok[li, j])[0, 1] ** 2 * len(li) if len(li) > 2 else np.nan})
        keep = v.r2_onek1k.values >= a.min_ok_r2
        LD.append(ld_fidelity(v, lat, ok, np.intersect1d(li, np.where(keep)[0]), a.ld_win))
        V.append(v[keep])
        log(f'{c}: {len(v)} variants, {keep.sum()} scored, {len(li)} lead SNPs matched, {len(donors)} donors')

    V = pd.concat(V, ignore_index=True)
    V['maf'] = np.minimum(V.raf, 1 - V.raf)
    V['maf_bin'] = pd.cut(V.maf, [0, 0.05, 0.2, 0.5])
    V['dist_target_bin'] = pd.cut(V.dist_target, [-1, 0, 5e3, 2e4, 1e5, 1e10],
                                  labels=['0', '<5kb', '5-20kb', '20-100kb', '>100kb'])
    V['dist_lead_bin'] = pd.cut(V.dist_lead, [-1, 0, 1e4, 1e5, 1e6, 1e10],
                                labels=['lead', '<10kb', '10-100kb', '100kb-1Mb', '>1Mb'])
    agg = {'r2': ['size', 'median', 'mean']}
    V.groupby(['dist_target_bin', 'maf_bin']).agg(agg).to_csv(f'{out}/summary_dist_target.tsv', sep='\t')
    V.groupby(['dist_lead_bin', 'dist_target_bin']).agg(agg).to_csv(f'{out}/summary_lead.tsv', sep='\t')
    LM = pd.concat(LM, ignore_index=True)
    LM = LM[LM.r2_onek1k >= a.min_ok_r2]
    LM.groupby(['CELL_ID', 'ROUND']).agg(agg).to_csv(f'{out}/summary_lead_ct.tsv', sep='\t')
    D = pd.DataFrame(donor_rows).groupby('donor')[['n_typed', 'n_lead', 'ss_typed', 'ss_lead']].sum()
    D['r2_typed'] = D.ss_typed / D.n_typed
    D['r2_lead'] = D.ss_lead / D.n_lead
    D[['n_typed', 'n_lead', 'r2_typed', 'r2_lead']].to_csv(f'{out}/donor_r2.tsv', sep='\t')
    LD = pd.concat(LD, ignore_index=True)
    LD.to_csv(f'{out}/lead_ld.tsv', sep='\t', index=False)
    V.drop(columns=['maf_bin', 'dist_target_bin', 'dist_lead_bin']).to_csv(
        f'{out}/variant_r2.tsv.gz', sep='\t', index=False, float_format='%.4g')

    log(f'scored variants {len(V)}; median r2 {V.r2.median():.3f}')
    log(V.groupby('dist_target_bin').r2.agg(['size', 'median', 'mean']).round(3).to_string())
    log(V.groupby('dist_lead_bin').r2.agg(['size', 'median', 'mean']).round(3).to_string())
    log(f'per-donor r2: typed median {D.r2_typed.median():.3f}, lead median {D.r2_lead.median():.3f}')
    log(f'LD fidelity over {len(LD)} leads: median corr_r2 {LD.corr_r2.median():.3f}, '
        f'median proxy Jaccard {LD.proxy_jaccard.median():.3f}')


if __name__ == '__main__':
    main()
