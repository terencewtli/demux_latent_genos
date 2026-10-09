"""A08b: one chromosome of matched dosage matrices for the OneK1K pilot donors, for eQTL-anchored scoring.

Two genotype sources on the same donors and the same variants (common biallelic SNVs, GRCh38):
  latent   GLIMPSE2 ambient_b2m5 dosages imputed from souporcell clusters (A07f), clusters renamed to
           OneK1K donors with the A07g preimpute assign.tsv (clean clusters only: ambiguous == 0, r >= --min-r).
  onek1k   Zenodo 7619796 OneK1K.noGP.vcf.gz DS (array + Minimac4 imputation; the "reference" genotypes),
           with INFO/R2 and a typed flag (site in A07c onek1k.typed.b38).
Plus per-variant context used by A08c:
  n_target_100kb / dist_target   souporcell GL target sites (A07f chrN.target.bcf, union over pools) within
                                 +-100 kb, and distance to the nearest one (the "local site density" of A06 Level B)

Pools and souporcell versions come from POOLS below (pool 70 / 11 / 19 use the A07d rerun "_v2"; pool 70 v1
merged two donors and is never used). Missing pools are skipped with a message, so the script can be rerun as
more pools finish.

Outputs (results/A08_eqtl_anchor/dosage/):
  chrN.variants.tsv.gz   chrom pos ref alt raf info_glimpse r2_onek1k typed n_target_100kb dist_target
  chrN.latent.npy / chrN.onek1k.npy   float32 [variant x donor], same row order as variants
  chrN.donors.txt        donor order (OneK1K_<n>)

usage: python scripts/A08b_extract_dosages.py --chrom chr22 [--min-raf 0.01] [--min-r 0.6]
"""
import argparse
import os
import subprocess
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A07'
ONEK1K = f'{SCR}/truth/OneK1K.noGP.vcf.gz'
TYPED = f'{SCR}/truth/onek1k.typed.b38.vcf.gz'
OUT = f'{PROJ}/results/A08_eqtl_anchor/dosage'
POOLS = {'70': '_v2', '1': '', '55': '', '11': '_v2', '19': '_v2'}


def log(msg: str) -> None:
    print(msg, flush=True)


def run(cmd: List[str]) -> str:
    return subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout


def sample_order(path: str, samples: List[str]) -> List[str]:
    hdr = subprocess.run(['bcftools', 'view', '-h', '-s', ','.join(samples), path], check=True,
                         stdout=subprocess.PIPE, universal_newlines=True).stdout
    return hdr.splitlines()[-1].split('\t')[9:]


def pool_donors(pool: str, sfx: str, min_r: float) -> Dict[str, str]:
    a = pd.read_csv(f'{PROJ}/results/A07_onek1k/preimpute/pool{pool}{sfx}/assign.tsv', sep='\t')
    a = a[(a.ambiguous == 0) & (a.r >= min_r)]
    return dict(zip(a.cluster, a.donor))


def read_latent(bcf: str, cmap: Dict[str, str], min_raf: float) -> pd.DataFrame:
    samples = run(['bcftools', 'query', '-l', bcf]).split()
    keep = [s for s in samples if s in cmap]
    expr = f'TYPE="snp" && N_ALT=1 && INFO/RAF>={min_raf} && INFO/RAF<={1 - min_raf}'
    txt = run(['bcftools', 'query', '-s', ','.join(keep), '-i', expr,
               '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF\t%INFO/INFO[\t%DS]\n', bcf])
    cols = ['chrom', 'pos', 'ref', 'alt', 'raf', 'info_glimpse'] + [cmap[s] for s in sample_order(bcf, keep)]
    rows = [l.split('\t') for l in txt.splitlines()]
    d = pd.DataFrame(rows, columns=cols)
    d['pos'] = d.pos.astype(int)
    for c in cols[4:]:
        d[c] = d[c].astype(np.float32)
    return d


def read_onek1k(chrom: str, donors: List[str], min_raf: float) -> pd.DataFrame:
    expr = f'TYPE="snp" && N_ALT=1 && INFO/AF>={min_raf} && INFO/AF<={1 - min_raf}'
    txt = run(['bcftools', 'query', '-r', chrom.replace('chr', ''), '-s', ','.join(donors), '-i', expr,
               '-f', '%POS\t%REF\t%ALT\t%INFO/R2[\t%DS]\n', ONEK1K])
    order = sample_order(ONEK1K, donors)
    d = pd.DataFrame([l.split('\t') for l in txt.splitlines()], columns=['pos', 'ref', 'alt', 'r2_onek1k'] + order)
    d['pos'] = d.pos.astype(int)
    for c in ['r2_onek1k'] + order:
        d[c] = d[c].astype(np.float32)
    return d


def target_context(pos: np.ndarray, targets: np.ndarray, win: int = 100000) -> Tuple[np.ndarray, np.ndarray]:
    targets = np.sort(np.unique(targets))
    lo = np.searchsorted(targets, pos - win, side='left')
    hi = np.searchsorted(targets, pos + win, side='right')
    i = np.clip(np.searchsorted(targets, pos), 1, len(targets) - 1)
    dist = np.minimum(np.abs(pos - targets[i - 1]), np.abs(targets[i] - pos))
    return (hi - lo).astype(np.int32), dist.astype(np.int64)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--min-raf', type=float, default=0.01)
    ap.add_argument('--min-r', type=float, default=0.6)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    c = a.chrom

    lat, targets = None, []
    for pool, sfx in POOLS.items():
        bcf = f'{PROJ}/results/A07_onek1k/impute/pool{pool}{sfx}/{c}.imputed.bcf'
        assign = f'{PROJ}/results/A07_onek1k/preimpute/pool{pool}{sfx}/assign.tsv'
        if not (os.path.exists(bcf + '.csi') and os.path.exists(assign)):
            log(f'  pool {pool}{sfx}: not ready, skipped')
            continue
        cmap = pool_donors(pool, sfx, a.min_r)
        d = read_latent(bcf, cmap, a.min_raf)
        log(f'  pool {pool}{sfx}: {len(cmap)} clean clusters, {len(d)} common SNVs')
        tgt = f'{PROJ}/results/A07_onek1k/impute/pool{pool}{sfx}/{c}.target.bcf'
        targets.append(np.array(run(['bcftools', 'query', '-f', '%POS\n', tgt]).split(), dtype=np.int64))
        if lat is None:
            lat = d
        else:
            dup = [x for x in d.columns[6:] if x in lat.columns]
            if dup:
                log(f'  WARNING: donors already seen in another pool, dropped from pool {pool}: {dup}')
                d = d.drop(columns=dup)
            lat = lat.merge(d.drop(columns=['raf', 'info_glimpse']), on=['chrom', 'pos', 'ref', 'alt'], how='inner')
    if lat is None:
        raise SystemExit('no pools ready')
    donors = list(lat.columns[6:])
    log(f'{c}: {len(donors)} donors, {len(lat)} common SNVs imputed in every pool')

    ok = read_onek1k(c, donors, a.min_raf)
    missing = [x for x in donors if x not in ok.columns]
    if missing:
        log(f'  WARNING: donors not in the OneK1K VCF, dropped: {missing}')
        donors = [x for x in donors if x in ok.columns]
    m = lat[['chrom', 'pos', 'ref', 'alt', 'raf', 'info_glimpse'] + donors].merge(
        ok[['pos', 'ref', 'alt', 'r2_onek1k'] + donors], on=['pos', 'ref', 'alt'], suffixes=('', '_ok'))
    log(f'  matched to OneK1K (exact alleles): {len(m)}')

    typed = set(run(['bcftools', 'query', '-r', c, '-f', '%POS:%REF:%ALT\n', TYPED]).split())
    m['typed'] = (m.pos.astype(str) + ':' + m.ref + ':' + m.alt).isin(typed).astype(np.int8)
    m['n_target_100kb'], m['dist_target'] = target_context(m.pos.values, np.concatenate(targets))

    vcols = ['chrom', 'pos', 'ref', 'alt', 'raf', 'info_glimpse', 'r2_onek1k', 'typed', 'n_target_100kb', 'dist_target']
    m[vcols].to_csv(f'{OUT}/{c}.variants.tsv.gz', sep='\t', index=False)
    np.save(f'{OUT}/{c}.latent.npy', m[donors].values.astype(np.float32))
    np.save(f'{OUT}/{c}.onek1k.npy', m[[x + '_ok' for x in donors]].values.astype(np.float32))
    with open(f'{OUT}/{c}.donors.txt', 'w') as f:
        f.write('\n'.join(donors) + '\n')
    log(f'{c}: wrote {len(m)} variants x {len(donors)} donors; typed {m.typed.sum()}')


if __name__ == '__main__':
    main()
