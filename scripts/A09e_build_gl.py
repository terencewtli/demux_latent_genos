#!/usr/bin/env python
"""A09e: donor-level genotype-likelihood target for GLIMPSE2 from one or more souporcell runs.

A target is <pools>_<arm>: pools = one pool (ys3a..ys3j) or 'all'; arm = gex | atac | multi (gex + atac).
Each run's clusters are renamed to donors with the A09d cluster map. Per run, PLs come from the cluster AO/RO
with the same ambient model as lib/souporcell_to_gl.py (that run's rho and pooled site alt fraction); log10
likelihoods are then summed over runs per donor (reads from different libraries are independent).

Runs whose genotype-based cluster map disagrees with the ambimux majority (A09d run_summary.tsv,
map_one_to_one_matches_ambi_majority False; ys3c GEX merges C29 + C39) are skipped.

Output: uncompressed VCF on stdout, samples C29 C37 C38 C39, FORMAT GT:DP:AD:PL, biallelic SNVs, one chromosome.
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict

import pandas as pd

sys.path.insert(0, '/u/project/cluo/terencew/claude/project_ideas/pool_design/scripts/latent_genos/lib')
from souporcell_to_gl import HEADER, alt_probs, read_ambient  # noqa: E402

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A09/souporcell_v2'
POOLS = [f'ys3{c}' for c in 'abcdefghij']
DONORS = ['C29', 'C37', 'C38', 'C39']


def runs_for(target: str) -> list[tuple[str, str]]:
    pools, arm = target.rsplit('_', 1)
    pl = POOLS if pools == 'all' else [pools]
    mods = ['gex', 'atac'] if arm == 'multi' else [arm]
    return [(p, m) for p in pl for m in mods]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--target', required=True)
    p.add_argument('--chrom', required=True)
    p.add_argument('--error', type=float, default=0.01)
    a = p.parse_args()

    cmap = pd.read_csv(f'{PROJ}/results/A09_multiome/match/cluster_map.tsv', sep='\t')
    rs = pd.read_csv(f'{PROJ}/results/A09_multiome/match/run_summary.tsv', sep='\t')
    bad = set(zip(rs.loc[~rs.map_one_to_one_matches_ambi_majority, 'pool'],
                  rs.loc[~rs.map_one_to_one_matches_ambi_majority, 'mod']))
    # (chrom, pos, ref, alt) -> donor -> [ll0, ll1, ll2, ro, ao]
    acc: dict = defaultdict(lambda: {d: [0.0, 0.0, 0.0, 0, 0] for d in DONORS})
    n = defaultdict(int)
    for pool, mod in runs_for(a.target):
        if (pool, mod) in bad:
            n[f'excluded_{pool}_{mod}'] = 1
            continue
        d_of = cmap[(cmap['pool'] == pool) & (cmap['mod'] == mod)].set_index('cluster')['donor'].to_dict()
        rho = read_ambient(f'{SCR}/{pool}_{mod}/ambient_rna.txt')
        with open(f'{SCR}/{pool}_{mod}/cluster_genotypes.vcf') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t', 9)
                if f[0] != a.chrom:
                    continue
                f = line.rstrip('\n').split('\t')
                if 'BACKGROUND' in f[6] or len(f[3]) != 1 or len(f[4]) != 1 or f[4] not in 'ACGT':
                    n['skipped'] += 1
                    continue
                keys = f[8].split(':')
                i_ao, i_ro = keys.index('AO'), keys.index('RO')
                counts = []
                for s in f[9:]:
                    v = s.split(':')
                    counts.append((int(v[i_ao]), int(v[i_ro])))
                tot = sum(x + y for x, y in counts)
                if tot == 0:
                    continue
                probs = alt_probs(a.error, rho, sum(x for x, _ in counts) / tot)
                rec = acc[(f[0], int(f[1]), f[3], f[4])]
                for ci, (ao, ro) in enumerate(counts):
                    if ao + ro == 0:
                        continue
                    r = rec[d_of[ci]]
                    for g in range(3):
                        r[g] += ao * math.log10(probs[g]) + ro * math.log10(1 - probs[g])
                    r[3] += ro
                    r[4] += ao
                n[f'records_{pool}_{mod}'] += 1

    out = sys.stdout
    out.write('\n'.join(HEADER) + '\n')
    out.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + DONORS) + '\n')
    for (c, pos, ref, alt) in sorted(acc, key=lambda k: (k[1], k[2], k[3])):
        calls = []
        for d in DONORS:
            ll0, ll1, ll2, ro, ao = acc[(c, pos, ref, alt)][d]
            if ro + ao == 0:
                pl = '0,0,0'
            else:
                top = max(ll0, ll1, ll2)
                pl = ','.join(str(int(round(-10 * (x - top)))) for x in (ll0, ll1, ll2))
            calls.append(f'./.:{ro + ao}:{ro},{ao}:{pl}')
        out.write('\t'.join([c, str(pos), '.', ref, alt, '.', 'PASS', '.', 'GT:DP:AD:PL'] + calls) + '\n')
    n['written'] = len(acc)
    sys.stderr.write(' '.join(f'{k}={v}' for k, v in sorted(n.items())) + '\n')


if __name__ == '__main__':
    main()
