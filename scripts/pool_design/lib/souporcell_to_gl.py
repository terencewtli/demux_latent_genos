"""Recompute genotype likelihoods from souporcell cluster allele counts, for GLIMPSE2.

souporcell's own GO/GN fields are natural-log values truncated to integers, so
they are too coarse at 1-5 reads. This recomputes them from the per-cluster
AO/RO counts (sum over the cells assigned to the cluster).

Models (--ambient off = the plain binomial):
  binomial:  P(alt read | g) = e + (1 - 2e) * g/2        g = 0, 1, 2 alt alleles
  ambient:   P(alt read | g) = (1 - rho) * [e + (1 - 2e) * g/2] + rho * a_s
             rho = souporcell's ambient RNA fraction (ambient_rna.txt), a_s =
             alt fraction at the site pooled over all clusters (the ambient pool)
  log10 L(g) = AO * log10 p_g + RO * log10(1 - p_g)   (binomial coefficient cancels)
  PL = round(-10 * (log10 L(g) - max_g log10 L(g)))    (min 0; 0,0,0 when no reads)

Output (uncompressed VCF on stdout; the caller converts to BCF):
  samples c0..c{k-1} (same names as souporcell_to_target.py and A03b assign.tsv)
  FORMAT GT:DP:AD:PL with GT = ./. (GLIMPSE2 --input-gl reads only PL)
  header contigs chr1-22,X,Y (souporcell's header carries stale b37 contig names
    while its records are hg38 chr-prefixed)
  only biallelic SNVs; BACKGROUND loci dropped; optionally one chromosome
Per-site counts go to stderr.
"""
import argparse
import math
import re
import sys
from typing import List, Optional, TextIO, Tuple

CONTIGS = [f'chr{c}' for c in list(range(1, 23)) + ['X', 'Y']]
HEADER = ['##fileformat=VCFv4.2', '##FILTER=<ID=PASS,Description="All filters passed">']
HEADER += [f'##contig=<ID={c}>' for c in CONTIGS]
HEADER += [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Unset; GLIMPSE2 reads PL">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="souporcell cluster depth RO+AO">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="souporcell cluster RO,AO">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods recomputed from RO/AO">',
]


def read_ambient(path: str) -> float:
    txt = open(path).read()
    m = re.search(r'([0-9.]+)%', txt)
    if not m:
        raise ValueError(f'cannot parse ambient fraction from {path}: {txt!r}')
    return float(m.group(1)) / 100


def alt_probs(e: float, rho: float, a_site: float) -> Tuple[float, float, float]:
    out: List[float] = []
    for g in (0, 1, 2):
        p = (1 - rho) * (e + (1 - 2 * e) * g / 2) + rho * a_site
        out.append(min(max(p, 1e-6), 1 - 1e-6))
    return out[0], out[1], out[2]


def phred_pl(ao: int, ro: int, probs: Tuple[float, float, float]) -> str:
    if ao + ro == 0:
        return '0,0,0'
    ll = [ao * math.log10(p) + ro * math.log10(1 - p) for p in probs]
    top = max(ll)
    return ','.join(str(int(round(-10 * (x - top)))) for x in ll)


def run(vcf_in: TextIO, out: TextIO, chrom: Optional[str], e: float, rho: float) -> None:
    n = {'records': 0, 'other_chrom': 0, 'background': 0, 'not_snv': 0, 'written': 0}
    for line in vcf_in:
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            k = len(f) - 9
            out.write('\n'.join(HEADER) + '\n')
            out.write('\t'.join(f[:9] + [f'c{i}' for i in range(k)]) + '\n')
            continue
        n['records'] += 1
        if chrom is not None and f[0] != chrom:
            n['other_chrom'] += 1
            continue
        if 'BACKGROUND' in f[6]:
            n['background'] += 1
            continue
        ref, alt = f[3], f[4]
        if len(ref) != 1 or len(alt) != 1 or alt not in 'ACGT':
            n['not_snv'] += 1
            continue
        keys = f[8].split(':')
        i_ao, i_ro = keys.index('AO'), keys.index('RO')
        counts = []
        for s in f[9:]:
            v = s.split(':')
            counts.append((int(v[i_ao]), int(v[i_ro])))
        tot = sum(a + r for a, r in counts)
        a_site = sum(a for a, _ in counts) / tot if tot else 0.0
        probs = alt_probs(e, rho, a_site)
        calls = [f'./.:{a + r}:{r},{a}:{phred_pl(a, r, probs)}' for a, r in counts]
        out.write('\t'.join([f[0], f[1], f[2], ref, alt, '.', 'PASS', '.', 'GT:DP:AD:PL'] + calls) + '\n')
        n['written'] += 1
    sys.stderr.write(' '.join(f'{k}={v}' for k, v in n.items()) + f' error={e} rho={rho:.4f}\n')


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('vcf', help='souporcell cluster_genotypes.vcf')
    p.add_argument('--chrom', default=None, help='keep only this chromosome (e.g. chr20)')
    p.add_argument('--error', type=float, default=0.01, help='per-read base error e')
    p.add_argument('--ambient', default=None, help='souporcell ambient_rna.txt; enables the ambient model')
    a = p.parse_args()
    rho = read_ambient(a.ambient) if a.ambient else 0.0
    with open(a.vcf) as fh:
        run(fh, sys.stdout, a.chrom, a.error, rho)


if __name__ == '__main__':
    main()
