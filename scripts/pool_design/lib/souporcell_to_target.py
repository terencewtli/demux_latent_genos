"""Convert souporcell cluster_genotypes.vcf into an imputation target VCF.

souporcell (consensus.py) writes one sample per cluster ("0".."k-1") with
FORMAT GT:AO:RO:T:E:GO:GN, where
  GO = natural-log genotype likelihoods, truncated to int, in the order [0/0, 1/1, 0/1]
  GN = log posteriors in the same order
  GT = argmax posterior if it exceeds 0.5, else ./.
  FILTER = BACKGROUND when the locus error weight exceeds 0.5

Output (uncompressed VCF on stdout; the caller bgzips it):
  samples renamed c0..c{k-1}
  FORMAT GT:AD:GL, with GL converted to log10 in VCF order [0/0, 0/1, 1/1] so
    GLIMPSE/Beagle can use it later (minimac4 reads only GT)
  INFO CR (call rate across clusters) and MONO (1 if the called genotypes are
    all hom-ref or all hom-alt), which A04b uses for the server-like QC
  only biallelic SNVs; BACKGROUND loci dropped
Per-site counts go to stderr.
"""
import argparse
import math
import sys
from typing import List, Optional, TextIO

LN10 = math.log(10)
HEADER_EXTRA = [
    '##INFO=<ID=CR,Number=1,Type=Float,Description="Fraction of clusters with a called GT">',
    '##INFO=<ID=MONO,Number=1,Type=Integer,Description="1 if called GTs are all 0/0 or all 1/1">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="souporcell hard call">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="souporcell RO,AO">',
    '##FORMAT=<ID=GL,Number=G,Type=Float,Description="log10 likelihoods from souporcell GO, reordered to 0/0,0/1,1/1">',
]


def convert_gl(go: str) -> str:
    vals: List[Optional[float]] = []
    for x in go.split(','):
        vals.append(None if x in ('NaN', 'nan', '.') else float(x) / LN10)
    if len(vals) != 3 or any(v is None for v in vals):
        return '.'
    hom_ref, hom_alt, het = vals
    return ','.join(f'{v:.3f}' for v in (hom_ref, het, hom_alt))


def run(vcf_in: TextIO, out: TextIO) -> None:
    n = {'records': 0, 'background': 0, 'not_snv': 0, 'written': 0}
    for line in vcf_in:
        if line.startswith('##'):
            if not line.startswith(('##FORMAT', '##INFO')):
                out.write(line)
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            k = len(f) - 9
            out.write('\n'.join(HEADER_EXTRA) + '\n')
            out.write('\t'.join(f[:9] + [f'c{i}' for i in range(k)]) + '\n')
            continue
        n['records'] += 1
        if 'BACKGROUND' in f[6]:
            n['background'] += 1
            continue
        ref, alt = f[3], f[4]
        if len(ref) != 1 or len(alt) != 1 or alt not in 'ACGT':
            n['not_snv'] += 1
            continue
        keys = f[8].split(':')
        i_gt, i_ao, i_ro, i_go = (keys.index(x) for x in ('GT', 'AO', 'RO', 'GO'))
        calls, gts = [], []
        for s in f[9:]:
            v = s.split(':')
            gt = v[i_gt]
            gts.append(gt)
            calls.append(f'{gt}:{v[i_ro]},{v[i_ao]}:{convert_gl(v[i_go])}')
        called = [g for g in gts if '.' not in g]
        cr = len(called) / len(gts)
        mono = int(len(called) == 0 or all(g == '0/0' for g in called) or all(g == '1/1' for g in called))
        out.write('\t'.join(f[:5] + ['.', 'PASS', f'CR={cr:.4f};MONO={mono}', 'GT:AD:GL'] + calls) + '\n')
        n['written'] += 1
    sys.stderr.write(' '.join(f'{k}={v}' for k, v in n.items()) + '\n')


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('vcf', help='souporcell cluster_genotypes.vcf')
    a = p.parse_args()
    with open(a.vcf) as fh:
        run(fh, sys.stdout)


if __name__ == '__main__':
    main()
