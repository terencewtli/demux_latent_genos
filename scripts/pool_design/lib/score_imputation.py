"""Score imputed dosages against true 1000G 30x genotypes.

Only sites that were NOT in the imputation input (--typed) are scored, and only
biallelic SNVs, matched on CHROM/POS/REF/ALT. Bins come from the panel's AF
over all 3,202 samples.

private: the variant's alt alleles in the 3,202-sample panel are all carried by
the scored samples (AC_panel == AC_samples > 0). If those samples are also in
the imputation reference, such variants are recovered almost perfectly; if not,
they should be very hard to impute. This is the leakage signal.

Writes:
  <out>.sites.tsv.gz  per-site: id, af, ac_panel, ac_samples, private, r2
  <out>.bins.tsv      per-bin aggregate r2 (squared Pearson over all
                      site x sample dosages in the bin) and the mean per-site r2
"""
import argparse
import subprocess
from typing import Dict, Iterator, List, Set, Tuple

import numpy as np
import pandas as pd

MAF_BINS = [0, 0.001, 0.005, 0.01, 0.05, 0.5000001]
LABELS = ['<0.1%', '0.1-0.5%', '0.5-1%', '1-5%', '5-50%']


def query(cmd: List[str], view: List[str] = None) -> Iterator[List[str]]:
    """bcftools query, optionally behind a bcftools view filter. bcftools 1.11 query has no -v / -m / -M (its -v is
    --vcf-list), so site filters go through view. Fails loudly on a non-zero exit (the old version silently
    returned 0 sites; fixed 2026-09-30)."""
    if view:
        v = subprocess.Popen(view + ['-Ou'], stdout=subprocess.PIPE)
        p = subprocess.Popen(cmd + ['-'], stdin=v.stdout, stdout=subprocess.PIPE, text=True)
        v.stdout.close()
    else:
        v, p = None, subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert p.stdout is not None
    for line in p.stdout:
        yield line.rstrip('\n').split('\t')
    for proc in (p, v):
        if proc is not None and proc.wait() != 0:
            raise RuntimeError(f'bcftools failed ({proc.returncode}): {" ".join(view if proc is v else cmd)}')


def gt_to_dosage(gt: str) -> float:
    if '.' in gt:
        return np.nan
    return float(gt.count('1'))


def read_ids(path: str, region: str) -> Set[str]:
    return {f'{r[0]}:{r[1]}:{r[2]}:{r[3]}' for r in
            query(['bcftools', 'query', '-r', region, '-f', '%CHROM\t%POS\t%REF\t%ALT\n', path])}


def read_dose(path: str, samples: List[str], region: str, typed: Set[str]) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {}
    view = ['bcftools', 'view', '-r', region, '-s', ','.join(samples), '-v', 'snps', path]
    cmd = ['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n']
    for r in query(cmd, view):
        vid = f'{r[0]}:{r[1]}:{r[2]}:{r[3]}'
        if vid in typed or len(r[2]) != 1 or len(r[3]) != 1:
            continue
        out[vid] = np.array([float(x) if x != '.' else np.nan for x in r[4:]])
    return out


def r2(x: np.ndarray, y: np.ndarray) -> float:
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if len(x) < 3 or x.std() == 0 or y.std() == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1] ** 2)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--dose', required=True)
    p.add_argument('--typed', required=True, help='imputation input VCF; its sites are not scored')
    p.add_argument('--truth', required=True, help='1000G 30x panel VCF for this chromosome')
    p.add_argument('--samples', required=True)
    p.add_argument('--region', required=True)
    p.add_argument('--out', required=True, help='output prefix')
    a = p.parse_args()

    samples = [x.strip() for x in open(a.samples) if x.strip()]
    typed = read_ids(a.typed, a.region)
    dose = read_dose(a.dose, samples, a.region, typed)
    print(f'imputed untyped SNVs: {len(dose)}; typed sites excluded: {len(typed)}')

    rows: List[Tuple] = []
    tru: Dict[str, np.ndarray] = {}
    view = ['bcftools', 'view', '-r', a.region, '-s', ','.join(samples), '-m2', '-M2', '-v', 'snps', a.truth]
    cmd = ['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC[\t%GT]\n']
    for r in query(cmd, view):
        vid = f'{r[0]}:{r[1]}:{r[2]}:{r[3]}'
        if vid not in dose:
            continue
        t = np.array([gt_to_dosage(g) for g in r[6:]])
        ac_s = int(np.nansum(t))
        ac_p = int(r[5])
        tru[vid] = t
        rows.append((vid, float(r[4]), ac_p, ac_s, int(ac_s > 0 and ac_s == ac_p), r2(t, dose[vid])))

    df = pd.DataFrame(rows, columns=['id', 'af', 'ac_panel', 'ac_samples', 'private', 'r2'])
    maf = np.minimum(df.af, 1 - df.af)
    df['bin'] = pd.cut(maf, MAF_BINS, labels=LABELS, right=False)
    df.to_csv(f'{a.out}.sites.tsv.gz', sep='\t', index=False)

    out = []
    groups = [('maf_' + str(b), df[df.bin == b]) for b in LABELS]
    groups += [('carried_by_samples', df[df.ac_samples > 0]),
               ('private_to_samples', df[df.private == 1]),
               ('private_singleton', df[(df.private == 1) & (df.ac_panel == 1)])]
    for name, g in groups:
        if len(g) == 0:
            out.append((name, 0, np.nan, np.nan))
            continue
        x = np.concatenate([tru[v] for v in g.id])
        y = np.concatenate([dose[v] for v in g.id])
        out.append((name, len(g), r2(x, y), g.r2.mean()))
    res = pd.DataFrame(out, columns=['group', 'n_sites', 'aggregate_r2', 'mean_site_r2'])
    res.to_csv(f'{a.out}.bins.tsv', sep='\t', index=False)
    print(res.to_string(index=False))


if __name__ == '__main__':
    main()
