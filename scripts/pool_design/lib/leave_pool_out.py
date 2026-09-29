"""Write the leave-pool-out sample lists for one pool.

The pool's donors are excluded together with their 1000G relatives from the
3,202-sample pedigree: parents, children, and siblings (anyone sharing a
non-zero parent). The 30x panel contains 602 trios, so leaving out only the
donors would still leak their haplotypes into the panel.

Writes <outdir>/exclude.txt and <outdir>/keep.txt (panel samples minus exclude).
"""
import argparse
from typing import Dict, List, Set, Tuple


def read_ped(path: str) -> Dict[str, Tuple[str, str]]:
    ped: Dict[str, Tuple[str, str]] = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.split()
            ped[f[1]] = (f[2], f[3])
    return ped


def relatives(donor: str, ped: Dict[str, Tuple[str, str]]) -> Set[str]:
    out: Set[str] = set()
    parents = {p for p in ped.get(donor, ('0', '0')) if p != '0'}
    out |= parents
    for s, (fa, mo) in ped.items():
        if donor in (fa, mo):
            out.add(s)
        if parents & {fa, mo}:
            out.add(s)
    out.discard(donor)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--donors', required=True, help='one donor ID per line')
    p.add_argument('--panel-samples', required=True, help='one panel sample ID per line')
    p.add_argument('--ped', required=True, help='20130606_g1k_3202_samples_ped_population.txt')
    p.add_argument('--outdir', required=True)
    a = p.parse_args()

    ped = read_ped(a.ped)
    donors: List[str] = [x.strip() for x in open(a.donors) if x.strip()]
    panel: List[str] = [x.strip() for x in open(a.panel_samples) if x.strip()]

    rel: Set[str] = set()
    for d in donors:
        rel |= relatives(d, ped)
    exclude = set(donors) | rel
    keep = [s for s in panel if s not in exclude]

    with open(f'{a.outdir}/exclude.txt', 'w') as fh:
        fh.write('\n'.join(sorted(exclude)) + '\n')
    with open(f'{a.outdir}/keep.txt', 'w') as fh:
        fh.write('\n'.join(keep) + '\n')
    print(f'donors={len(donors)} relatives_in_ped={len(rel)} '
          f'excluded_from_panel={len(exclude & set(panel))} keep={len(keep)}')


if __name__ == '__main__':
    main()
