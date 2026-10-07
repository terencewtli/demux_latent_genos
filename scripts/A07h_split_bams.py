"""A07h step 1: split a pool's STARsolo BAM into one BAM per donor-like group, for read-based imputation
(GLIMPSE2 --bam-list, QUILT2), as in Kockelbergh 2026 (scTAPAS) but without known genotypes.

Two labelings of the same cells:
  soup    souporcell singlets, one BAM per cluster (c0, c1, ...): the genotype-free setting
  oracle  OneK1K's published per-cell donor labels (GEO Individual_Barcodes), one BAM per donor, named by the
          OneK1K VCF sample ID (convention resolved from A07g assign.tsv): the known-genotype setting of scTAPAS. oracle vs soup = the cost of demultiplexing
          without genotypes.
Reads kept: primary (no secondary / supplementary), MAPQ >= MIN_MAPQ (STAR unique = 255), valid CB and UB.
UMI dedup: one read per (CB, UB, start, strand) per group (STARsolo does not mark duplicates; GLIMPSE2 / QUILT
would otherwise count PCR copies of one molecule as independent evidence).
Output: $SCR/A07/h_bams/pool<P>/<mode>/<name>.bam(.bai), names.txt, bamlist.txt (GLIMPSE2: "<bam> <name>"),
        quilt_bamlist.txt (paths only), stats.tsv
usage: python scripts/A07h_split_bams.py --pool 70 --mode soup --chroms chr6,chr22
"""
import argparse
import glob
import os

import pandas as pd
import pysam

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
SCR = '/u/project/cluo_scratch/terencew/claude/latent_genos/A07'
GEO = f'{PROJ}/reference/onek1k_geo'
MIN_MAPQ = 255


def soup_labels(pool: str) -> dict:
    c = pd.read_csv(f'{SCR}/souporcell/pool{pool}/clusters.tsv', sep='\t')
    c = c[c.status == 'singlet']
    return dict(zip(c.barcode.str.replace(r'-1$', '', regex=True), 'c' + c.assignment.astype(str)))


def oracle_labels(pool: str) -> dict:
    """GEO 'Individual ID' is '<a>_<b>' and both OneK1K_<a> and OneK1K_<b> are VCF samples, so the ID alone does not
    say which is the donor. A07g resolves this per pool from the souporcell genotype match (assign.tsv); only that
    naming convention is taken from it, the per-cell labels are GEO's."""
    pools = pd.read_csv(f'{PROJ}/txt/onek1k_pilot_pools.txt', sep='\t', header=None, names=['pool', 'gsm', 'k'], dtype=str)
    gsm = pools.loc[pools.pool == pool, 'gsm'].item()
    g = pd.read_csv(glob.glob(f'{GEO}/{gsm}_*Individual_Barcodes.csv.gz')[0], dtype=str)
    assigned = set(pd.read_csv(f'{PROJ}/results/A07_onek1k/score/pool{pool}/assign.tsv', sep='\t').donor)
    ids = g['Individual ID'].unique()
    first = sum(f'OneK1K_{x.split("_")[0]}' in assigned for x in ids)
    second = sum(f'OneK1K_{x.split("_")[1]}' in assigned for x in ids)
    k = 0 if first >= second else 1
    print(f'GEO ids: {len(ids)}; matching A07g-assigned donors: first half {first}, second {second}; '
          f'using {"first" if k == 0 else "second"}', flush=True)
    g['name'] = 'OneK1K_' + g['Individual ID'].str.split('_').str[k]
    return dict(zip(g.Barcode.str.replace(r'-1$', '', regex=True), g.name))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--pool', required=True)
    p.add_argument('--mode', choices=['soup', 'oracle'], required=True)
    p.add_argument('--chroms', default='chr6,chr22')
    p.add_argument('--outdir', default=None, help='override (tests)')
    a = p.parse_args()
    out = a.outdir or f'{SCR}/h_bams/pool{a.pool}/{a.mode}'
    if os.path.exists(f'{out}/stats.tsv'):
        print(f'{out}/stats.tsv exists, skipping'); return
    os.makedirs(out, exist_ok=True)
    lab = soup_labels(a.pool) if a.mode == 'soup' else oracle_labels(a.pool)
    names = sorted(set(lab.values()))
    print(f'pool {a.pool} {a.mode}: {len(lab)} labelled cells, {len(names)} groups', flush=True)
    src = pysam.AlignmentFile(f'{SCR}/starsolo/pool{a.pool}/Aligned.sortedByCoord.out.bam')
    outs = {n: pysam.AlignmentFile(f'{out}/{n}.tmp.bam', 'wb', template=src) for n in names}
    kept = {n: 0 for n in names}
    seen = set()
    n_in = n_dup = n_nolab = 0
    for chrom in a.chroms.split(','):
        seen.clear()
        for r in src.fetch(chrom):
            n_in += 1
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < MIN_MAPQ:
                continue
            if not (r.has_tag('CB') and r.has_tag('UB')):
                continue
            cb, ub = r.get_tag('CB'), r.get_tag('UB')
            if cb == '-' or ub == '-':
                continue
            n = lab.get(cb)
            if n is None:
                n_nolab += 1
                continue
            key = (cb, ub, r.reference_start, r.is_reverse)
            if key in seen:
                n_dup += 1
                continue
            seen.add(key)
            outs[n].write(r)
            kept[n] += 1
        print(f'  {chrom}: {n_in:,} reads scanned, {sum(kept.values()):,} kept, {n_dup:,} UMI duplicates, '
              f'{n_nolab:,} from unlabelled cells', flush=True)
    for n, f in outs.items():
        f.close()
        os.replace(f'{out}/{n}.tmp.bam', f'{out}/{n}.bam')
        pysam.index(f'{out}/{n}.bam')
    with open(f'{out}/names.txt', 'w') as fh:
        fh.write('\n'.join(names) + '\n')
    with open(f'{out}/bamlist.txt', 'w') as fh:
        fh.write(''.join(f'{out}/{n}.bam {n}\n' for n in names))
    with open(f'{out}/quilt_bamlist.txt', 'w') as fh:
        fh.write(''.join(f'{out}/{n}.bam\n' for n in names))
    cells = pd.Series(list(lab.values())).value_counts()
    pd.DataFrame({'name': names, 'cells': [int(cells[n]) for n in names], 'reads': [kept[n] for n in names]}).to_csv(
        f'{out}/stats.tsv', sep='\t', index=False)   # written last: marks completion
    print(pd.read_csv(f'{out}/stats.tsv', sep='\t').to_string(index=False))


if __name__ == '__main__':
    main()
