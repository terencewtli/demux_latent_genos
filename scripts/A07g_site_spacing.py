import gzip, sys, numpy as np, pandas as pd
def load(p):
    op = gzip.open if p.endswith('gz') else open
    rows = []
    with op(p, 'rt') as f:
        for l in f:
            if l[0] == '#': continue
            c, pos = l.split('\t', 2)[:2]
            c = c if c.startswith('chr') else 'chr' + c
            rows.append((c, int(pos)))
    d = pd.DataFrame(rows, columns=['chrom', 'pos'])
    return d[d.chrom.isin([f'chr{i}' for i in range(1, 23)])].drop_duplicates()
def stats(name, d):
    gaps = np.concatenate([np.diff(np.sort(g.pos.values)) for _, g in d.groupby('chrom')])
    out = {'set': name, 'n': len(d), 'gap_median_kb': np.median(gaps)/1e3, 'gap_p90_kb': np.percentile(gaps, 90)/1e3,
           'gap_p99_kb': np.percentile(gaps, 99)/1e3, 'n_gaps_gt_500kb': int((gaps > 5e5).sum()),
           'Mb_in_gaps_gt_250kb': gaps[gaps > 2.5e5].sum()/1e6}
    for w in [1e5, 1e6]:
        occ = set(zip(d.chrom, d.pos // int(w)))
        tot = 0
        for c, g in d.groupby('chrom'):
            tot += g.pos.max() // int(w) + 1
        out[f'frac_{int(w/1e3)}kb_bins_hit'] = len(occ) / tot
    return out
soup = load(sys.argv[1]); arr = load(sys.argv[2])
R = pd.DataFrame([stats('souporcell_pool1', soup), stats('array_typed', arr)])
print(R.round(3).T.to_string())
