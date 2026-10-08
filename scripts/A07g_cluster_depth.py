import sys, numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
def one(arg):
    label, d = arg
    dps, chrom = [], []
    with open(f'{d}/cluster_genotypes.vcf') as f:
        for l in f:
            if l[0] == '#': continue
            x = l.rstrip('\n').split('\t')
            dps.append([int(s.split(':')[1]) + int(s.split(':')[2]) for s in x[9:]])
    D = np.array(dps)
    cl = pd.read_csv(f'{d}/clusters.tsv', sep='\t')
    sing = cl[cl.status == 'singlet'].assignment.astype(int).value_counts()
    rows = []
    for k in range(D.shape[1]):
        dp = D[:, k]; cov = dp[dp > 0]
        rows.append({'cells': sing.get(k, 0), 'sites_cov': (dp > 0).sum(), 'sites_ge5': (dp >= 5).sum(),
                     'sites_ge10': (dp >= 10).sum(), 'median_dp_cov': np.median(cov), 'mean_dp_cov': cov.mean(),
                     'frac_dp1_2': ((cov <= 2).sum()) / len(cov)})
    R = pd.DataFrame(rows)
    R = R[R.cells >= 100]
    m = R.median()
    return {'run': label, 'n_sites': D.shape[0], 'clusters_used': len(R), 'cells_per_donor': m.cells,
            'sites_cov_per_donor': m.sites_cov, 'sites_ge5': m.sites_ge5, 'sites_ge10': m.sites_ge10,
            'median_dp_at_cov': m.median_dp_cov, 'mean_dp_at_cov': m.mean_dp_cov, 'frac_cov_sites_1_2_reads': m.frac_dp1_2}
args = [l.rstrip('\n').split('\t') for l in open(sys.argv[1])]
with ProcessPoolExecutor(8) as ex:
    out = pd.DataFrame(list(ex.map(one, args)))
pd.set_option('display.width', 250)
print(out.round(2).to_string(index=False))
out.to_csv(sys.argv[2], sep='\t', index=False)
