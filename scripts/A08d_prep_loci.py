"""A08d prep: per fine-mapping locus, published OneK1K z-scores + genotype matrices from four LD sources.

Loci: A08a loci.tsv (round-1 lead per cell type x gene, |z| >= 5), one cell type per run; --max-loci random loci
(seed 1) so the run size is bounded and unbiased with respect to signal strength above the threshold.
For each locus (lead +- --win, 1000G MAF >= --min-maf, variants present in A08b for this chromosome):
  z.tsv           variant, z per ALT allele (OneK1K Spearman P -> z, sign of rho; A2 is the assessed allele,
                  flipped when A2 == REF), GRCh37 positions lifted to GRCh38
  geno_<arm>.tsv.gz   samples x variants ALT dosage, columns in z.tsv order:
    onek1k_all     every OneK1K sample (~1,098; the in-sample LD the z-scores came from: "gold")
    onek1k_pilot   OneK1K dosages, pilot donors only (A08b)        } same donors: the difference is imputation
    latent_pilot   souporcell -> GLIMPSE2 dosages, pilot donors (A08b) }
    kg_eur         1000G 30x EUR founders (out-of-sample reference panel LD)
Output: results/A08_eqtl_anchor/finemap/<ct>/<gene>/

usage: python scripts/A08d_prep_loci.py --ct cd4nc [--win 250000] [--max-loci 200] [--min-maf 0.05] [--chroms 22]
"""
import argparse
import os
import subprocess
from typing import Dict, List

import numpy as np
import pandas as pd
from pyliftover import LiftOver

PROJ = '/u/project/cluo/terencew/claude/project_ideas/latent_genos'
OUT = f'{PROJ}/results/A08_eqtl_anchor'
TABLES = f'{PROJ}/reference/onek1k/published_eqtl'
ONEK1K = '/u/project/cluo_scratch/terencew/claude/latent_genos/A07/truth/OneK1K.noGP.vcf.gz'
KG = '/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom/1000G.{c}.vcf.gz'
KG_CHR1 = f'{PROJ}/reference/1000G_30x/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
PED = '/u/project/cluo/terencew/demux_benchmark/pool_design/csv/1000G/meta/20130606_g1k_3202_samples_ped_population.txt'
CHAIN = '/u/project/cluo/terencew/reference/liftover/hg19ToHg38.over.chain.gz'
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
GT = {'0|0': 0, '0|1': 1, '1|0': 1, '1|1': 2, '0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2}


def log(msg: str) -> None:
    print(msg, flush=True)


def run(cmd: List[str]) -> str:
    return subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout


def read_z(ct: str, genes: set) -> pd.DataFrame:
    cols = ['GENE_ID', 'CHR', 'POS', 'A1', 'A2', 'SPEARMANS_RHO', 'P_VALUE', 'ROUND']
    keep = []
    for ch in pd.read_csv(f'{TABLES}/{ct}_eqtl_table.tsv.gz', sep='\t', usecols=cols, chunksize=2_000_000):
        keep.append(ch[(ch.ROUND == 1) & ch.GENE_ID.isin(genes)])
    z = pd.concat(keep, ignore_index=True)
    from scipy.stats import norm
    z['z'] = np.sign(z.SPEARMANS_RHO) * norm.isf(z.P_VALUE.clip(lower=1e-300) / 2)
    return z


def lift(z: pd.DataFrame) -> pd.DataFrame:
    lo = LiftOver(CHAIN)
    out = []
    for c, p in zip(z.CHR.astype(int), z.POS.astype(int)):
        r = lo.convert_coordinate(f'chr{c}', p - 1)
        out.append((r[0][0], r[0][1] + 1, r[0][2]) if r and len(r) == 1 else (None, -1, None))
    z['chrom'], z['pos'], z['strand'] = zip(*out)
    z = z[z.chrom == 'chr' + z.CHR.astype(int).astype(str)].copy()
    flip = z.strand == '-'
    z.loc[flip, 'A1'] = z.loc[flip, 'A1'].map(COMP)
    z.loc[flip, 'A2'] = z.loc[flip, 'A2'].map(COMP)
    return z


def kg_dosage(c: str, pos: np.ndarray, ref: np.ndarray, alt: np.ndarray, samples: List[str]) -> np.ndarray:
    path = KG_CHR1 if c == 'chr1' else KG.format(c=c)
    region = f'{c}:{pos.min()}-{pos.max()}'
    txt = run(['bcftools', 'query', '-r', region, '-s', ','.join(samples), '-i', 'TYPE="snp" && N_ALT=1',
               '-f', '%POS\t%REF\t%ALT[\t%GT]\n', path])
    idx = {f'{p}:{r}:{a}': i for i, (p, r, a) in enumerate(zip(pos, ref, alt))}
    X = np.full((len(samples), len(pos)), np.nan, dtype=np.float32)
    for line in txt.splitlines():
        f = line.split('\t')
        i = idx.get(f'{f[0]}:{f[1]}:{f[2]}')
        if i is not None:
            X[:, i] = [GT.get(g, np.nan) for g in f[3:]]
    return X


def onek1k_all(c: str, pos: np.ndarray, ref: np.ndarray, alt: np.ndarray) -> np.ndarray:
    txt = run(['bcftools', 'query', '-r', f'{c[3:]}:{pos.min()}-{pos.max()}', '-i', 'TYPE="snp" && N_ALT=1',
               '-f', '%POS\t%REF\t%ALT[\t%DS]\n', ONEK1K])
    idx = {f'{p}:{r}:{a}': i for i, (p, r, a) in enumerate(zip(pos, ref, alt))}
    n = len(run(['bcftools', 'query', '-l', ONEK1K]).split())
    X = np.full((n, len(pos)), np.nan, dtype=np.float32)
    for line in txt.splitlines():
        f = line.split('\t')
        i = idx.get(f'{f[0]}:{f[1]}:{f[2]}')
        if i is not None:
            X[:, i] = np.array(f[3:], dtype=np.float32)
    return X


def write_geno(path: str, X: np.ndarray, rows: List[str], cols: List[str]) -> None:
    pd.DataFrame(X, index=rows, columns=cols).to_csv(path, sep='\t', float_format='%.3f', compression='gzip')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--ct', required=True)
    ap.add_argument('--win', type=int, default=250000)
    ap.add_argument('--max-loci', type=int, default=200)
    ap.add_argument('--min-maf', type=float, default=0.05)
    ap.add_argument('--chroms', default='', help='comma-separated GRCh37 CHR values to restrict to (testing)')
    a = ap.parse_args()

    loci = pd.read_csv(f'{OUT}/loci.tsv', sep='\t')
    loci = loci[loci.CELL_ID == a.ct]
    if a.chroms:
        loci = loci[loci.CHR.astype(str).isin(a.chroms.split(','))]
    if len(loci) > a.max_loci:
        loci = loci.sample(a.max_loci, random_state=1)
    log(f'{a.ct}: {len(loci)} loci')
    z = lift(read_z(a.ct, set(loci.GENE_ID)))
    log(f'  round-1 z rows for these genes: {len(z)} (lifted)')
    ped = pd.read_csv(PED, sep=' ')
    eur = list(ped[(ped.Superpopulation == 'EUR') & (ped.FatherID == '0') & (ped.MotherID == '0')].SampleID)
    ok_samples = run(['bcftools', 'query', '-l', ONEK1K]).split()

    cache: Dict[str, tuple] = {}
    for _, L in loci.iterrows():
        lead_lift = z[(z.GENE_ID == L.GENE_ID) & (z.CHR == L.CHR) & (z.POS == L.POS)]
        if len(lead_lift) != 1:
            log(f'  {L.GENE_ID}: lead not lifted, skipped')
            continue
        c, lead_pos = lead_lift.chrom.iat[0], int(lead_lift.pos.iat[0])
        if c not in cache:
            p = f'{OUT}/dosage/{c}.variants.tsv.gz'
            if not os.path.exists(p):
                log(f'  {L.GENE_ID}: A08b {c} missing, skipped')
                continue
            cache = {c: (pd.read_csv(p, sep='\t'), np.load(f'{OUT}/dosage/{c}.latent.npy'),
                         np.load(f'{OUT}/dosage/{c}.onek1k.npy'), open(f'{OUT}/dosage/{c}.donors.txt').read().split())}
        v, lat, ok, donors = cache[c]
        zg = z[(z.GENE_ID == L.GENE_ID) & (abs(z.pos - lead_pos) <= a.win)]
        m = v.reset_index().merge(zg[['pos', 'A1', 'A2', 'z']], on='pos')
        same = (m.A2 == m.alt) & (m.A1 == m.ref)
        flip = (m.A2 == m.ref) & (m.A1 == m.alt)
        m = m[same | flip].copy()
        m.loc[flip[same | flip].values, 'z'] *= -1
        m = m[np.minimum(m.raf, 1 - m.raf) >= a.min_maf].drop_duplicates('pos').sort_values('pos')
        if len(m) < 20 or not (m.pos == lead_pos).any():
            log(f'  {L.GENE_ID}: {len(m)} variants, lead present {(m.pos == lead_pos).any()}; skipped')
            continue
        d = f'{OUT}/finemap/{a.ct}/{L.GENE_ID}'
        os.makedirs(d, exist_ok=True)
        vid = [f'{c}:{p}:{r}:{t}' for p, r, t in zip(m.pos, m.ref, m.alt)]
        pd.DataFrame({'variant': vid, 'z': m.z.values, 'is_lead': (m.pos == lead_pos).astype(int).values,
                      'n_target_100kb': m.n_target_100kb.values, 'typed': m.typed.values}).to_csv(
            f'{d}/z.tsv', sep='\t', index=False)
        ix = m['index'].values
        write_geno(f'{d}/geno_latent_pilot.tsv.gz', lat[ix].T, donors, vid)
        write_geno(f'{d}/geno_onek1k_pilot.tsv.gz', ok[ix].T, donors, vid)
        P, R, A = m.pos.values, m.ref.values, m.alt.values
        write_geno(f'{d}/geno_onek1k_all.tsv.gz', onek1k_all(c, P, R, A), ok_samples, vid)
        write_geno(f'{d}/geno_kg_eur.tsv.gz', kg_dosage(c, P, R, A, eur), eur, vid)
        log(f'  {L.GENE_ID} {c}:{lead_pos}: {len(m)} variants, max |z| {m.z.abs().max():.1f}')


if __name__ == '__main__':
    main()
