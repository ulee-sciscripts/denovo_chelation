#!/usr/bin/env python3
"""Diagnose the length-matched proteome baseline + do a LENGTH-STRATIFIED
H-x-H incidence comparison (intergenic ORFs vs proteome vs neORFs), so length
is controlled within bins rather than by matching to a single median.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, csv
import numpy as np
CODE = os.path.join(_HERE, "..")
sys.path.insert(0, CODE)
import denovo_hisspacing_analysis as A
from hs_common import read_fasta
from scipy import stats

TMP = sys.argv[1]; A.SCRATCH = TMP
GENOME = f"{_HERE}/../../data/analysis/ref/dm6_genome.fa"
INTER = os.path.join(_HERE, "..", "..", "data", "intermediate")
CHRS = ["2L", "2R", "3L", "3R"]
B2I, LUT, H, STAR = A.BASE2INT, A.AA_LUT, A.ORD_H, A.ORD_STAR
NPY = f"{_HERE}/data/interg_LH.npy"

import os
if os.path.exists(NPY):
    LH = np.load(NPY); iL, iH = LH[0], LH[1]
    print(f"[loaded cached intergenic scan: {iL.size:,} ORFs]")
else:
    seqs = {}
    for hdr, s in read_fasta(GENOME):
        c = hdr.split()[0]
        if c in CHRS: seqs[c] = s.upper()
    genes, reps = A.load_intervals()
    def scan_clean(seqb):
        out = []
        x = B2I[seqb]; xr = np.where(x == 4, 4, 3 - x).astype(np.uint8)[::-1]
        for src in (x, xr):
            n = src.size
            for fr in range(3):
                m = (n - fr) // 3
                if m < 10: continue
                cod = src[fr:fr + 3 * m].reshape(m, 3).astype(np.int64)
                valid = (cod < 4).all(axis=1)
                idx = np.where(valid, cod[:, 0]*16 + cod[:, 1]*4 + cod[:, 2], 0)
                aa = np.where(valid, LUT[idx], STAR).astype(np.uint8)
                stop = np.flatnonzero(aa == STAR); bnd = np.concatenate(([-1], stop, [aa.size]))
                for k in range(bnd.size - 1):
                    s0, e0 = bnd[k]+1, bnd[k+1]; Ln = e0 - s0
                    if Ln < 10: continue
                    hm = (aa[s0:e0] == H)
                    out.append((Ln, int(np.count_nonzero(hm[:-2] & hm[2:])) if Ln > 2 else 0))
        return out
    orfs = []
    for c in CHRS:
        seqb = np.frombuffer(seqs[c].encode(), np.uint8); Lc = len(seqs[c])
        clean = (~A.mask_from_intervals(genes[c], Lc)) & (~A.mask_from_intervals(reps[c], Lc))
        for s, e in A.runs_of_true(clean):
            if e - s >= 30: orfs.extend(scan_clean(seqb[s:e]))
    iL = np.array([o[0] for o in orfs]); iH = np.array([o[1] for o in orfs])
    np.save(NPY, np.vstack([iL, iH]))
    print(f"[scanned intergenic: {iL.size:,} ORFs, cached]")

def load(f):
    r = list(csv.DictReader(open(f)))
    return np.array([int(x['length']) for x in r]), np.array([int(x['HxH']) for x in r])
pL, pH = load(f"{INTER}/chelation_percluster_proteome.csv")
nL, nH = load(f"{INTER}/chelation_percluster_denovo.csv")

print("\n=== Q1: what IS the length-matched proteome for 22-aa intergenic ORFs? ===")
for cut in (20, 30, 40, 50, 70, 100):
    sub = pL <= cut
    if sub.sum():
        print(f"  proteome proteins <= {cut} aa: n={sub.sum():5d}  "
              f"His/100aa={100*(pH[sub].sum()/pL[sub].sum() if pL[sub].sum() else 0):.3f}  "
              f"H-x-H incidence={100*(pH[sub]>=1).mean():.2f}%")
print(f"  (proteome total n={pL.size:,}, median len={int(np.median(pL))} aa)")

print("\n=== Q2: LENGTH-STRATIFIED H-x-H incidence (length controlled within bins) ===")
bins = [(10,20),(20,30),(30,40),(40,60),(60,100),(100,10**9)]
print(f"  {'bin(aa)':>10} | {'interg n / inc%':>22} | {'proteome n / inc%':>22} | {'neORF n / inc%':>20} | Fisher(int vs prot)")
for lo,hi in bins:
    im = (iL>=lo)&(iL<hi); pm=(pL>=lo)&(pL<hi); nm=(nL>=lo)&(nL<hi)
    ii=100*(iH[im]>=1).mean() if im.sum() else float('nan')
    pi=100*(pH[pm]>=1).mean() if pm.sum() else float('nan')
    ni=100*(nH[nm]>=1).mean() if nm.sum() else float('nan')
    if im.sum()>=20 and pm.sum()>=20:
        iw=int((iH[im]>=1).sum()); pw=int((pH[pm]>=1).sum())
        OR,pf=stats.fisher_exact([[iw,im.sum()-iw],[pw,pm.sum()-pw]])
        fish=f"OR={OR:4.2f} p={pf:.1e}"
    else:
        fish="(insufficient n)"
    lab=f"[{lo},{hi if hi<10**9 else '∞'})"
    print(f"  {lab:>10} | {im.sum():8d} / {ii:6.2f} | {pm.sum():8d} / {pi:6.2f} | {nm.sum():7d} / {ni:6.2f} | {fish}")
print("\n  (intergenic median 22 aa; proteome median 469 aa; neORF median ~70 aa —")
print("   overlap where BOTH are well-represented is ~20-100 aa)")
