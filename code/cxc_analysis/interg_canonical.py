#!/usr/bin/env python3
"""Reconcile with fig1's canonical proteome (dmel_pep.fa.gz), then compute the
intergenic-vs-proteome enrichment against the SAME baseline fig1 committed (4.39%).
All three sets matched to the neORF length distribution with fig1's exact method.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, csv
import numpy as np
from scipy import stats
CODE = os.path.join(_HERE, "..")
sys.path.insert(0, CODE)
from hs_common import load_proteins
SP = f"{_HERE}/data"
DMEL = f"{_HERE}/data/data/dmel_pep.fa.gz"
INTER = os.path.join(_HERE, "..", "..", "data", "intermediate")
SEED = 20260812

def hxh1(s):
    a = np.frombuffer(s.encode(), np.uint8); m = a == 72       # 'H'
    return 1 if (len(s) > 2 and bool(np.any(m[:-2] & m[2:]))) else 0

prot = load_proteins(DMEL, min_aa=10)
pL = np.array([len(s) for s in prot]); pH = np.array([hxh1(s) for s in prot])
print(f"canonical proteome (dmel_pep.fa.gz, min_aa=10): n={len(prot):,}  "
      f"median len={int(np.median(pL))}  overall HxH incidence={100*pH.mean():.2f}%")
print(f"  (percluster_proteome.csv had n=30,801 — different set/filtering)")
for cut in (30, 50, 70, 100):
    s = pL <= cut
    print(f"  proteome <= {cut} aa: n={int(s.sum()):5d}  incidence={100*pH[s].mean() if s.sum() else 0:.2f}%")

def load_csv(f):
    r = list(csv.DictReader(open(f)))
    return (np.array([int(x['length']) for x in r]),
            np.array([1 if int(x['HxH']) >= 1 else 0 for x in r]))
nL, nH = load_csv(f"{INTER}/chelation_percluster_denovo.csv")
LH = np.load(f"{SP}/interg_LH.npy"); iL = LH[0]; iH = (LH[1] >= 1).astype(int)

def match_to(tgtL, srcL, srcH, seed):
    o = np.argsort(srcL, kind="stable"); sL, sH = srcL[o], srcH[o]
    rng = np.random.default_rng(seed)
    out = np.empty(len(tgtL), int)
    for i, x in enumerate(tgtL):
        j = int(np.searchsorted(sL, int(x)))
        cand = [k for k in (j - 1, j, j + 1) if 0 <= k < len(sL)]
        d = [abs(int(sL[k]) - int(x)) for k in cand]
        best = [cand[t] for t in range(len(cand)) if d[t] == min(d)]
        out[i] = sH[rng.choice(best)]
    return out

mp = match_to(nL, pL, pH, SEED)          # canonical proteome @ neORF lengths (should ~4.39%)
mi = match_to(nL, iL, iH, SEED)          # intergenic @ neORF lengths
inc = lambda h: 100 * h.mean()
print(f"\n=== all @ neORF length distribution (N={len(nL)}) ===")
print(f"  neORFs                         : {inc(nH):.2f}%")
print(f"  intergenic ORFs @ neORF lengths: {inc(mi):.2f}%")
print(f"  canonical proteome (baseline)  : {inc(mp):.2f}%   <-- validate vs fig1's 4.39%")

def fish(a, b):
    aw, bw = int(a.sum()), int(b.sum())
    OR, p = stats.fisher_exact([[aw, len(a) - aw], [bw, len(b) - bw]])
    return OR, p
OR, p = fish(mi, mp); ORn, pn = fish(nH, mp)
print(f"\n  intergenic vs proteome : fold={inc(mi)/inc(mp):.2f}x  OR={OR:.2f}  p={p:.2e}")
print(f"  neORF      vs proteome : fold={inc(nH)/inc(mp):.2f}x  OR={ORn:.2f}  p={pn:.2e}   <-- should ~= committed 2.2x")
# stability over seeds
fi = [inc(match_to(nL, iL, iH, 900 + s)) / inc(match_to(nL, pL, pH, 900 + s)) for s in range(8)]
print(f"  intergenic/proteome fold across 8 seeds: median {np.median(fi):.2f}x  range {min(fi):.2f}-{max(fi):.2f}")
