#!/usr/bin/env python3
"""Put clean-intergenic ORFs on the neORF length distribution, then compare to
the SAME length-matched proteome baseline the neORFs use. Three sets at one
(neORF) length distribution: neORF, length-matched intergenic, length-matched
proteome. Answers: do intergenic ORFs (at neORF lengths) show the neORF-like
H-x-H enrichment vs known proteins, using the well-populated ~64-aa baseline
(not the degenerate 22-aa one)?
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import csv
import numpy as np
from scipy import stats

SP = f"{_HERE}/data"
INTER = os.path.join(_HERE, "..", "..", "data", "intermediate")

LH = np.load(f"{SP}/interg_LH.npy"); iL, iH = LH[0], LH[1]        # 1.47M clean-intergenic ORFs
def load(f):
    r = list(csv.DictReader(open(f)))
    return np.array([int(x['length']) for x in r]), np.array([int(x['HxH']) for x in r])
nL, nH = load(f"{INTER}/chelation_percluster_denovo.csv")          # neORFs
pL, pH = load(f"{INTER}/chelation_percluster_proteome.csv")        # proteome

def match_to(tgtL, srcL, srcH, seed):
    """For each target length, sample (with replacement) a source of nearest length."""
    o = np.argsort(srcL); sL, sH = srcL[o], srcH[o]
    rng = np.random.default_rng(seed)
    out = np.empty(len(tgtL), int)
    for i, x in enumerate(tgtL):
        j = int(np.searchsorted(sL, x))
        cand = [k for k in (j - 1, j, j + 1) if 0 <= k < len(sL)]
        d = [abs(int(sL[k]) - int(x)) for k in cand]
        best = [cand[t] for t in range(len(cand)) if d[t] == min(d)]
        out[i] = sH[rng.choice(best)]
    return out

# match intergenic and proteome to the neORF length distribution
miH = match_to(nL, iL, iH, 20260812)          # intergenic @ neORF lengths
mpH = match_to(nL, pL, pH, 20260813)          # proteome  @ neORF lengths (reproduces fig1 ~4.4%)

def inc(h): return 100 * (h >= 1).mean()
N = len(nL)
print(f"All three sets on the neORF length distribution (N={N}, median len={int(np.median(nL))} aa):\n")
print(f"  neORFs (realized)                : H-x-H incidence = {inc(nH):5.2f}%")
print(f"  intergenic ORFs @ neORF lengths  : H-x-H incidence = {inc(miH):5.2f}%")
print(f"  proteome @ neORF lengths (baseline): H-x-H incidence = {inc(mpH):5.2f}%   [fig1 committed = 4.39%]")

def fisher(a, b):
    aw, bw = int((a >= 1).sum()), int((b >= 1).sum())
    OR, p = stats.fisher_exact([[aw, len(a) - aw], [bw, len(b) - bw]])
    return OR, p, aw, bw

OR, p, iw, pw = fisher(miH, mpH)
print(f"\n  >>> length-matched INTERGENIC vs length-matched PROTEOME:")
print(f"      {inc(miH):.2f}% vs {inc(mpH):.2f}%  fold={inc(miH)/inc(mpH):.2f}x  OR={OR:.2f}  Fisher p={p:.2e}")
OR2, p2, _, _ = fisher(miH, nH)
print(f"  >>> length-matched INTERGENIC vs neORFs:")
print(f"      {inc(miH):.2f}% vs {inc(nH):.2f}%  fold={inc(miH)/inc(nH):.2f}x  OR={OR2:.2f}  Fisher p={p2:.2e}")
OR3, p3, _, _ = fisher(nH, mpH)
print(f"  >>> neORFs vs length-matched PROTEOME (the manuscript's 2.2x, re-derived):")
print(f"      {inc(nH):.2f}% vs {inc(mpH):.2f}%  fold={inc(nH)/inc(mpH):.2f}x  OR={OR3:.2f}  Fisher p={p3:.2e}")

# stability: repeat the intergenic match over several seeds
folds = []
for s in range(10):
    m = match_to(nL, iL, iH, 900 + s)
    folds.append(inc(m) / inc(mpH))
folds = np.array(folds)
print(f"\n  [intergenic-vs-proteome fold across 10 match seeds: "
      f"median {np.median(folds):.2f}x, range {folds.min():.2f}-{folds.max():.2f}]")
