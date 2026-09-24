#!/usr/bin/env python3
"""Cys% and His% (and Lys%) in neORFs (Grandchamp gc2023, n=7,867 -- the Fig 1 set), the Fig 1
length-matched proteome draw, the whole proteome (all dmel_pep entries >=10 aa, as in Fig 1), and known
metal-binders. Clean-intergenic values are taken from cxc_analysis.log (same machinery)."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, glob
import numpy as np
sys.path.insert(0, os.path.join(_HERE, ".."))
from hs_common import load_denovo, load_proteins
D = _HERE
S = f"{D}/data/md_hp6_scratch"
dn = [d["prot"] for d in load_denovo(sorted(glob.glob(f"{S}/gc2023/*_DeNovoGenesNuc.fa")), min_aa=10)]
prot = load_proteins(f"{S}/data/dmel_pep.fa.gz", min_aa=10)
meta = load_proteins(f"{S}/data/dmel_metalbinding.fasta", min_aa=10)
dl = np.array([len(s) for s in dn]); plen = np.array([len(s) for s in prot]); order = np.argsort(plen, kind="stable")
ps = [prot[i] for i in order]; pls = plen[order]; rng = np.random.default_rng(20260812); matched = []
for L in dl:   # Fig 1 nearest-length whole-protein draw (denovo_fig1_wholeprotein.py)
    j = int(np.searchsorted(pls, int(L))); cand = [k for k in (j - 1, j, j + 1) if 0 <= k < len(ps)]
    d = np.array([abs(int(pls[k]) - int(L)) for k in cand])
    matched.append(ps[int(rng.choice([cand[t] for t in range(len(cand)) if d[t] == d.min()]))])
def f(seqs, r): return 100 * sum(s.count(r) for s in seqs) / sum(map(len, seqs))
rows = [("neORFs", dn), ("length-matched proteome (Fig 1 draw)", matched), ("whole proteome", prot), ("known metal-binders", meta)]
lm = {r: f(matched, r) for r in "CHK"}; wp = {r: f(prot, r) for r in "CHK"}
for r, name in [("C", "Cys"), ("H", "His"), ("K", "Lys")]:
    print(f"\n{name}:  {'set':38s} {'n':>6s} {'%':>6s} {'vs whole':>9s} {'vs length-matched':>18s}")
    for k, seqs in rows:
        v = f(seqs, r); print(f"      {k:38s} {len(seqs):6d} {v:6.2f} {v/wp[r]:8.2f}x {v/lm[r]:17.2f}x")
print(f"\nlength-matched draw: {len(matched)} draws from {len(set(matched))} distinct proteins")
