#!/usr/bin/env python3
"""Order-fold (observed / within-sequence scramble expectation, analytic; paper's method) for C-x-C,
C-x-x-C, H-x-H in annotated FlyBase protein-coding genes by NCBI gene-age stratum (gene_age_ncbi.json),
longest isoform per gene; neORFs as reference. Bootstrap 95% CI over genes."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, gzip, re, json
import numpy as np
sys.path.insert(0, os.path.join(_HERE, ".."))
from hs_common import read_fasta, translate
D = _HERE
PEP = f"{_HERE}/data/md_hp6_scratch/data/dmel_pep.fa.gz"
AGE = json.load(open(f"{D}/data/gene_age_ncbi.json"))
best = {}
for h, s in read_fasta(PEP):
    s = s.rstrip("*")
    if len(s) < 10: continue
    g = re.search(r"gene:(FBgn\d+)", h).group(1)
    if g not in best or len(s) > len(best[g]): best[g] = s
seen = set(); neo = []
for h, s in read_fasta(f"{_HERE}/../../data/analysis/neORF.fa"):
    p, _ = translate(s.upper(), truncate_at_stop=True)
    if len(p) >= 10 and p not in seen: seen.add(p); neo.append(p)
def oe(s, r, d):
    L = len(s)
    if L <= d: return 0.0, 0.0
    m = np.frombuffer(s.encode(), np.uint8) == ord(r); n = m.sum()
    return float(np.count_nonzero(m[:L-d] & m[d:])), float(n * (n - 1) * (L - d) / (L * (L - 1)))
M = [("C-x-C", "C", 2), ("C-x-x-C", "C", 3), ("H-x-H", "H", 2)]
rng = np.random.default_rng(20260923)
def fold(seqs):
    cells = []
    for _, r, d in M:
        a = np.array([oe(s, r, d) for s in seqs]); o, e = a[:, 0], a[:, 1]; n = len(seqs)
        bs = [o[ix].sum() / e[ix].sum() for ix in (rng.integers(0, n, n) for _ in range(1000))]
        cells.append(f"{o.sum()/e.sum():5.2f} [{np.percentile(bs,2.5):4.2f}-{np.percentile(bs,97.5):4.2f}]")
    return cells
STR = ["D.mel only / no ortholog", "melanogaster subgroup/group", "Sophophora/Drosophila",
       "Diptera", "Insecta", "Arthropoda or older"]
groups = [(f"{i+1} {s}", [best[g] for g in best if AGE.get(g) == i]) for i, s in enumerate(STR)]
groups.insert(0, ("YOUNG (strata 1-3, Drosophila-restricted)", [best[g] for g in best if AGE.get(g, 9) <= 2]))
groups.insert(1, ("OLD (strata 5-6, Insecta or older)", [best[g] for g in best if AGE.get(g, -1) >= 4]))
groups.append(("neORFs (reference)", neo))
print(f"{'group':44s} {'n':>6s} {'med aa':>6s} " + " ".join(f"{m[0]:>18s}" for m in M) + "   Cys%  His%")
for k, seqs in groups:
    tot = sum(map(len, seqs))
    print(f"{k:44s} {len(seqs):6d} {int(np.median([len(s) for s in seqs])):6d} " + " ".join(f"{c:>18s}" for c in fold(seqs))
          + f"  {100*sum(s.count('C') for s in seqs)/tot:5.2f} {100*sum(s.count('H') for s in seqs)/tot:5.2f}")

# all annotated genes (longest isoform per gene), same method -- the like-for-like "whole proteome" for age comparisons
print(f"{'ALL genes (longest isoform)':44s} {len(best):6d} {int(np.median([len(s) for s in best.values()])):6d} "
      + " ".join(f"{c:>18s}" for c in fold(list(best.values()))))
