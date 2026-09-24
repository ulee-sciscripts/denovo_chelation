#!/usr/bin/env python3
"""Val at the central (x) position of C-x-C in neORFs (the (TG)n reading: TGT.GTG.TGT = Cys-Val-Cys),
parallel to thr_middle_binomial.py. Binomial exact vs background Val frequency."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, os.path.join(_HERE, ".."))
from hs_common import read_fasta, translate
from scipy import stats
seen, prots = set(), []
for h, s in read_fasta(f"{_HERE}/../../data/analysis/neORF.fa"):
    p, _ = translate(s.upper(), truncate_at_stop=True)
    if len(p) >= 10 and p not in seen:
        seen.add(p); prots.append(p)
tot = sum(len(p) for p in prots); f_V = sum(p.count("V") for p in prots) / tot
mid = [p[i+1] for p in prots for i in range(len(p)-2) if p[i] == "C" and p[i+2] == "C"]
k, n = mid.count("V"), len(mid)
bt = stats.binomtest(k, n, f_V); ci = bt.proportion_ci(0.95)
print(f"n neORF={len(prots)}  background Val = {100*f_V:.2f}%")
print(f"C-x-C middle = Val: {k}/{n} = {100*k/n:.1f}% [{100*ci.low:.1f}-{100*ci.high:.1f}]  -> {k/n/f_V:.2f}x background, binomial p = {bt.pvalue:.2e}")
from collections import Counter
print("top middle residues:", Counter(mid).most_common(6))
