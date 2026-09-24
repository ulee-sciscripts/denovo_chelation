#!/usr/bin/env python3
"""C-x-x-C (C at i, i+3) prevalence: incidence (>=1 motif), density (/100 aa), length-stratified
incidence, in clean intergenic ORFs (6-frame, repeat-masked, non-genic, 2L/2R/3L/3R), neORFs
(zinc_analysis/neORF.fa, n=7,849), distinct proteome (dmel_pep.fa.gz), and metal-binders."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys
import numpy as np
sys.path.insert(0, os.path.join(_HERE, ".."))
import denovo_hisspacing_analysis as A
from hs_common import read_fasta, translate, load_proteins
A.SCRATCH = sys.argv[1]
S = f"{_HERE}/data/md_hp6_scratch/data"
def cxxc(p):
    a = np.frombuffer(p.encode(), np.uint8) == 67
    return int(np.count_nonzero(a[:-3] & a[3:])) if len(p) > 3 else 0
def stats_of(LC):
    L, c = LC[:, 0], LC[:, 1]
    return len(L), int(np.median(L)), 100 * (c >= 1).mean(), 100 * c.sum() / L.sum()
seen = set(); neo = []
for h, s in read_fasta(f"{_HERE}/../../data/analysis/neORF.fa"):
    p, _ = translate(s.upper(), truncate_at_stop=True)
    if len(p) >= 10 and p not in seen: seen.add(p); neo.append((len(p), cxxc(p)))
prot = list(dict.fromkeys(load_proteins(f"{S}/dmel_pep.fa.gz", min_aa=10)))
meta = list(dict.fromkeys(load_proteins(f"{S}/dmel_metalbinding.fasta", min_aa=10)))
B2I, LUT, STAR, C = A.BASE2INT, A.AA_LUT, A.ORD_STAR, ord("C")
seqs = {}
for h, s in read_fasta(f"{_HERE}/../../data/analysis/ref/dm6_genome.fa"):
    c = h.split()[0]
    if c in ("2L", "2R", "3L", "3R"): seqs[c] = s.upper()
genes, reps = A.load_intervals(); out = []
for c, sq in seqs.items():
    seqb = np.frombuffer(sq.encode(), np.uint8); Lc = len(sq)
    clean = (~A.mask_from_intervals(genes[c], Lc)) & (~A.mask_from_intervals(reps[c], Lc))
    for s0, e0 in A.runs_of_true(clean):
        if e0 - s0 < 30: continue
        x = B2I[seqb[s0:e0]]; xr = np.where(x == 4, 4, 3 - x).astype(np.uint8)[::-1]
        for src in (x, xr):
            for fr in range(3):
                m = (src.size - fr) // 3
                if m < 10: continue
                cod = src[fr:fr + 3 * m].reshape(m, 3).astype(np.int64)
                valid = (cod < 4).all(axis=1)
                aa = np.where(valid, LUT[np.where(valid, cod[:, 0]*16 + cod[:, 1]*4 + cod[:, 2], 0)], STAR).astype(np.uint8)
                mc = np.zeros(aa.size, np.int64); isC = aa == C
                mc[:-3] = isC[:-3] & isC[3:]            # motif start at i (i..i+3 must lie in one ORF: checked below)
                cs = np.concatenate(([0], np.cumsum(mc)))
                bnd = np.concatenate(([-1], np.flatnonzero(aa == STAR), [aa.size]))
                a0 = bnd[:-1] + 1; b0 = bnd[1:]; ok = (b0 - a0) >= 10
                a0, b0 = a0[ok], b0[ok]
                out.append(np.stack([b0 - a0, cs[np.maximum(b0 - 3, a0)] - cs[a0]], 1))
ig = np.concatenate(out)
sets = {"clean intergenic ORFs": ig, "neORFs": np.array(neo),
        "whole proteome (distinct)": np.array([(len(p), cxxc(p)) for p in prot]),
        "known metal-binders": np.array([(len(p), cxxc(p)) for p in meta])}
print(f"{'set':28s} {'n':>9s} {'median aa':>9s} {'incidence':>10s} {'per 100 aa':>11s}")
for k, v in sets.items():
    n, md, inc, den = stats_of(v); print(f"{k:28s} {n:9d} {md:9d} {inc:9.2f}% {den:11.3f}")
print("\nlength-stratified incidence (% with >=1 C-x-x-C; n in parentheses)")
bins = [(10, 30), (30, 60), (60, 100), (100, 150), (150, 300), (300, 10**6)]
print(f"{'set':28s} " + " ".join(f"{f'[{lo},{hi if hi<10**6 else chr(8734)})':>18s}" for lo, hi in bins))
for k, v in sets.items():
    cells = []
    for lo, hi in bins:
        s = v[(v[:, 0] >= lo) & (v[:, 0] < hi)]
        cells.append(f"{100*(s[:,1]>=1).mean():6.2f}% ({len(s):>7d})" if len(s) >= 20 else f"{'-':>18s}")
    print(f"{k:28s} " + " ".join(f"{c:>18s}" for c in cells))
