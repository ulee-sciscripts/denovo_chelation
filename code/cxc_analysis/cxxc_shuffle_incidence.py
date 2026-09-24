#!/usr/bin/env python3
"""C-x-x-C incidence observed vs within-sequence shuffle (length and Cys count preserved).
P(>=1 CxxC | L, nC) under a random permutation depends only on (L, nC): Monte Carlo (2,000 draws)
for L<=300, Poisson approximation 1-exp(-nC(nC-1)(L-3)/(L(L-1))) above. Expected incidence = mean
over ORFs. Sets: clean intergenic ORFs, neORFs, distinct proteome, known metal-binders."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys
import numpy as np
sys.path.insert(0, os.path.join(_HERE, ".."))
import denovo_hisspacing_analysis as A
from hs_common import read_fasta, translate, load_proteins
A.SCRATCH = sys.argv[1]
S = f"{_HERE}/data/md_hp6_scratch/data"
rng = np.random.default_rng(20260923)
def obs(p):
    a = np.frombuffer(p.encode(), np.uint8) == 67
    return len(p), int(a.sum()), int(np.any(a[:-3] & a[3:])) if len(p) > 3 else 0
cache = {}
def pexp(L, n):
    if n < 2 or L < 4: return 0.0
    if (L, n) in cache: return cache[(L, n)]
    if L > 300:
        v = 1 - np.exp(-n * (n - 1) * (L - 3) / (L * (L - 1)))
    else:
        pos = np.argsort(rng.random((2000, L)), axis=1)[:, :n]
        m = np.zeros((2000, L), bool); np.put_along_axis(m, pos, True, axis=1)
        v = float(np.any(m[:, :-3] & m[:, 3:], axis=1).mean())
    cache[(L, n)] = v; return v
seen = set(); neo = []
for h, s in read_fasta(f"{_HERE}/../../data/analysis/neORF.fa"):
    p, _ = translate(s.upper(), truncate_at_stop=True)
    if len(p) >= 10 and p not in seen: seen.add(p); neo.append(obs(p))
prot = [obs(p) for p in dict.fromkeys(load_proteins(f"{S}/dmel_pep.fa.gz", min_aa=10))]
meta = [obs(p) for p in dict.fromkeys(load_proteins(f"{S}/dmel_metalbinding.fasta", min_aa=10))]
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
                isC = aa == C; mc = np.zeros(aa.size, np.int64); mc[:-3] = isC[:-3] & isC[3:]
                cs = np.concatenate(([0], np.cumsum(mc))); cc = np.concatenate(([0], np.cumsum(isC)))
                bnd = np.concatenate(([-1], np.flatnonzero(aa == STAR), [aa.size]))
                a0 = bnd[:-1] + 1; b0 = bnd[1:]; ok = (b0 - a0) >= 10; a0, b0 = a0[ok], b0[ok]
                out.append(np.stack([b0 - a0, cc[b0] - cc[a0], (cs[np.maximum(b0 - 3, a0)] - cs[a0]) > 0], 1))
ig = np.concatenate(out)
print(f"{'set':28s} {'n':>9s} {'observed':>9s} {'shuffled':>9s} {'obs/shuf':>9s}")
for k, v in [("clean intergenic ORFs", ig), ("neORFs", np.array(neo)),
             ("whole proteome (distinct)", np.array(prot)), ("known metal-binders", np.array(meta))]:
    L, n, o = v[:, 0], v[:, 1], v[:, 2]
    e = np.array([pexp(int(a), int(b)) for a, b in zip(L, n)])
    print(f"{k:28s} {len(v):9d} {100*o.mean():8.2f}% {100*e.mean():8.2f}% {o.mean()/e.mean():9.2f}")
