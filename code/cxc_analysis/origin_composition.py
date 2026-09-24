#!/usr/bin/env python3
"""Origin-composition control: Lys%/Cys%/His% in neORFs (zinc_analysis/neORF.fa, n=7,849) and clean
(repeat-masked, non-genic) intergenic ORFs (>=10 aa, 6-frame stop-to-stop, arms 2L/2R/3L/3R), to compare
with young vs old cytosolic short FlyBase proteins (cys_short_age.py). Intergenic scan reuses the paper's
denovo_hisspacing_analysis machinery (same definition as gbe_revision/intergenic_analysis)."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys
import numpy as np
sys.path.insert(0, os.path.join(_HERE, ".."))
import denovo_hisspacing_analysis as A
from hs_common import read_fasta, translate
A.SCRATCH = sys.argv[1]
AAS = "KCHRDE"
def comp(counts, tot): return "  ".join(f"{a} {100*counts[a]/tot:5.2f}%" for a in AAS)

seen = set(); cnt = {a: 0 for a in AAS}; tot = 0
for h, s in read_fasta(f"{_HERE}/../../data/analysis/neORF.fa"):
    p, _ = translate(s.upper(), truncate_at_stop=True)
    if len(p) < 10 or p in seen: continue
    seen.add(p); tot += len(p)
    for a in AAS: cnt[a] += p.count(a)
print(f"neORFs (n={len(seen)}):             {comp(cnt, tot)}")

B2I, LUT, STAR = A.BASE2INT, A.AA_LUT, A.ORD_STAR
seqs = {}
for h, s in read_fasta(f"{_HERE}/../../data/analysis/ref/dm6_genome.fa"):
    c = h.split()[0]
    if c in ("2L", "2R", "3L", "3R"): seqs[c] = s.upper()
genes, reps = A.load_intervals()
hist = np.zeros(256, np.int64); ntot = 0; norf = 0
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
                bnd = np.concatenate(([-1], np.flatnonzero(aa == STAR), [aa.size]))
                for k in range(bnd.size - 1):
                    a0, b0 = bnd[k] + 1, bnd[k + 1]
                    if b0 - a0 < 10: continue
                    hist += np.bincount(aa[a0:b0], minlength=256); ntot += b0 - a0; norf += 1
cnt = {a: int(hist[ord(a)]) for a in AAS}
print(f"clean intergenic ORFs (n={norf}): {comp(cnt, ntot)}")
