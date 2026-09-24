#!/usr/bin/env python3
"""Clean-intergenic ORFs (genome \\ (FlyBase-gene UNION RepeatMasker)) vs the
established proteome — H-x-H incidence & density, length-matched. Reuses the
paper's exact masking (denovo_hisspacing_analysis.load_intervals / masks) with a
freshly downloaded UCSC dm6 rmsk track; own 6-frame ORF scan records length +
H-x-H for EVERY ORF (>=10 aa), not just His>=2.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, csv
import numpy as np
CODE = os.path.join(_HERE, "..")
sys.path.insert(0, CODE)
import denovo_hisspacing_analysis as A     # main() is guarded -> import-safe
from hs_common import read_fasta
from scipy import stats

TMP = sys.argv[1]
A.SCRATCH = TMP                            # -> TMP/dmel.gtf.gz + TMP/dm6_rmsk.txt.gz
GENOME = f"{_HERE}/../../data/analysis/ref/dm6_genome.fa"
PROT = os.path.join(_HERE, "..", "..", "data", "intermediate", "chelation_percluster_proteome.csv")
CHRS = ["2L", "2R", "3L", "3R"]
B2I, LUT, H, STAR = A.BASE2INT, A.AA_LUT, A.ORD_H, A.ORD_STAR

seqs = {}
for hdr, s in read_fasta(GENOME):
    c = hdr.split()[0]
    if c in CHRS:
        seqs[c] = s.upper()

genes, reps = A.load_intervals()
print(f"masked FlyBase genes: {sum(len(v) for v in genes.values()):,} | "
      f"repeat intervals: {sum(len(v) for v in reps.values()):,}")


def scan_clean(seqb):
    """6-frame stop-to-stop ORFs >=10 aa in one clean region; return (L, hxh) list."""
    out = []
    x = B2I[seqb]
    xr = np.where(x == 4, 4, 3 - x).astype(np.uint8)[::-1]
    for src in (x, xr):
        n = src.size
        for fr in range(3):
            m = (n - fr) // 3
            if m < 10:
                continue
            cod = src[fr:fr + 3 * m].reshape(m, 3).astype(np.int64)
            valid = (cod < 4).all(axis=1)
            idx = np.where(valid, cod[:, 0] * 16 + cod[:, 1] * 4 + cod[:, 2], 0)
            aa = np.where(valid, LUT[idx], STAR).astype(np.uint8)
            stop = np.flatnonzero(aa == STAR)
            bounds = np.concatenate(([-1], stop, [aa.size]))
            for k in range(bounds.size - 1):
                s0, e0 = bounds[k] + 1, bounds[k + 1]
                L = e0 - s0
                if L < 10:
                    continue
                seg = aa[s0:e0]
                hm = seg == H
                hxh = int(np.count_nonzero(hm[:-2] & hm[2:])) if L > 2 else 0
                out.append((L, hxh))
    return out


orfs = []
for c in CHRS:
    seq = seqs[c]
    Lc = len(seq)
    seqb = np.frombuffer(seq.encode(), np.uint8)
    gmask = A.mask_from_intervals(genes[c], Lc)
    rmask = A.mask_from_intervals(reps[c], Lc)
    clean = (~gmask) & (~rmask)
    for s, e in A.runs_of_true(clean):
        if e - s >= 30:
            orfs.extend(scan_clean(seqb[s:e]))
    print(f"  chr{c}: clean_bp={int(clean.sum()):,}  cum ORFs>=10aa={len(orfs):,}")

L = np.array([o[0] for o in orfs]); hx = np.array([o[1] for o in orfs]); n = len(orfs)
inc = (hx >= 1).mean(); dens = 100 * hx.sum() / L.sum()
print(f"\n== CLEAN INTERGENIC ORFs ==")
print(f"  n={n:,}  median len={int(np.median(L))} aa (paper reports ~22)  "
      f"H-x-H incidence={100*inc:.2f}%  density={dens:.3f}/100aa")

# proteome
pr = list(csv.DictReader(open(PROT)))
pL = np.array([int(r['length']) for r in pr]); pH = np.array([int(r['HxH']) for r in pr])
print(f"== PROTEOME (all, n={len(pr):,}) ==  median len={int(np.median(pL))}  "
      f"incidence={100*(pH>=1).mean():.2f}%  density={100*pH.sum()/pL.sum():.3f}/100aa")

# length-matched proteome (nearest length, with replacement)
o = np.argsort(pL); pLs, pHs = pL[o], pH[o]
rng = np.random.default_rng(20260812)
mi = np.empty(n, int)
for i, x in enumerate(L):
    j = int(np.searchsorted(pLs, x))
    cand = [k for k in (j - 1, j, j + 1) if 0 <= k < len(pLs)]
    d = [abs(int(pLs[k]) - int(x)) for k in cand]
    best = [cand[t] for t in range(len(cand)) if d[t] == min(d)]
    mi[i] = rng.choice(best)
mL, mH = pLs[mi], pHs[mi]
m_inc = (mH >= 1).mean(); m_dens = 100 * mH.sum() / mL.sum()
print(f"== LENGTH-MATCHED PROTEOME (matched to intergenic) ==  median len={int(np.median(mL))}  "
      f"incidence={100*m_inc:.2f}%  density={m_dens:.3f}/100aa")

iw, mw = int((hx >= 1).sum()), int((mH >= 1).sum())
OR, pf = stats.fisher_exact([[iw, n - iw], [mw, n - mw]])
print(f"\n>>> CLEAN INTERGENIC vs LENGTH-MATCHED PROTEOME (H-x-H incidence):")
print(f"    {100*inc:.2f}% vs {100*m_inc:.2f}%  fold={inc/m_inc:.2f}x  OR={OR:.2f}  Fisher p={pf:.2e}")
print(f"    density {dens:.3f} vs {m_dens:.3f}/100aa  fold={dens/m_dens:.2f}x")
print(f"\n(reference: neORFs 9.7% incidence / 0.194 density; neORF-matched proteome 4.4% / 0.080)")
