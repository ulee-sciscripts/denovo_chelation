#!/usr/bin/env python3
"""Statistics for the two summary figures (fig_motif_enrichment, fig_prevalence).
Categories: clean intergenic ORFs; neORFs (canonical neORF.fa, n=7,849); young de novo (Heames et al. 2020
D. mel de novo genes); young non-de novo (NCBI Drosophila-restricted minus Heames de novo); all young (union);
older (NCBI Insecta or older); known metal-binders (dmel_metalbinding.fasta, distinct). Annotated genes: longest current FlyBase isoform per FBgn.
(1) Order-fold (obs / analytic within-sequence scramble) for H-x-H, C-x-C, C-x-x-C; 95% bootstrap CI over
    sequences; two-sided p vs 1.  (2) His% and Cys% (pooled); 95% bootstrap CI; two-sided p vs intergenic
    (independent bootstraps).  Exact tests (used for the figure stars): (1) Fisher 2x2 of position pairs at
    spacing d, motif vs not, observed vs scramble-expected (rounded); (2) binomial exact of residue count vs the
    intergenic residue frequency (residues treated as independent).  B = 2,000 (intergenic: Poisson bootstrap weights).
Usage: fig_stats.py <interval dir>   (data/interg_intervals)"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, re, csv, json
import numpy as np
from scipy.stats import fisher_exact, binomtest
sys.path.insert(0, os.path.join(_HERE, ".."))
import denovo_hisspacing_analysis as A
from hs_common import read_fasta, translate, load_proteins
A.SCRATCH = sys.argv[1]
D = _HERE
B = 2000; rng = np.random.default_rng(20260924)
MOT = [("HxH", "H", 2), ("CxC", "C", 2), ("CxxC", "C", 3)]

def feats(seqs):
    """per-sequence array: L, nH, nC, obs/exp for each motif."""
    out = np.zeros((len(seqs), 9))
    for i, s in enumerate(seqs):
        a = np.frombuffer(s.encode(), np.uint8); L = len(s); out[i, 0] = L
        mh, mc = a == 72, a == 67; out[i, 1], out[i, 2] = mh.sum(), mc.sum()
        for j, (_, r, d) in enumerate(MOT):
            m = mh if r == "H" else mc; n = m.sum()
            if L > d:
                out[i, 3 + 2*j] = np.count_nonzero(m[:L-d] & m[d:]); out[i, 4 + 2*j] = n*(n-1)*(L-d)/(L*(L-1))
    return out

# ---- annotated gene sets ----
best = {}
for h, s in read_fasta(f"{D}/data/md_hp6_scratch/data/dmel_pep.fa.gz"):
    s = s.rstrip("*")
    if len(s) < 10: continue
    g = re.search(r"gene:(FBgn\d+)", h).group(1)
    if g not in best or len(s) > len(best[g]): best[g] = s
AGE = json.load(open(f"{D}/data/gene_age_ncbi.json"))
dn = {r["gene"] for r in csv.DictReader(open(f"{D}/data/heames2020/heames2020_dmel_all_orphans.csv"))
      if r["mechanism"] in ("denovo", "denovo-intron")}
young_ndn = {g for g in best if AGE.get(g, 9) <= 2 and g not in dn}
young_dn = {g for g in dn if g in best}
old = {g for g in best if AGE.get(g, -1) >= 4}
seen = set(); neo = []
for h, s in read_fasta(f"{_HERE}/../../data/analysis/neORF.fa"):
    p, _ = translate(s.upper(), truncate_at_stop=True)
    if len(p) >= 10 and p not in seen: seen.add(p); neo.append(p)

# ---- intergenic scan: per-ORF features, vectorised ----
import os
IGF = f"{D}/data/interg_features.npy"
if os.path.exists(IGF):
    IG = np.load(IGF)
else:
    B2I, LUT, STAR = A.BASE2INT, A.AA_LUT, A.ORD_STAR
    seqs = {}
    for h, s in read_fasta(f"{_HERE}/../../data/analysis/ref/dm6_genome.fa"):
        c = h.split()[0]
        if c in ("2L", "2R", "3L", "3R"): seqs[c] = s.upper()
    genes, reps = A.load_intervals(); chunks = []
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
                    cod = src[fr:fr + 3*m].reshape(m, 3).astype(np.int64)
                    valid = (cod < 4).all(axis=1)
                    aa = np.where(valid, LUT[np.where(valid, cod[:, 0]*16 + cod[:, 1]*4 + cod[:, 2], 0)], STAR).astype(np.uint8)
                    bnd = np.concatenate(([-1], np.flatnonzero(aa == STAR), [aa.size]))
                    a0 = bnd[:-1] + 1; b0 = bnd[1:]; ok = (b0 - a0) >= 10; a0, b0 = a0[ok], b0[ok]
                    if not a0.size: continue
                    isH, isC = aa == 72, aa == 67
                    cH = np.concatenate(([0], np.cumsum(isH))); cC = np.concatenate(([0], np.cumsum(isC)))
                    L = (b0 - a0).astype(float); nH = (cH[b0] - cH[a0]).astype(float); nC = (cC[b0] - cC[a0]).astype(float)
                    cols = [L, nH, nC]
                    for _, r, d in MOT:
                        iso = isH if r == "H" else isC
                        pm = np.zeros(aa.size, np.int64); pm[:-d] = iso[:-d] & iso[d:]
                        cp = np.concatenate(([0], np.cumsum(pm)))
                        n = nH if r == "H" else nC
                        cols += [(cp[np.maximum(b0 - d, a0)] - cp[a0]).astype(float), n*(n-1)*(L-d)/(L*(L-1))]
                    chunks.append(np.stack(cols, 1))
    IG = np.concatenate(chunks); np.save(IGF, IG)


CATS = [("intergenic", "intergenic ORFs", IG), ("neORF", "neORFs", feats(neo)),
        ("young_dn", "young de novo", feats([best[g] for g in young_dn])),
        ("young_ndn", "young non-de novo", feats([best[g] for g in young_ndn])),
        ("young_all", "all young", feats([best[g] for g in young_dn | young_ndn])),
        ("older", "older", feats([best[g] for g in old])),
        ("metal", "metal-binders", feats(list(dict.fromkeys(load_proteins(f"{D}/data/md_hp6_scratch/data/dmel_metalbinding.fasta", min_aa=10)))))]

def boots(F):
    """bootstrap matrix of column sums (B x 9); Poisson weights for very large sets."""
    n = len(F); out = np.empty((B, F.shape[1]))
    for b in range(B):
        if n > 100000: out[b] = rng.poisson(1.0, n) @ F
        else: out[b] = F[rng.integers(0, n, n)].sum(0)
    return out
def star(p): return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
def p2(x):   # two-sided bootstrap p that x differs from 0, floored at 1/(B+1)
    return max(min(1.0, 2*min((x <= 0).mean(), (x >= 0).mean())), 1/(B + 1))

mrows, prow = [], []; bs = {}
for key, lab, F in CATS:
    bs[key] = boots(F); tot = F.sum(0)
    for j, (m, _, _) in enumerate(MOT):
        f = tot[3+2*j] / tot[4+2*j]; bf = bs[key][:, 3+2*j] / bs[key][:, 4+2*j]
        p = p2(bf - 1)
        d = MOT[j][2]; npos = int(np.clip(F[:, 0] - d, 0, None).sum())      # position pairs at spacing d
        o, e = int(tot[3+2*j]), int(round(tot[4+2*j]))
        OR, pf = fisher_exact([[o, npos - o], [e, npos - e]])               # observed vs scramble-expected
        mrows.append([key, lab, len(F), m, o, round(tot[4+2*j], 1), f,
                      np.percentile(bf, 2.5), np.percentile(bf, 97.5), p, star(p), npos, OR, pf, star(pf)])
for key, lab, F in CATS:
    tot = F.sum(0)
    for res, col in [("His", 1), ("Cys", 2)]:
        v = 100*tot[col]/tot[0]; bv = 100*bs[key][:, col]/bs[key][:, 0]
        bi = 100*bs["intergenic"][:, col]/bs["intergenic"][:, 0]
        p = np.nan if key == "intergenic" else p2(bv - bi)
        k, N = int(tot[col]), int(tot[0]); p0 = CATS[0][2][:, col].sum() / CATS[0][2][:, 0].sum()
        pb = np.nan if key == "intergenic" else binomtest(k, N, p0).pvalue      # residues vs intergenic freq
        prow.append([key, lab, len(F), res, v, np.percentile(bv, 2.5), np.percentile(bv, 97.5), p,
                     "ref" if key == "intergenic" else star(p), k, N, pb, "ref" if key == "intergenic" else star(pb)])
with open(f"{D}/results/fig_motif_enrichment.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["key", "label", "n", "motif", "obs", "exp", "fold", "lo", "hi", "p_vs1", "star", "n_pospairs", "fisher_OR", "fisher_p", "fisher_star"]); w.writerows(mrows)
with open(f"{D}/results/fig_prevalence.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["key", "label", "n", "residue", "pct", "lo", "hi", "p_vs_intergenic", "star", "count", "n_residues", "binom_p", "binom_star"]); w.writerows(prow)
for r in mrows: print("  ".join(str(x) if not isinstance(x, float) else f"{x:.4g}" for x in r))
for r in prow: print("  ".join(str(x) if not isinstance(x, float) else f"{x:.4g}" for x in r))
