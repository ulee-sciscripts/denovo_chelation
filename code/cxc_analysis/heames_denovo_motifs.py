#!/usr/bin/env python3
"""Motif order-fold (C-x-C, C-x-x-C, H-x-H; within-sequence scramble, analytic) in de novo genes
(Heames et al. 2020, D. melanogaster: mechanism 'denovo' + 'denovo-intron') vs other orphan classes
('divergent' = orphan by divergence from a homolog; 'putative'), NCBI-young non-de-novo genes, old genes,
and neORFs. Longest current FlyBase isoform per FBgn (dmel_pep); bootstrap 95% CI over genes."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, re, csv, json
import numpy as np
sys.path.insert(0, os.path.join(_HERE, ".."))
from hs_common import read_fasta, translate
D = _HERE
PEP = f"{D}/data/md_hp6_scratch/data/dmel_pep.fa.gz"
AGE = json.load(open(f"{D}/data/gene_age_ncbi.json"))
best = {}
for h, s in read_fasta(PEP):
    s = s.rstrip("*")
    if len(s) < 10: continue
    g = re.search(r"gene:(FBgn\d+)", h).group(1)
    if g not in best or len(s) > len(best[g]): best[g] = s
orph = list(csv.DictReader(open(f"{D}/data/heames2020/heames2020_dmel_all_orphans.csv")))
mech = {}
for r in orph:
    m = "de novo" if r["mechanism"] in ("denovo", "denovo-intron") else r["mechanism"]
    mech.setdefault(r["gene"], set()).add(m)
def genes(m): return {g for g, ms in mech.items() if m in ms}
dn_genes = genes("de novo")
seen = set(); neo = []
for h, s in read_fasta(f"{_HERE}/../../data/analysis/neORF.fa"):
    p, _ = translate(s.upper(), truncate_at_stop=True)
    if len(p) >= 10 and p not in seen: seen.add(p); neo.append(p)
def oe(s, r, d):
    L = len(s); m = np.frombuffer(s.encode(), np.uint8) == ord(r); n = m.sum()
    return (np.count_nonzero(m[:L-d] & m[d:]), n*(n-1)*(L-d)/(L*(L-1))) if L > d else (0, 0)
rng = np.random.default_rng(20260923)
M = [("C-x-C", "C", 2), ("C-x-x-C", "C", 3), ("H-x-H", "H", 2)]
def row(k, seqs):
    cells = []
    for _, r, d in M:
        a = np.array([oe(s, r, d) for s in seqs], float); n = len(seqs)
        bs = [a[ix, 0].sum()/a[ix, 1].sum() for ix in (rng.integers(0, n, n) for _ in range(2000)) if a[ix, 1].sum() > 0]
        cells.append(f"{a[:,0].sum()/a[:,1].sum():5.2f} [{np.percentile(bs,2.5):4.2f}-{np.percentile(bs,97.5):4.2f}] ({int(a[:,0].sum())})")
    tot = sum(map(len, seqs))
    print(f"{k:44s} {len(seqs):6d} {int(np.median([len(s) for s in seqs])):6d} " + " ".join(f"{c:>24s}" for c in cells)
          + f"  {100*sum(s.count('C') for s in seqs)/tot:5.2f} {100*sum(s.count('H') for s in seqs)/tot:5.2f}")
print(f"Heames orphan genes: {len(mech)}; de novo genes: {len(dn_genes)}; found in current FlyBase pep: "
      f"{sum(g in best for g in dn_genes)}")
print(f"{'group':44s} {'n':>6s} {'med aa':>6s} " + " ".join(f"{m[0]+' (obs pairs)':>24s}" for m in M) + "   Cys%  His%")
row("Heames DE NOVO (denovo + denovo-intron)", [best[g] for g in dn_genes if g in best])
for m in ["divergent", "putative"]:
    row(f"Heames {m} orphans (not de novo)", [best[g] for g in genes(m) - dn_genes if g in best])
row("NCBI-young, excluding Heames de novo", [best[g] for g in best if AGE.get(g, 9) <= 2 and g not in dn_genes])
row("NCBI-old (Insecta or older)", [best[g] for g in best if AGE.get(g, -1) >= 4])
row("neORFs (reference)", neo)

# small-n check for the de novo set: Poisson p (obs >= expected; ignores clustering) and how many genes carry the pairs
from scipy.stats import poisson
dn = [best[g] for g in dn_genes if g in best]
print("\nHeames de novo, small-n check:")
for lab, r, d in M:
    a = np.array([oe(s, r, d) for s in dn], float); o, e = a[:, 0].sum(), a[:, 1].sum()
    carriers = sorted(a[a[:, 0] > 0, 0].astype(int), reverse=True)
    print(f"  {lab:8s} obs {int(o):3d}  exp {e:5.1f}  Poisson P(>=obs) = {poisson.sf(o-1, e):.3f}   "
          f"genes carrying pairs: {len(carriers)} (pairs per gene {carriers})")
