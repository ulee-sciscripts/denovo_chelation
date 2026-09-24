#!/usr/bin/env python3
"""(CA)n vs (TG)n tracts on the coding strand (>=4 units; cxc_registers.py definition) in
neORFs, clean intergenic ORFs, the Fig 1 length-matched proteome draw, and the whole proteome.
Proteome CDSs: zinc_analysis/ref/dmel_cds.fa.gz, matched to the canonical dmel_pep.fa.gz proteins
by translated sequence. (CA)n on coding strand -> His/Thr (H-x-H); (TG)n -> Cys/Val (C-x-C)."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, glob, gzip
import numpy as np
from scipy.stats import fisher_exact
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))
import cxc_registers as X
from hs_common import load_proteins, read_fasta, translate
SEED = 20260812
S = f"{_HERE}/data/md_hp6_scratch"

def count_cds(cds_list):
    counts, tot = {}, 0
    for dna in cds_list:
        arr = X.BASE2INT[np.frombuffer(dna.encode(), np.uint8)]; arr = arr[arr != 4]
        tot += arr.size; X.tally(arr, counts)
    return counts, tot

cds_of = {}
for h, s in read_fasta(f"{_HERE}/../../data/analysis/ref/dmel_cds.fa.gz"):
    p, c = translate(s.upper(), truncate_at_stop=True)
    if len(p) >= 10: cds_of.setdefault(p, c)
prot = load_proteins(f"{S}/data/dmel_pep.fa.gz", min_aa=10)
# Fig 1 length-matched draw (verbatim method, gc2023 neORF lengths)
den = X.load_denovo(sorted(glob.glob(f"{S}/gc2023/*_DeNovoGenesNuc.fa")), min_aa=10)
dl = np.array([len(d["prot"]) for d in den])
plen = np.array([len(s) for s in prot]); order = np.argsort(plen, kind="stable")
ps = [prot[i] for i in order]; pls = plen[order]; rng = np.random.default_rng(SEED); matched = []
for L in dl:
    j = int(np.searchsorted(pls, int(L)))
    cand = [k for k in (j - 1, j, j + 1) if 0 <= k < len(ps)]
    d = np.array([abs(int(pls[k]) - int(L)) for k in cand])
    matched.append(ps[int(rng.choice([cand[t] for t in range(len(cand)) if d[t] == d.min()]))])
hxh = np.mean([bool(np.any((np.frombuffer(s.encode(), np.uint8)[:-2] == 72) & (np.frombuffer(s.encode(), np.uint8)[2:] == 72))) for s in matched])
print(f"anchor: matched-proteome H-x-H incidence {100*hxh:.2f}% (Fig 1: 4.39%)")
uniq = list(dict.fromkeys(prot))
miss_w = sum(p not in cds_of for p in uniq); miss_m = sum(p not in cds_of for p in matched)
print(f"CDS found: whole proteome {len(uniq)-miss_w}/{len(uniq)}, matched {len(matched)-miss_m}/{len(matched)}")

res = {}
X.log("neORF"); res["neORFs"] = X.neorf_counts()[:2]
X.log("intergenic"); res["clean intergenic ORFs"] = X.intergenic_counts()
res["length-matched proteome"] = count_cds([cds_of[p] for p in matched if p in cds_of])
res["whole proteome"] = count_cds([cds_of[p] for p in uniq if p in cds_of])

print(f"\n{'set':26s} {'coding Mb':>9s} {'(CA)n':>7s} {'(TG)n':>7s} {'CA:TG':>6s} {'% CA':>6s} {'CA/Mb':>7s} {'TG/Mb':>7s}   Fisher CA:TG vs intergenic")
ig = res["clean intergenic ORFs"][0]
igc = (sum(v for k, v in ig.items() if k[0] == "H-x-H"), sum(v for k, v in ig.items() if k[0] == "C-x-C"))
for k, (c, tot) in res.items():
    ca = sum(v for kk, v in c.items() if kk[0] == "H-x-H"); tg = sum(v for kk, v in c.items() if kk[0] == "C-x-C")
    mb = tot / 1e6
    fx = "" if k.startswith("clean") else "OR %.2f p=%.1e" % fisher_exact([[ca, tg], list(igc)])
    print(f"{k:26s} {mb:9.2f} {ca:7d} {tg:7d} {ca/tg:6.2f} {100*ca/(ca+tg):5.1f}% {ca/mb:7.1f} {tg/mb:7.1f}   {fx}")

# resampling check: matched draw reuses proteins; recount on distinct proteins and list top contributors
import collections
um = list(dict.fromkeys(p for p in matched if p in cds_of))
c, tot = count_cds([cds_of[p] for p in um])
ca = sum(v for k, v in c.items() if k[0] == "H-x-H"); tg = sum(v for k, v in c.items() if k[0] == "C-x-C")
print(f"\nlength-matched, DISTINCT proteins only: n={len(um)} of {len(matched)} draws  (CA)n {ca}  (TG)n {tg}")
per = collections.Counter()
for p in matched:
    if p in cds_of:
        cc, _ = count_cds([cds_of[p]]); per[p] += sum(v for k, v in cc.items() if k[0] == "C-x-C")
print("top (TG)n contributors among draws (tracts, length, times drawn):")
mc = collections.Counter(matched)
for p, v in per.most_common(6): print(f"  {v:3d}  len {len(p)}  drawn {mc[p]}x  {p[:40]}...")

# length-restricted proteome (distinct proteins in neORF length range) as the usable length control
for lo, hi in [(30, 150), (10, 150)]:
    sub = [p for p in uniq if lo <= len(p) <= hi and p in cds_of]
    c, tot = count_cds([cds_of[p] for p in sub])
    ca = sum(v for k, v in c.items() if k[0] == "H-x-H"); tg = sum(v for k, v in c.items() if k[0] == "C-x-C")
    print(f"proteome {lo}-{hi} aa distinct n={len(sub)}  {tot/1e6:.2f} Mb  (CA)n {ca}  (TG)n {tg}  "
          f"CA:TG {ca/tg if tg else float('nan'):.2f}  CA/Mb {ca/(tot/1e6):.1f}  TG/Mb {tg/(tot/1e6):.1f}  "
          f"Fisher vs intergenic OR {fisher_exact([[ca, tg], list(igc)])[0]:.2f} p={fisher_exact([[ca, tg], list(igc)])[1]:.1e}")
