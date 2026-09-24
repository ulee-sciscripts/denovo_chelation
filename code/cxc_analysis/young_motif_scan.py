#!/usr/bin/env python3
"""Test 2: do young (Drosophila-restricted), unannotated, cytosolic, Cys-rich short proteins carry
metal-binding spacing beyond composition? Order-fold = observed / within-sequence scramble expectation
(analytic, as in the paper) for C-x-C, C-x-x-C, H-x-H, mixed C/H at d=3 (C-x-x-H / H-x-x-C);
'paired C-x-x-C' (>=2 CxxC with starts 8-40 aa apart; zinc-ribbon/treble-clef-like) vs 200 shuffles.
Groups: young unannotated (test), zinc-domain Cys-rich (positive control), secreted Cys-rich short
(disulfide comparator), neORFs (reference). Bootstrap 95% CI over proteins."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, json, re, runpy, io, contextlib
import numpy as np
sys.path.insert(0, os.path.join(_HERE, ".."))
from hs_common import read_fasta, translate
D = _HERE
with contextlib.redirect_stdout(io.StringIO()): g = runpy.run_path(f"{D}/code/cys_short_zinc.py")
rows, hi = g["rows"], g["hi"]
AGE = json.load(open(f"{D}/data/gene_age_ncbi.json"))
ZN2 = re.compile(g["ZN"].pattern + r"|LITAF|LPS-induced tumour|Yippee|Transcription elongation factor 1|Siva|"
                 r"Zinc-binding ribosomal|eL37|eL40|eS27|eS31|zinc ribbon", re.I)
q = np.percentile([r["cp"] for r in rows], 66.7)
def uniq(rs):
    o = {}
    for r in rs: o.setdefault(r["s"], r)
    return [r["s"] for r in o.values()]
G = {"young unannotated cyto Cys-rich": uniq([r for r in hi if not r["ipr"] and AGE.get(r["g"], 9) <= 2]),
     "zinc-domain Cys-rich (+ ctrl)":  uniq([r for r in hi if any(ZN2.search(x) for x in r["ipr"])]),
     "secreted Cys-rich short":        uniq([r for r in rows if r["cp"] > q and r["sec"]])}
seen = set(); neo = []
for h, s in read_fasta(f"{_HERE}/../../data/analysis/neORF.fa"):
    p, _ = translate(s.upper(), truncate_at_stop=True)
    if len(p) >= 10 and p not in seen: seen.add(p); neo.append(p)
G["neORFs (reference)"] = neo

def pair_oe(s, a, b, d):
    L = len(s)
    if L <= d: return 0.0, 0.0
    x = np.frombuffer(s.encode(), np.uint8)
    ma, mb = x == ord(a), x == ord(b)
    o = np.count_nonzero(ma[:L-d] & mb[d:]); na, nb = ma.sum(), mb.sum()
    e = (na * (na - 1) if a == b else na * nb) * (L - d) / (L * (L - 1))
    if a != b:   # count both orientations (a..b and b..a)
        o += np.count_nonzero(mb[:L-d] & ma[d:]); e *= 2
    return float(o), float(e)
CXXC = re.compile(r"(?=C..C)")
def paired(s):
    st = [m.start() for m in CXXC.finditer(s)]
    return int(any(8 <= b - a <= 40 for i, a in enumerate(st) for b in st[i+1:]))
rng = np.random.default_rng(20260923)
def paired_oe(seqs, nsh=200):
    o = np.array([paired(s) for s in seqs], float); e = np.zeros(len(seqs))
    for i, s in enumerate(seqs):
        if s.count("C") < 4: continue
        a = np.frombuffer(s.encode(), np.uint8).copy()
        e[i] = np.mean([paired(rng.permutation(a).tobytes().decode()) for _ in range(nsh)])
    return o, e
def boot(o, e, B=2000):
    n = len(o); f = o.sum() / e.sum() if e.sum() else np.nan
    bs = [o[ix].sum() / e[ix].sum() for ix in (rng.integers(0, n, n) for _ in range(B)) if e[ix].sum() > 0]
    return f, np.percentile(bs, 2.5), np.percentile(bs, 97.5)
M = [("C-x-C", "C", "C", 2), ("C-x-x-C", "C", "C", 3), ("H-x-H", "H", "H", 2), ("C/H d=3 (CxxH/HxxC)", "C", "H", 3)]
print(f"{'group':34s} {'n':>5s} " + " ".join(f"{m[0]:>22s}" for m in M) + f" {'paired C-x-x-C':>22s}")
for k, seqs in G.items():
    cells = []
    for lab, a, b, d in M:
        oe = np.array([pair_oe(s, a, b, d) for s in seqs])
        f, lo, hi_ = boot(oe[:, 0], oe[:, 1]); cells.append(f"{f:5.2f} [{lo:4.2f}-{hi_:4.2f}]")
    if k.startswith("neORF"):
        cells.append("(skipped)")
    else:
        o, e = paired_oe(seqs); f, lo, hi_ = boot(o, e)
        cells.append(f"{f:5.2f} [{lo:4.2f}-{hi_:4.2f}] ({int(o.sum())})")
    cys = 100 * sum(s.count("C") for s in seqs) / sum(map(len, seqs)); his = 100 * sum(s.count("H") for s in seqs) / sum(map(len, seqs))
    print(f"{k:34s} {len(seqs):5d} " + " ".join(f"{c:>22s}" for c in cells) + f"   Cys {cys:.1f}% His {his:.1f}%")
