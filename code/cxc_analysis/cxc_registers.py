#!/usr/bin/env python3
"""ADDENDUM: reading-frame/strand register prevalence of (CA/TG)n tracts inside ORFs.

A (CA)n microsatellite read on the coding strand appears as:
  CA-family text (CACAC / ACACA) -> His/Thr alternation -> H-x-H register (+ strand of tract)
  TG-family text (TGTGT / GTGTG) -> Cys/Val alternation -> C-x-C register (- strand of tract)
Frame only changes which residue LEADS (His-lead vs Thr-lead; Cys-lead vs Val-lead);
motif class is fixed by the coding-strand family. So the 6 (3 frame x 2 strand)
registers collapse to 2 motif classes, each with 2 distinct leading-residue phases
(a period-2 repeat in period-3 codons repeats every 6 nt -> frame0==frame2).

For each ORF (neORF coding strand/frame is defined; intergenic = every 6-frame
stop-to-stop ORF, both strands scanned as separate ORFs = the neutral 50/50 baseline)
we count embedded (CA/TG)n tracts (>= MINUNITS units) by register.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import os, sys, glob, gzip, time
import numpy as np
CODE = os.path.join(_HERE, "..")
sys.path.insert(0, CODE)
from hs_common import CODON2AA, read_fasta, load_denovo

SC   = f"{_HERE}/data/md_hp6_scratch"
DATA = f"{SC}/data"; GC = f"{SC}/gc2023"
CHRS = ["2L", "2R", "3L", "3R"]
MINUNITS = 4          # >=4 CA units = >=8 nt
MINMASK  = 2*MINUNITS - 2

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)

BASE2INT = np.full(256, 4, np.uint8)
for i, b in enumerate("ACGT"): BASE2INT[ord(b)] = i
AA_LUT = np.zeros(64, np.uint8)
for i0, a in enumerate("ACGT"):
    for i1, b in enumerate("ACGT"):
        for i2, c in enumerate("ACGT"):
            AA_LUT[i0*16+i1*4+i2] = ord(CODON2AA[a+b+c])
ORD_STAR = ord("*")
C,A_,G,T = (BASE2INT[ord(x)] for x in "CAGT")
CA_SET = frozenset((int(C), int(A_)))   # C,A
TG_SET = frozenset((int(T), int(G)))    # T,G

# families & phases -> register label
def classify_tract(codarr, start, length):
    """codarr: uint8 base ints (0..3), coding strand frame0. Return (motif, lead) or None."""
    b0 = int(codarr[start]); b1 = int(codarr[start+1])
    fam = {b0, b1}
    if fam == CA_SET:   motif = "H-x-H"     # His/Thr
    elif fam == TG_SET: motif = "C-x-C"     # Cys/Val
    else: return None
    p0 = start if start % 3 == 0 else start + (3 - start % 3)  # first full codon
    if p0 + 3 > codarr.size: p0 = start - (start % 3)          # fallback: prior boundary
    ci = int(codarr[p0])*16 + int(codarr[p0+1])*4 + int(codarr[p0+2])
    lead = chr(AA_LUT[ci])
    return motif, lead

def find_tracts(codarr):
    """Yield (start,length) maximal period-2 runs with >=MINMASK matches."""
    n = codarr.size
    if n < MINUNITS*2: return
    m = codarr[:-2] == codarr[2:]        # period-2 match
    # runs of True in m
    if not m.any(): return
    d = np.diff(m.astype(np.int8))
    st = list(np.flatnonzero(d == 1)+1); en = list(np.flatnonzero(d == -1)+1)
    if m[0]: st = [0]+st
    if m[-1]: en = en+[m.size]
    for s, e in zip(st, en):
        if e - s >= MINMASK:
            yield s, e - s + 2   # nt span

def tally(codarr, counts):
    for s, length in find_tracts(codarr):
        r = classify_tract(codarr, s, length)
        if r: counts[r] = counts.get(r, 0) + 1

# ---- neORF: coding strand = the CDS itself (frame 0) ----
def neorf_counts():
    paths = sorted(glob.glob(f"{GC}/*_DeNovoGenesNuc.fa"))
    den = load_denovo(paths, min_aa=10)
    counts = {}; totnt = 0; norf = len(den)
    for rec in den:
        dna = rec["dna"]
        if set(dna) - set("ACGT"):
            arr = BASE2INT[np.frombuffer(dna.encode(), np.uint8)]
            arr = arr[arr != 4]  # drop non-ACGT (rare); keeps period detection sane
        else:
            arr = BASE2INT[np.frombuffer(dna.encode(), np.uint8)]
        totnt += arr.size
        tally(arr, counts)
    return counts, totnt, norf

# ---- intergenic: every 6-frame stop-to-stop ORF, both strands ----
def load_intervals():
    genes = {c: [] for c in CHRS}; reps = {c: [] for c in CHRS}; chrset = set(CHRS)
    with gzip.open(f"{DATA}/dmel.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.split("\t")
            if len(f) < 5 or f[2] != "gene": continue
            if f[0] in chrset: genes[f[0]].append((int(f[3])-1, int(f[4])))
    with gzip.open(f"{DATA}/dm6_rmsk.txt.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12: continue
            gn = f[5]
            if gn.startswith("chr") and gn[3:] in chrset: reps[gn[3:]].append((int(f[6]), int(f[7])))
    return genes, reps
def mask_iv(iv, L):
    m = np.zeros(L, bool)
    for s, e in iv:
        if e <= 0 or s >= L: continue
        m[max(0,s):min(L,e)] = True
    return m
def runs_true(mask):
    if mask.size == 0: return
    d = np.diff(mask.astype(np.int8)); st = list(np.flatnonzero(d==1)+1); en = list(np.flatnonzero(d==-1)+1)
    if mask[0]: st=[0]+st
    if mask[-1]: en=en+[mask.size]
    for s,e in zip(st,en): yield s,e

def scan_region_orfs(dna_bytes, counts, tot):
    x = BASE2INT[dna_bytes]
    comp = np.where(x == 4, 4, 3-x).astype(np.uint8)
    for src in (x, comp[::-1]):
        n = src.size
        for frame in range(3):
            m = (n-frame)//3
            if m < 10: continue
            cod = src[frame:frame+3*m].reshape(m,3).astype(np.int64)
            bad = (cod==4).any(axis=1)
            ci = cod[:,0]*16+cod[:,1]*4+cod[:,2]; ci[bad]=0
            aa = AA_LUT[ci].copy(); aa[bad]=ORD_STAR
            is_stop = aa == ORD_STAR
            idx = np.flatnonzero(is_stop)
            b = np.empty(idx.size+2, np.int64); b[0]=-1; b[-1]=m; b[1:-1]=idx
            codarr_full = src[frame:frame+3*m]  # coding nt, frame0 aligned to codons
            for k in range(b.size-1):
                s=b[k]+1; e=b[k+1]
                if e-s < 10: continue
                sub = codarr_full[3*s:3*e]
                tot[0] += sub.size
                if (sub[:-2]==sub[2:]).any():
                    tally(sub, counts)

def intergenic_counts():
    genes, reps = load_intervals()
    counts = {}; tot=[0]
    for c in CHRS:
        seq=None
        for h,s in read_fasta(f"{DATA}/dna_{c}.fa.gz"): seq=s.upper(); break
        Lc=len(seq); sb=np.frombuffer(seq.encode(),np.uint8)
        gm=mask_iv(genes[c],Lc); rm=mask_iv(reps[c],Lc)
        clean=(~gm)&(~rm)
        for s,e in runs_true(clean):
            if e-s>=30: scan_region_orfs(sb[s:e], counts, tot)
        log(f"  chr{c} done, tracts so far={sum(counts.values()):,}")
    return counts, tot[0]

def report(name, counts, totnt, norf=None):
    order=[("H-x-H","H"),("H-x-H","T"),("C-x-C","C"),("C-x-C","V")]
    labels={("H-x-H","H"):"His-lead (H-x-H)",("H-x-H","T"):"Thr-lead (H-x-H)",
            ("C-x-C","C"):"Cys-lead (C-x-C)",("C-x-C","V"):"Val-lead (C-x-C)"}
    tot=sum(counts.values())
    hxh=sum(v for k,v in counts.items() if k[0]=="H-x-H")
    cxc=sum(v for k,v in counts.items() if k[0]=="C-x-C")
    other=tot-hxh-cxc
    mb=totnt/1e6
    print(f"\n### {name}: total (CA/TG)n tracts>= {MINUNITS} units inside ORFs = {tot}"
          + (f"  ({norf} ORFs," if norf else "  (") + f" {mb:.2f} Mb coding)")
    print(f"  {'register':22s} {'count':>8s} {'% of tracts':>12s} {'per Mb':>10s}")
    for key in order:
        v=counts.get(key,0)
        print(f"  {labels[key]:22s} {v:8d} {100*v/tot if tot else 0:11.2f}% {v/mb:10.2f}")
    print(f"  {'-- H-x-H (His) class':22s} {hxh:8d} {100*hxh/tot if tot else 0:11.2f}% {hxh/mb:10.2f}")
    print(f"  {'-- C-x-C (Cys) class':22s} {cxc:8d} {100*cxc/tot if tot else 0:11.2f}% {cxc/mb:10.2f}")
    if other: print(f"  (other period-2 repeats, non CA/TG, excluded from classes: {other})")
    return dict(hxh=hxh,cxc=cxc,tot=tot,mb=mb,hxh_permb=hxh/mb,cxc_permb=cxc/mb,
                his_frac=hxh/(hxh+cxc) if (hxh+cxc) else float('nan'))

def main():
    log("neORF tracts ..."); cn, ntn, norf = neorf_counts()
    log("intergenic tracts (genome pass) ..."); ci, nti = intergenic_counts()
    dn=report("neORF", cn, ntn, norf)
    ig=report("intergenic ORFs (both strands = neutral baseline)", ci, nti)
    print("\n=== DISCRIMINATOR ===")
    print(f"His-register fraction of CA/TG tracts:  neORF={100*dn['his_frac']:.1f}%   intergenic={100*ig['his_frac']:.1f}%   neutral=50.0%")
    print(f"H-x-H (His) tracts per Mb: neORF={dn['hxh_permb']:.1f}  intergenic={ig['hxh_permb']:.1f}   ratio={dn['hxh_permb']/ig['hxh_permb']:.2f}x")
    print(f"C-x-C (Cys) tracts per Mb: neORF={dn['cxc_permb']:.1f}  intergenic={ig['cxc_permb']:.1f}   ratio={dn['cxc_permb']/ig['cxc_permb']:.2f}x")
    # 2x2 chi-square: (His,Cys) x (neORF,intergenic)
    from scipy.stats import chi2_contingency, fisher_exact
    tab=np.array([[dn['hxh'],dn['cxc']],[ig['hxh'],ig['cxc']]])
    chi2,p,_,_=chi2_contingency(tab)
    print(f"His:Cys ratio neORF vs intergenic: chi2={chi2:.2f} p={p:.2e}  (table {tab.tolist()})")

if __name__=="__main__":
    main()
