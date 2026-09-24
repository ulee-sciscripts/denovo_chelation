#!/usr/bin/env python3
"""C-x-C register-asymmetry analysis, reusing the paper's His-spacing machinery.

Order-fold = observed_spacing / analytic-scramble-expected  (== MC limit),
identical method to denovo_hisspacing_analysis.py fig1_profile.
Spacing convention: 'n' = number of spacer residues between the two residues.
  n=1 spacer  ->  X-x-X    ->  count_spacings offset d = 2
  n=k spacers ->  d = k+1
We compute d=1..7 so n=0..6 are covered; report n=1..6 (d=2..7).
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import os, sys, glob, gzip, time
import numpy as np

CODE = os.path.join(_HERE, "..")
sys.path.insert(0, CODE)
from hs_common import CODON2AA, AAS, read_fasta, translate, load_denovo, load_proteins

SC   = f"{_HERE}/data/md_hp6_scratch"
DATA = f"{SC}/data"
GC   = f"{SC}/gc2023"
SEED = 20260812
DMAX = 7
CHRS = ["2L", "2R", "3L", "3R"]

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

# ---- per-peptide spacing rows (L,h,o1..oDMAX) for a given residue ----
def rows_from_prots(prots, res):
    R = ord(res); rows = []
    for p in prots:
        aa = np.frombuffer(p.encode(), np.uint8)
        L = aa.size
        if L < 10: continue
        m = aa == R; h = int(m.sum())
        if h < 2: continue
        o = np.zeros(DMAX, np.int64)
        for d in range(1, DMAX+1):
            if L > d: o[d-1] = int(np.count_nonzero(m[:-d] & m[d:]))
        rows.append((L, h, *o.tolist()))
    return np.array(rows, np.int64) if rows else np.zeros((0, 2+DMAX), np.int64)

def analytic_exp_rows(tab, d):
    L = tab[:,0].astype(float); h = tab[:,1].astype(float)
    return np.where(L > d, (L-d)*h*(h-1)/(L*(L-1)), 0.0)

def pooled_fold(tab, d):
    o = tab[:,1+d].astype(float).sum(); e = analytic_exp_rows(tab, d).sum()
    return (o/e if e > 0 else np.nan), o, e

def boot_fold(tab, d, nboot, rng):
    n = tab.shape[0]
    if n == 0: return np.full(nboot, np.nan)
    o = tab[:,1+d].astype(float); e = analytic_exp_rows(tab, d)
    out = np.empty(nboot)
    for b in range(nboot):
        ix = rng.integers(0, n, n); se = e[ix].sum()
        out[b] = (o[ix].sum()/se) if se > 0 else np.nan
    return out

def boot_ci(tab, d, nboot=3000, seed=SEED):
    b = boot_fold(tab, d, nboot, np.random.default_rng(seed)); b = b[~np.isnan(b)]
    return [np.nan, np.nan] if b.size == 0 else np.percentile(b, [2.5, 97.5]).tolist()

def p_vs1(tab, d, nboot=4000, seed=SEED):
    b = boot_fold(tab, d, nboot, np.random.default_rng(seed)); b = b[~np.isnan(b)]
    if b.size == 0: return np.nan
    p = 2*min(np.mean(b <= 1), np.mean(b >= 1))
    return min(1.0, max(p, 1.0/(b.size+1)))

def boot_diff_p(tabA, tabB, d, nboot=3000, seed=SEED):
    rng = np.random.default_rng(seed)
    fa = boot_fold(tabA, d, nboot, rng); fb = boot_fold(tabB, d, nboot, rng)
    diff = fa - fb; diff = diff[~np.isnan(diff)]
    if diff.size == 0: return np.nan
    p = 2*min(np.mean(diff <= 0), np.mean(diff >= 0))
    return min(1.0, max(p, 1.0/(diff.size+1)))

# ---- intergenic scan (one genome pass; accumulate H and C rows) ----
class Accum:
    def __init__(self):
        self.rowsH = []; self.rowsC = []
        self.compH = 0; self.compC = 0; self.compTot = 0
        self.norf = 0
        self.nH = 0   # ORFs (>=10aa) with >=1 H-x-H (d=2)
        self.nC = 0   # ORFs (>=10aa) with >=1 C-x-C (d=2)
    def tabH(self): return np.array(self.rowsH, np.int64) if self.rowsH else np.zeros((0,2+DMAX),np.int64)
    def tabC(self): return np.array(self.rowsC, np.int64) if self.rowsC else np.zeros((0,2+DMAX),np.int64)

ORD_H, ORD_C = ord("H"), ord("C")
def _scan_frame(aa, acc):
    is_stop = aa == ORD_STAR
    nonstop = aa[~is_stop]
    if nonstop.size:
        acc.compTot += nonstop.size
        acc.compH += int(np.count_nonzero(nonstop == ORD_H))
        acc.compC += int(np.count_nonzero(nonstop == ORD_C))
    idx = np.flatnonzero(is_stop)
    bounds = np.empty(idx.size+2, np.int64); bounds[0] = -1; bounds[-1] = aa.size; bounds[1:-1] = idx
    for k in range(bounds.size-1):
        s = bounds[k]+1; e = bounds[k+1]; L = e-s
        if L < 10: continue
        acc.norf += 1
        seg = aa[s:e]
        for R, rows, cntattr in ((ORD_H,'rowsH','nH'), (ORD_C,'rowsC','nC')):
            m = seg == R; h = int(m.sum())
            if h < 2: continue
            o = np.zeros(DMAX, np.int64)
            for d in range(1, DMAX+1):
                if L > d: o[d-1] = int(np.count_nonzero(m[:-d] & m[d:]))
            getattr(acc, rows).append((L, h, *o.tolist()))
            if o[1] > 0:  # d=2 == X-x-X, n=1 spacer
                setattr(acc, cntattr, getattr(acc, cntattr)+1)

def scan_region_seq(dna_bytes, acc):
    x = BASE2INT[dna_bytes]
    comp = np.where(x == 4, 4, 3-x).astype(np.uint8)
    for src in (x, comp[::-1]):
        n = src.size
        for frame in range(3):
            m = (n-frame)//3
            if m < 10: continue
            cod = src[frame:frame+3*m].reshape(m,3).astype(np.int64)
            bad = (cod == 4).any(axis=1)
            ci = cod[:,0]*16+cod[:,1]*4+cod[:,2]; ci[bad] = 0
            aa = AA_LUT[ci].copy(); aa[bad] = ORD_STAR
            _scan_frame(aa, acc)

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

def mask_from_intervals(intervals, L):
    m = np.zeros(L, bool)
    for s, e in intervals:
        if e <= 0 or s >= L: continue
        m[max(0,s):min(L,e)] = True
    return m

def runs_of_true(mask):
    if mask.size == 0: return
    d = np.diff(mask.astype(np.int8))
    starts = list(np.flatnonzero(d == 1)+1); ends = list(np.flatnonzero(d == -1)+1)
    if mask[0]: starts = [0]+starts
    if mask[-1]: ends = ends+[mask.size]
    for s, e in zip(starts, ends): yield s, e

def build_intergenic():
    genes, reps = load_intervals()
    clean = Accum(); repint = Accum()
    for c in CHRS:
        seq = None
        for hdr, s in read_fasta(f"{DATA}/dna_{c}.fa.gz"): seq = s.upper(); break
        Lc = len(seq); seqb = np.frombuffer(seq.encode(), np.uint8)
        gmask = mask_from_intervals(genes[c], Lc); rmask = mask_from_intervals(reps[c], Lc)
        clean_mask = (~gmask) & (~rmask); repint_mask = rmask & (~gmask)
        for s, e in runs_of_true(clean_mask):
            if e-s >= 30: scan_region_seq(seqb[s:e], clean)
        for s, e in runs_of_true(repint_mask):
            if e-s >= 30: scan_region_seq(seqb[s:e], repint)
        log(f"  chr{c}: cleanORF={clean.norf:,} repintORF={repint.norf:,}")
    return clean, repint

# ---- incidence + composition helpers for peptide sets ----
def incidence_and_comp(prots):
    nH = nC = norf = 0; cH = cC = ctot = 0
    for p in prots:
        aa = np.frombuffer(p.encode(), np.uint8); L = aa.size
        if L < 10: continue
        norf += 1
        ctot += L; cH += int(np.count_nonzero(aa == ORD_H)); cC += int(np.count_nonzero(aa == ORD_C))
        for R, flag in ((ORD_H,'H'), (ORD_C,'C')):
            m = aa == R
            if L > 2 and np.count_nonzero(m[:-2] & m[2:]) > 0:
                if flag == 'H': nH += 1
                else: nC += 1
    return dict(norf=norf, incH=nH/norf, incC=nC/norf,
                cyspct=100*cC/ctot, hispct=100*cH/ctot)

def spectrum(tab, label):
    log(f"--- {label}: n_rows(>=2 res,>=10aa)={tab.shape[0]:,}")
    res = {}
    for d in range(1, DMAX+1):
        f, o, e = pooled_fold(tab, d)
        res[d] = dict(fold=f, o=o, e=e)
    return res

def main():
    paths = sorted(glob.glob(f"{GC}/*_DeNovoGenesNuc.fa"))
    den = load_denovo(paths, min_aa=10)
    prots_dn = [d["prot"] for d in den]
    log(f"neORF peptides: {len(prots_dn)}")
    prot = load_proteins(f"{DATA}/dmel_pep.fa.gz", min_aa=10)
    meta = load_proteins(f"{DATA}/dmel_metalbinding.fasta", min_aa=10)
    log(f"proteome {len(prot)}  metal {len(meta)}")

    tabH_dn = rows_from_prots(prots_dn, "H"); tabC_dn = rows_from_prots(prots_dn, "C")
    tabH_pr = rows_from_prots(prot, "H");     tabC_pr = rows_from_prots(prot, "C")
    tabH_mt = rows_from_prots(meta, "H");     tabC_mt = rows_from_prots(meta, "C")

    log("scanning intergenic ...")
    clean, repint = build_intergenic()
    tabH_cl = clean.tabH(); tabC_cl = clean.tabC()
    log(f"intergenic clean ORFs>=10aa={clean.norf:,}  H-rows={tabH_cl.shape[0]:,}  C-rows={tabC_cl.shape[0]:,}")

    # ========= DELIVERABLE 2: full spectra n=1..6 (d=2..7) =========
    specH_dn = spectrum(tabH_dn, "neORF H"); specC_dn = spectrum(tabC_dn, "neORF C")
    specH_cl = spectrum(tabH_cl, "interg H"); specC_cl = spectrum(tabC_cl, "interg C")
    specH_pr = spectrum(tabH_pr, "prot H");   specC_pr = spectrum(tabC_pr, "prot C")

    print("\n=== SPECTRUM: fold vs analytic scramble null (order-fold) ===")
    print("n=spacers  d   | H neORF  H interg  H prot | C neORF  C interg  C prot")
    for n in range(0, 7):
        d = n+1
        print(f"  n={n:<2} d={d:<2} | "
              f"{specH_dn[d]['fold']:7.3f} {specH_cl[d]['fold']:8.3f} {specH_pr[d]['fold']:7.3f} | "
              f"{specC_dn[d]['fold']:7.3f} {specC_cl[d]['fold']:8.3f} {specC_pr[d]['fold']:7.3f}")

    # ========= DELIVERABLE 1: C-x-C (n=1,d=2) with CI + p =========
    d = 2
    print("\n=== DELIVERABLE 1: C-x-C (n=1 spacer, d=2) ===")
    for lab, tab in (("neORF", tabC_dn), ("intergenic", tabC_cl), ("proteome", tabC_pr), ("metal", tabC_mt)):
        f, o, e = pooled_fold(tab, d)
        ci = boot_ci(tab, d, seed=SEED+11)
        pv = p_vs1(tab, d, seed=SEED+21)
        print(f"  C-x-C {lab:11s}: fold={f:.3f}  CI[{ci[0]:.3f},{ci[1]:.3f}]  p(vs1)={pv:.2e}  (o={o:.0f}, e={e:.1f}, rows={tab.shape[0]})")
    p_dnvcl_C = boot_diff_p(tabC_dn, tabC_cl, d, seed=SEED+31)
    p_dnvpr_C = boot_diff_p(tabC_dn, tabC_pr, d, seed=SEED+41)
    print(f"  C-x-C p(neORF != intergenic) = {p_dnvcl_C:.2e}")
    print(f"  C-x-C p(neORF != proteome)   = {p_dnvpr_C:.2e}")

    print("\n=== Reference: H-x-H (n=1 spacer, d=2) ===")
    for lab, tab in (("neORF", tabH_dn), ("intergenic", tabH_cl), ("proteome", tabH_pr), ("metal", tabH_mt)):
        f, o, e = pooled_fold(tab, d)
        ci = boot_ci(tab, d, seed=SEED+12)
        pv = p_vs1(tab, d, seed=SEED+22)
        print(f"  H-x-H {lab:11s}: fold={f:.3f}  CI[{ci[0]:.3f},{ci[1]:.3f}]  p(vs1)={pv:.2e}")
    p_dnvcl_H = boot_diff_p(tabH_dn, tabH_cl, d, seed=SEED+32)
    print(f"  H-x-H p(neORF != intergenic) = {p_dnvcl_H:.2e}")

    # C-x-x-C (n=2 spacer, d=3) reference
    print("\n=== Reference: C-x-x-C (n=2 spacers, d=3) ===")
    for lab, tab in (("neORF", tabC_dn), ("intergenic", tabC_cl), ("proteome", tabC_pr), ("metal", tabC_mt)):
        f, o, e = pooled_fold(tab, 3); ci = boot_ci(tab, 3, seed=SEED+13)
        print(f"  C-x-x-C {lab:11s}: fold={f:.3f}  CI[{ci[0]:.3f},{ci[1]:.3f}]")

    # ========= DELIVERABLE 3: incidence =========
    inc_dn = incidence_and_comp(prots_dn)
    inc_pr = incidence_and_comp(prot)
    inc_cl = dict(norf=clean.norf, incH=clean.nH/clean.norf, incC=clean.nC/clean.norf,
                  cyspct=100*clean.compC/clean.compTot, hispct=100*clean.compH/clean.compTot)
    print("\n=== DELIVERABLE 3: raw incidence (fraction of ORFs with >=1 motif) ===")
    print(f"  neORF:      >=1 H-x-H = {inc_dn['incH']*100:.2f}%   >=1 C-x-C = {inc_dn['incC']*100:.2f}%   (n={inc_dn['norf']})")
    print(f"  intergenic: >=1 H-x-H = {inc_cl['incH']*100:.2f}%   >=1 C-x-C = {inc_cl['incC']*100:.2f}%   (n={inc_cl['norf']})")
    print(f"  proteome:   >=1 H-x-H = {inc_pr['incH']*100:.2f}%   >=1 C-x-C = {inc_pr['incC']*100:.2f}%   (n={inc_pr['norf']})")

    # ========= DELIVERABLE 4: baseline composition =========
    print("\n=== DELIVERABLE 4: baseline residue composition ===")
    print(f"  neORF:      Cys%={inc_dn['cyspct']:.3f}   His%={inc_dn['hispct']:.3f}")
    print(f"  intergenic: Cys%={inc_cl['cyspct']:.3f}   His%={inc_cl['hispct']:.3f}")
    print(f"  proteome:   Cys%={inc_pr['cyspct']:.3f}   His%={inc_pr['hispct']:.3f}")

if __name__ == "__main__":
    main()
