#!/usr/bin/env python3
"""(CA)n fingerprint of H-x-H motifs: neORF vs clean-intergenic.
A (CA)n-derived H-x-H reads His-Thr-His with His codons CAC (CAC.ACA.CAC).
So per set, among H-x-H motifs, measure: Thr-middle fraction, CAC-codon fraction,
and full (CA)n signature (Thr-middle AND both His = CAC). Plus His density.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import sys, csv
import numpy as np
from scipy import stats
CODE=os.path.join(_HERE, "..")
sys.path.insert(0, CODE)
import denovo_hisspacing_analysis as A
from hs_common import read_fasta, translate
TMP=sys.argv[1]; A.SCRATCH=TMP
GENOME=f"{_HERE}/../../data/analysis/ref/dm6_genome.fa"
NEORF=f"{_HERE}/../../data/analysis/neORF.fa"
CHRS=["2L","2R","3L","3R"]
B2I,LUT,H,STAR=A.BASE2INT,A.AA_LUT,A.ORD_H,A.ORD_STAR
T=ord('T')

def summarize(name, n_hxh, thr_mid, cac, cat, tot_his, tot_len, full_can):
    thr_f = thr_mid/n_hxh if n_hxh else 0
    cac_f = cac/(cac+cat) if (cac+cat) else 0
    full_f = full_can/n_hxh if n_hxh else 0
    print(f"\n== {name} ==")
    print(f"  H-x-H motifs: {n_hxh:,}   His density: {100*tot_his/tot_len:.2f} His/100aa")
    print(f"  Thr-middle fraction : {100*thr_f:5.1f}%   (His-Thr-His; the (CA)n reading)")
    print(f"  CAC-codon fraction  : {100*cac_f:5.1f}%   (of His in H-x-H; (CA)n forces CAC)")
    print(f"  FULL (CA)n signature: {100*full_f:5.1f}%   (Thr-middle AND both His=CAC)")
    return thr_f, cac_f, full_f

# ---------- neORFs ----------
seen=set(); n_hxh=thr=cac=cat=full=this=0; tot_his=tot_len=0
for hdr,seq in read_fasta(NEORF):
    p,dna=translate(seq.upper(),truncate_at_stop=True)
    if len(p)<10 or p in seen: continue
    seen.add(p); L=len(p)
    a=np.frombuffer(p.encode(),np.uint8); m=a==H
    tot_his+=int(m.sum()); tot_len+=L
    if L>2:
        for i in np.flatnonzero(m[:-2]&m[2:]):
            n_hxh+=1
            if p[i+1]=='T': thr+=1
            cA=dna[3*i:3*i+3]=="CAC"; cB=dna[3*(i+2):3*(i+2)+3]=="CAC"
            cac+= cA+cB; cat+= (dna[3*i:3*i+3]=="CAT")+(dna[3*(i+2):3*(i+2)+3]=="CAT")
            if p[i+1]=='T' and cA and cB: full+=1
neo=summarize("neORFs", n_hxh, thr, cac, cat, tot_his, tot_len, full)
NEO=(n_hxh, thr, cac, cat)   # (H-x-H motifs, Thr-middle, His-CAC, His-CAT)

# ---------- clean intergenic (re-scan, tracking codons at H-x-H) ----------
seqs={}
for h,s in read_fasta(GENOME):
    c=h.split()[0]
    if c in CHRS: seqs[c]=s.upper()
genes,reps=A.load_intervals()
n_hxh=thr=cac=cat=full=0; tot_his=tot_len=0
def scan(seqb):
    global n_hxh,thr,cac,cat,full,tot_his,tot_len
    x=B2I[seqb]; xr=np.where(x==4,4,3-x).astype(np.uint8)[::-1]
    for src in (x,xr):
        n=src.size
        for fr in range(3):
            mm=(n-fr)//3
            if mm<10: continue
            cod=src[fr:fr+3*mm].reshape(mm,3).astype(np.int64)
            valid=(cod<4).all(axis=1)
            idx=np.where(valid,cod[:,0]*16+cod[:,1]*4+cod[:,2],0)
            aa=np.where(valid,LUT[idx],STAR).astype(np.uint8)
            stop=np.flatnonzero(aa==STAR); bnd=np.concatenate(([-1],stop,[aa.size]))
            for k in range(bnd.size-1):
                s0,e0=bnd[k]+1,bnd[k+1]; Ln=e0-s0
                if Ln<10: continue
                seg=aa[s0:e0]; hm=seg==H
                tot_his+=int(hm.sum()); tot_len+=Ln
                if Ln>2:
                    for li in np.flatnonzero(hm[:-2]&hm[2:]):
                        gi=s0+li; n_hxh+=1
                        if aa[gi+1]==T: thr+=1
                        cA=(cod[gi,0]==1 and cod[gi,1]==0 and cod[gi,2]==1)
                        cB=(cod[gi+2,0]==1 and cod[gi+2,1]==0 and cod[gi+2,2]==1)
                        tA=(cod[gi,0]==1 and cod[gi,1]==0 and cod[gi,2]==3)
                        tB=(cod[gi+2,0]==1 and cod[gi+2,1]==0 and cod[gi+2,2]==3)
                        # int(): cA..tB are numpy bools here; np.bool_ '+' is logical-OR,
                        # not arithmetic, so a both-CAC motif would wrongly add 1 not 2.
                        cac+=int(cA)+int(cB); cat+=int(tA)+int(tB)
                        if aa[gi+1]==T and cA and cB: full+=1
for c in CHRS:
    seqb=np.frombuffer(seqs[c].encode(),np.uint8); Lc=len(seqs[c])
    clean=(~A.mask_from_intervals(genes[c],Lc))&(~A.mask_from_intervals(reps[c],Lc))
    for s,e in A.runs_of_true(clean):
        if e-s>=30: scan(seqb[s:e])
    print(f"  ...chr{c} done, cum H-x-H={n_hxh:,}")
ig=summarize("clean intergenic (repeat-masked)", n_hxh, thr, cac, cat, tot_his, tot_len, full)
IG=(n_hxh, thr, cac, cat)

print("\n=== VERDICT (Thr-middle / CAC / full-(CA)n) ===")
print(f"  neORF:      {100*neo[0]:.0f}% / {100*neo[1]:.0f}% / {100*neo[2]:.0f}%")
print(f"  intergenic: {100*ig[0]:.0f}% / {100*ig[1]:.0f}% / {100*ig[2]:.0f}%")
print(f"  (baseline: Thr is ~{100*1/20:.0f}% of residues by chance; CAC ~50% of His codons at neutrality)")

# ---------- Fisher exact: neORF vs intergenic (composition difference, NOT masking) ----------
# The (CA)n reading-frame signature is carried by neORF H-x-H far more than intergenic H-x-H,
# despite equal His DENSITY -> equal incidence. This is what inverts the Fig 1A order-fold.
print("\n=== neORF vs intergenic Fisher exact (does the (CA)n fingerprint differ?) ===")
orT, pT = stats.fisher_exact([[NEO[1], NEO[0]-NEO[1]], [IG[1], IG[0]-IG[1]]])
print(f"  Thr-middle (motif-level): neORF {NEO[1]}/{NEO[0]} ({100*NEO[1]/NEO[0]:.1f}%) vs "
      f"intergenic {IG[1]}/{IG[0]} ({100*IG[1]/IG[0]:.1f}%)  ->  OR={orT:.2f}  p={pT:.2e}")
orC, pC = stats.fisher_exact([[NEO[2], NEO[3]], [IG[2], IG[3]]])
print(f"  CAC/CAT (His-codon lvl) : neORF {NEO[2]}/{NEO[2]+NEO[3]} ({100*NEO[2]/(NEO[2]+NEO[3]):.1f}%) vs "
      f"intergenic {IG[2]}/{IG[2]+IG[3]} ({100*IG[2]/(IG[2]+IG[3]):.1f}%)  ->  OR={orC:.2f}  p={pC:.2e}")
print("  [His-codon test is pseudoreplicated (2 His per motif); Thr-middle is the clean motif-level test.]")
