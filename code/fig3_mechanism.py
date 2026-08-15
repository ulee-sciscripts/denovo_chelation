#!/usr/bin/env python3
"""Figure 3 -- the (CA)n microsatellite -> His-x-H mechanism, worked through the
ZMEG de novo flagship.  Non-duplicative with Fig 1 (H-x-H enrichment/incidence vs
comparators) and Fig 2 (His-codon CAC bias + Thr middle-residue data): this figure
is the reading-frame *mechanism* made visual, plus the poly-(H-x-H) multiplicity tail.

  A  (CA)n reading-frame codon logic: (+) strand CAC.ACA.CAC -> His-Thr-His (poly H-x-H,
     beta-strand register, metal-coordination-competent); (-) strand -> Cys-Val-Cys.
  B  ZMEG worked example: (CA)9 tract -> His32/34/36 triple; apo His32-His36 = 4.63 A.
  C  de novo H-x-H multiplicity across the neORF set (the poly-HxH tail); ZMEG marked.

Recomputes A/B/C from neORF.fa so the figure is self-contained.  Austere greyscale."""
import os
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from hs_mpl import apply_style, panel_label, save, INK, GREY1, GREY2, GREY3, GREY4, GREYbg, WHITE

ROOT = os.environ.get("ZROOT", "../../zinc_analysis")  # analysis data tree (sibling of submission/)
OUT  = os.environ.get("ZOUT", "figures")               # release figures/ dir
os.makedirs(OUT, exist_ok=True)
apply_style()

# ------------------------------------------------------------------ data
CODON = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def translate(nt):
    aa=[]
    for i in range(0,len(nt)-2,3):
        a=CODON.get(nt[i:i+3],'X')
        if a=='*': break
        aa.append(a)
    return ''.join(aa)
def read_fa(path):
    h=None; s=[]; d={}
    for ln in open(path):
        ln=ln.rstrip('\n')
        if ln.startswith('>'):
            if h: d[h]=''.join(s)
            h=ln[1:].split()[0]; s=[]
        else: s.append(ln)
    if h: d[h]=''.join(s)
    return d
nt = read_fa(os.path.join(ROOT,"neORF.fa"))
prot={}; seen=set(); uniq=[]
for h,s in nt.items():
    aa=translate(s.upper()); prot[h]=aa
    if aa not in seen: seen.add(aa); uniq.append(h)
def hpos(aa): return [i for i,c in enumerate(aa) if c=='H']
def nhxh(aa):
    hp=set(hpos(aa)); return sum(1 for i in hpos(aa) if i+2 in hp)
def is_poly(aa):
    hp=set(hpos(aa)); return any(i+2 in hp and i+4 in hp for i in hpos(aa))
N=len(uniq)
mult=Counter(nhxh(prot[h]) for h in uniq)
carriers=sum(v for k,v in mult.items() if k>=1)
npoly=sum(1 for h in uniq if is_poly(prot[h]))
ZID='UA_STRG.15519.1_1_193_357'
zaa=prot[ZID]
zhis=[i+1 for i in hpos(zaa)]                       # 1-based His positions
n_z=nhxh(zaa)

# ================================================================== figure
BAND="#f4f4f4"          # faint His-codon emphasis wash (kept within the grey ramp)

# ---- one His icon + one Zn icon, used identically everywhere (Panels A & B) ----
HIS_MS = 15.0                          # single His pentagon size (points)
ZN_MS  = 25.0                          # Zn circle: sized so "Zn2+" sits fully inside
def draw_his(ax, x, y, ms=HIS_MS, fill=INK, zorder=6):
    """The His icon: a vertex-DOWN black pentagon, identical size / orientation / fill
    in the sequence row, all four Panel-A cartoons, and the Panel-B lollipops."""
    ax.plot(x, y, marker=(5,0,180), ms=ms, mfc=fill, mec=INK, mew=1.0, zorder=zorder)
def draw_zn(ax, x, y, zorder=6):
    """Bound Zn(2+): filled grey circle large enough for the white label to fit inside."""
    ax.plot(x, y, marker='o', ms=ZN_MS, mfc=GREY2, mec=INK, mew=1.2, zorder=zorder)
    ax.text(x, y, "Zn$^{2+}$", ha="center", va="center", fontsize=8.0,
            fontweight="bold", color=WHITE, zorder=zorder+1)

def hx_cartoon(ax, xc, yb, n_between, name, coordinating, step=3.0):
    """One His-spacing cartoon.
    coordinating=True : the two His sit on the same strand face and converge on a
        bound Zn (favourable, e.g. i,i+2 and i,i+4).
    coordinating=False: the His fall on opposite faces (one above, one below the
        backbone), so they point away from a bound Zn and cannot both reach one
        ion (geometrically disfavoured, e.g. i,i+1 and i,i+3).  The Zn is drawn at
        the same spot; no bonds are drawn -- the point reads from geometry alone.
    His = flipped imidazole pentagon icon (single draw_his routine, identical
    everywhere); intervening residues = 'x' spacer beads; flanking (x) notation
    shows the chain continues either way."""
    nres = 2 + n_between
    half = (nres - 1) * step / 2.0
    xs   = [xc - half + i*step for i in range(nres)]
    flank= step*0.8
    his  = (0, nres-1)
    metal= (xc, yb-8.6)
    ax.plot([xs[0]-flank, xs[-1]+flank],[yb,yb],color=INK,lw=1.4,
            solid_capstyle="round",zorder=3)                       # backbone
    for i,xx in enumerate(xs):                                     # spacer residues
        if i not in his:
            ax.plot(xx,yb,marker='o',ms=8,mfc=GREY4,mec=INK,mew=0.9,zorder=5)
            ax.text(xx,yb+3.3,"x",ha="center",va="bottom",fontsize=10.5,
                    style="italic",color=GREY1,zorder=6)
    if coordinating:
        for i in his:                                             # coordination bonds
            ax.plot([xs[i],metal[0]],[yb,metal[1]],color=INK,lw=1.0,ls=(0,(1,1.4)),zorder=4)
        for i in his:                                             # His on-backbone, same face
            draw_his(ax, xs[i], yb)
            ax.text(xs[i],yb+3.3,"H",ha="center",va="bottom",fontsize=12,
                    fontweight="bold",color=INK,zorder=6)
        draw_zn(ax, *metal)
    else:
        off_up, off_dn = 5.5, 3.0        # up His clears the ion; down His kept high off it
        for i in his:                                            # His on OPPOSITE faces
            up = (i != 0)                                        # left His below, right His above
            hy = yb + (off_up if up else -off_dn)
            ax.plot([xs[i],xs[i]],[yb,hy],color=INK,lw=1.2,zorder=4)          # His stem to its face
            draw_his(ax, xs[i], hy)
            if up:
                ax.text(xs[i],hy+2.6,"H",ha="center",va="bottom",fontsize=12,
                        fontweight="bold",color=INK,zorder=6)
            else:
                ax.text(xs[i]-2.4,hy,"H",ha="right",va="center",fontsize=12,
                        fontweight="bold",color=INK,zorder=6)
        draw_zn(ax, *metal)
    ax.text(xs[0]-flank-0.7,yb,"…(x)–",ha="right",va="center",fontsize=10,color=GREY1,zorder=6)
    ax.text(xs[-1]+flank+0.7,yb,"–(x)…",ha="left",va="center",fontsize=10,color=GREY1,zorder=6)
    ax.text(xc,yb-13.5,name,ha="center",va="top",fontsize=12,fontweight="bold",color=INK,zorder=6)

fig=plt.figure(figsize=(13.0,8.8))
gs=fig.add_gridspec(2,2,height_ratios=[1.5,1.0],width_ratios=[1.3,1.0],
                    hspace=0.2,wspace=0.24,left=0.05,right=0.985,top=0.955,bottom=0.07)
axA=fig.add_subplot(gs[0,:]); axB=fig.add_subplot(gs[1,0]); axC=fig.add_subplot(gs[1,1])

# ---------------------------------------------------------- PANEL A : mechanism
axA.set_xlim(0,100); axA.set_ylim(26,92); axA.axis("off")
panel_label(axA,"A",x=0.0,y=1.0,fs=18)

plus="CA"*6                                   # 12 nt (+) strand, His-Thr reading frame
x0=2.0; pitch=3.7; cw=3.2; yd=64.0; hd=11.0
cellx = lambda k: x0+k*pitch
cellc = lambda k: x0+k*pitch+cw/2.0
codonc= lambda c: cellc(3*c+1)

# grey wash behind the two His codons -> conveys the codon grouping (no dashes needed);
# spans only the two close rows, marking the His codon in both the DNA and protein registers
for c in (0,2):
    xl=cellx(3*c)-0.3; xr=cellx(3*c+2)+cw+0.3
    axA.add_patch(Rectangle((xl,40),xr-xl,35,facecolor=BAND,edgecolor="none",zorder=0))

# (+) strand DNA ribbon
axA.text(x0-2.2,yd+hd/2,"(+)",ha="right",va="center",fontsize=12.5,fontweight="bold")
axA.text(cellc(0),yd+hd+2.4,"5'",ha="center",va="bottom",fontsize=12)
axA.text(cellc(11),yd+hd+2.4,"3'",ha="center",va="bottom",fontsize=12)
for k,ch in enumerate(plus):
    axA.add_patch(Rectangle((cellx(k),yd),cw,hd,facecolor=WHITE,edgecolor=INK,lw=0.8,zorder=2))
    axA.text(cellc(k),yd+hd/2,ch,ha="center",va="center",fontsize=12.5,zorder=3)

# translated residues: His=flipped imidazole pentagon (hero), Thr=spacer bead
# sits just below the DNA row -- a compact codon-over-residue block, vertically
# centred so its midpoint aligns with the middle of the 2x2 cartoon grid
yr=51.0
axA.plot([codonc(0)-2.8,codonc(3)+2.8],[yr,yr],color=INK,lw=1.4,
         solid_capstyle="round",zorder=2)
for c in range(4):
    aa=CODON[plus[c*3:c*3+3]]; xc=codonc(c)
    if aa=='H':
        draw_his(axA, xc, yr, zorder=5)
        axA.text(xc,yr-7.6,"His",ha="center",va="top",fontsize=14.5,fontweight="bold",color=INK)
    else:
        axA.plot(xc,yr,marker='o',ms=11,mfc=GREY4,mec=INK,mew=1.0,zorder=5)
        axA.text(xc,yr-7.6,"Thr",ha="center",va="top",fontsize=13.5,color=GREY1)
axA.text((codonc(1)+codonc(2))/2.0,yr-13.0,"poly (H-x-H)",ha="center",va="top",
         fontsize=11,style="italic",color=GREY1)

# ---- four His-spacing cartoons: coordinating (left) vs geometrically disfavoured (right)
cxL,cxR=61.0,87.0; yb_t,yb_b=70.0,44.0
axA.plot([74,74],[30,85],color=GREY4,lw=0.8,zorder=1)                 # group divider
axA.text(cxL,88.0,"coordinates",ha="center",va="bottom",fontsize=11.5,fontweight="bold",color=INK)
axA.text(cxL,84.3,"His same face",ha="center",va="bottom",fontsize=9.5,style="italic",color=GREY2)
axA.text(cxR,88.0,"disfavoured",ha="center",va="bottom",fontsize=11.5,fontweight="bold",color=GREY1)
axA.text(cxR,84.3,"His opposite faces",ha="center",va="bottom",fontsize=9.5,style="italic",color=GREY2)
hx_cartoon(axA,cxL,yb_t,1,"H–x–H",True)      # i, i+2   (x = Thr)
hx_cartoon(axA,cxL,yb_b,3,"H–x–x–x–H",True)  # i, i+4   (ZMEG His32·His36)
hx_cartoon(axA,cxR,yb_t,0,"H–H",False)       # i, i+1   adjacent
hx_cartoon(axA,cxR,yb_b,2,"H–x–x–H",False)   # i, i+3

# ---------------------------------------------------------- PANEL B : ZMEG example
axB.set_xlim(0,60); axB.set_ylim(3.0,8.1); axB.axis("off")
panel_label(axB,"B",x=0.0,y=1.04,fs=18)
xs_,xe_,barlen=13.0,58.0,45.0
xp=lambda p: xs_+(p-1)/54.0*barlen
by=3.9; bh=1.15
axB.add_patch(Rectangle((xs_,by),barlen,bh,facecolor=GREYbg,edgecolor=INK,lw=1.1,zorder=2))
t0,t1=xp(31),xp(37)
axB.add_patch(Rectangle((t0,by),t1-t0,bh,facecolor=GREY4,edgecolor=INK,lw=1.0,zorder=3))
axB.text((t0+t1)/2,by-0.30,"(CA)$_9$",ha="center",va="top",fontsize=11,style="italic",color=GREY1)
axB.text(xs_,by-0.30,"1",ha="center",va="top",fontsize=11,color=GREY1)
axB.text(xe_,by-0.30,"55 aa",ha="center",va="top",fontsize=11,color=GREY1)
triple=(32,34,36)
for p in zhis:                                                # His lollipops (pentagon heads)
    x=xp(p); tri=p in triple
    top=by+bh+(2.05 if tri else 1.15)
    axB.plot([x,x],[by+bh,top],color=INK,lw=(2.0 if tri else 1.1),
             solid_capstyle="round",zorder=4)
    if tri:
        draw_his(axB, x, top, zorder=5)
    else:
        draw_his(axB, x, top, ms=8.5, fill=WHITE, zorder=5)
axB.text(xp(34),by+bh+2.05+0.55,"His32 · His34 · His36",ha="center",va="bottom",
         fontsize=12.5,fontweight="bold")
axB.text(1.0,by+bh/2,"ZMEG",ha="left",va="center",fontsize=23,fontweight="bold",color=INK)

# ---------------------------------------------------------- PANEL C : multiplicity
panel_label(axC,"C",x=-0.15,y=1.05,fs=18)
axC.text(0.0,1.05,"de novo H-x-H multiplicity",transform=axC.transAxes,
         fontsize=14.5,fontweight="bold",va="bottom",ha="left")
labels=["1","2","3","4","5","≥6"]
vals=[mult[1],mult[2],mult[3],mult[4],mult[5],sum(mult[k] for k in mult if k>=6)]
xpos=list(range(len(vals))); ymax=max(vals)
axC.set_ylim(0,ymax*1.16)
cols=[GREY3]*len(vals); cols[2]=INK          # emphasise ZMEG's bin (n=3)
axC.bar(xpos,vals,width=0.7,color=cols,edgecolor=INK,lw=0.9,zorder=3)
for x,v in zip(xpos,vals):
    axC.text(x,v+ymax*0.018,str(v),ha="center",va="bottom",fontsize=11.5)
axC.set_xticks(xpos); axC.set_xticklabels(labels,fontsize=12.5)
axC.tick_params(axis='y',labelsize=12)
axC.set_xlabel("number of H-x-H motifs in the neORF",fontsize=13)
axC.set_ylabel("neORFs (carriers)",fontsize=13)
axC.annotate("ZMEG",xy=(2.36,vals[2]+ymax*0.008),xytext=(3.5,ymax*0.40),
             ha="center",va="center",fontsize=13.5,fontweight="bold",
             arrowprops=dict(arrowstyle="-",color=INK,lw=1.2,shrinkA=6,shrinkB=3))

save(fig,os.path.join(OUT,"fig3"))
print("wrote fig3.png/.pdf | N=%d carriers=%d poly=%d | ZMEG n_hxh=%d"%(N,carriers,npoly,n_z))
print("multiplicity:",dict(sorted(mult.items())))
