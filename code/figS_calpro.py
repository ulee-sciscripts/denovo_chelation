#!/usr/bin/env python3
"""Supplementary figure -- apo-model calibration control for the ZMEG bis-His site.

apo ESMFold reconstructs calprotectin's (S100A8/S100A9) two known transition-metal
sites at HIGH confidence and with the donor atoms converged onto one ion, whereas the
de novo ZMEG His32/His36 pair sits at LOW pLDDT -- as does calprotectin's own flexible
C-terminal His103 (an internal negative control).  So the apo prediction is calibrated
where the answer is known, and is inconclusive for ZMEG.

  A  minimum coordinating donor (N/O) distance to a site partner  vs  per-residue pLDDT
  B  per-residue pLDDT for every site residue, with the 0.70 confidence line

All numbers recomputed from the coordinate files.  Austere greyscale."""
import os
import matplotlib.pyplot as plt
from hs_mpl import apply_style, panel_label, save, INK, GREY1, GREY2, GREY3, GREY4, GREYbg, WHITE

ROOT= os.environ.get("ZROOT", "../../zinc_analysis")   # analysis data tree (sibling of submission/)
OUT = os.environ.get("ZOUT", "figures")                # release figures/ dir
CALPRO = os.path.join(ROOT,"structures","calprotectin","calprotectin_heterodimer.pdb")
CAND3  = os.path.join(ROOT,"structures","cand3.pdb")
os.makedirs(OUT, exist_ok=True)
apply_style()

def parse_pdb(path):
    A=[]
    for ln in open(path):
        if ln.startswith(("ATOM","HETATM")):
            A.append(dict(name=ln[12:16].strip(),resn=ln[17:20].strip(),chain=ln[21].strip(),
                resi=int(ln[22:26]),x=float(ln[30:38]),y=float(ln[38:46]),z=float(ln[46:54]),b=float(ln[60:66])))
    return A
def donors(A,c,r):
    return [x for x in A if x['chain']==c and x['resi']==r and x['name'] in ('ND1','NE2','OD1','OD2')]
def plddt(A,c,r):
    v=[x['b'] for x in A if x['chain']==c and x['resi']==r and x['name']=='CA']; return v[0]
def mind(g1,g2):
    b=1e9
    for p in g1:
        for q in g2:
            b=min(b,((p['x']-q['x'])**2+(p['y']-q['y'])**2+(p['z']-q['z'])**2)**.5)
    return b
cp=parse_pdb(CALPRO)
zm=parse_pdb(CAND3)

# sites: (label, chain, resi, restype)
site1=[('A83','A',83,'His'),('A87','A',87,'His'),('B20','B',20,'His'),('B30','B',30,'Asp')]
site2=[('A17','A',17,'His'),('B91','B',91,'His'),('B95','B',95,'His')]
zmeg =[('H32','A',32,'His'),('H36','A',36,'His')]

def nearest_in(A,members,c,r):
    return min(mind(donors(A,c,r),donors(A,oc,orr)) for (_,oc,orr,_) in members if (oc,orr)!=(c,r))
# build point list: (label, group, dist, pLDDT, restype)
pts=[]
for lab,c,r,ty in site1: pts.append((lab,'s1',nearest_in(cp,site1,c,r),plddt(cp,c,r),ty))
for lab,c,r,ty in site2: pts.append((lab,'s2',nearest_in(cp,site2,c,r),plddt(cp,c,r),ty))
# His103: distance to centroid of the site-2 core donor atoms (its splay from the confident triad)
import numpy as np
core=[]
for _,c,r,_ in site2: core+=[(a['x'],a['y'],a['z']) for a in donors(cp,c,r)]
cen=np.mean(np.array(core),axis=0)
h103=min(np.linalg.norm(np.array([a['x'],a['y'],a['z']])-cen) for a in donors(cp,'B',103))
pts.append(('H103','h103',float(h103),plddt(cp,'B',103),'His'))
for lab,c,r,ty in zmeg: pts.append((lab,'zmeg',nearest_in(zm,zmeg,c,r),plddt(zm,c,r),ty))

# group styles (shape only; austere greyscale)
GST={'s1':dict(marker='o',mfc=INK,  mec=INK,  lab='Calprotectin Site 1 (His$_3$Asp) — known'),
     's2':dict(marker='s',mfc=GREY2,mec=INK,  lab='Calprotectin Site 2 (His triad) — known'),
     'h103':dict(marker='D',mfc=WHITE,mec=INK, lab='Calprotectin His103 — low-conf. control'),
     'zmeg':dict(marker='^',mfc=WHITE,mec=INK, lab='ZMEG His32/His36 — de novo')}

fig,(axA,axB)=plt.subplots(1,2,figsize=(11.6,4.9),gridspec_kw=dict(width_ratios=[1.05,1.0],wspace=0.28))
fig.subplots_adjust(left=0.075,right=0.985,top=0.86,bottom=0.14)

# ---------------- Panel A: distance vs pLDDT ----------------
CONF=0.70; COORD=4.5
axA.axvspan(0,CONF,color=GREYbg,zorder=0)                       # low-confidence zone
axA.axhspan(0,COORD,xmin=0,xmax=1,color='none')
axA.axhline(COORD,color=GREY3,lw=1.0,ls=(0,(4,3)),zorder=1)
axA.axvline(CONF,color=INK,lw=1.1,ls=(0,(5,3)),zorder=1)
for lab,g,d,p,ty in pts:
    st=GST[g]
    axA.scatter(p,d,s=95,marker=st['marker'],facecolor=st['mfc'],edgecolor=st['mec'],
                linewidth=1.2,zorder=4)
    if g in ('zmeg','h103'):                      # label only the low-confidence points
        off=(8,3) if lab!='H36' else (-8,3)
        ha='left' if lab!='H36' else 'right'
        axA.annotate(lab,(p,d),textcoords="offset points",xytext=off,ha=ha,fontsize=8.0,color=GREY1)
# zone labels split to opposite empty corners (clear of legend + point cluster)
axA.text(CONF-0.012,2.30,"apo prediction\ninconclusive  ←",ha="right",va="bottom",fontsize=8.4,color=GREY1,style="italic")
axA.text(CONF+0.015,8.55,"→ confident\n(pLDDT ≥ 0.70)",ha="left",va="top",fontsize=8.4,color=GREY1,style="italic")
axA.set_xlim(0.35,1.0); axA.set_ylim(2,9)
axA.set_xlabel("per-residue pLDDT (apo ESMFold)")
axA.set_ylabel("min. coordinating donor distance to\na site partner (Å)")
panel_label(axA,"A",x=-0.13,y=1.04)
# legend (shape key)
from matplotlib.lines import Line2D
handles=[Line2D([0],[0],marker=GST[g]['marker'],mfc=GST[g]['mfc'],mec=GST[g]['mec'],
                ls='none',ms=8,mew=1.2,label=GST[g]['lab']) for g in ('s1','s2','h103','zmeg')]
axA.legend(handles=handles,loc='upper left',fontsize=7.7,handlelength=1.4,
           handletextpad=0.5,borderpad=0.5,labelspacing=1.0)

# ---------------- Panel B: per-residue pLDDT ----------------
order=[('A83','s1'),('A87','s1'),('B20','s1'),('B30','s1'),
       ('A17','s2'),('B91','s2'),('B95','s2'),
       ('H103','h103'),('H32','zmeg'),('H36','zmeg')]
pd={(l):p for l,g,d,p,ty in pts}
xs=range(len(order)); vals=[pd[l] for l,_ in order]
cols={'s1':INK,'s2':GREY2,'h103':GREY4,'zmeg':GREY3}
bars=axB.bar(xs,vals,width=0.72,color=[cols[g] for _,g in order],edgecolor=INK,lw=0.9,zorder=3)
axB.axhline(CONF,color=INK,lw=1.1,ls=(0,(5,3)),zorder=4)
axB.text(len(order)-0.4,CONF+0.015,"confidence 0.70",ha="right",va="bottom",fontsize=8.2)
for x,v in zip(xs,vals): axB.text(x,v+0.012,f"{v:.2f}",ha="center",va="bottom",fontsize=7.7)
axB.set_xticks(list(xs)); axB.set_xticklabels([l for l,_ in order],rotation=45,ha='right',fontsize=8.6)
axB.set_ylim(0,1.0); axB.set_ylabel("per-residue pLDDT")
# group brackets under axis
def grp_bracket(i0,i1,txt):
    axB.annotate('',xy=(i0-0.35,-0.20),xytext=(i1+0.35,-0.20),xycoords=('data','axes fraction'),
                 arrowprops=dict(arrowstyle='-',color=GREY2,lw=1.0))
    axB.text((i0+i1)/2,-0.255,txt,ha='center',va='top',fontsize=7.6,color=GREY1,
             transform=axB.get_xaxis_transform())
grp_bracket(0,3,"Calprotectin\nSite 1 (known)")
grp_bracket(4,6,"Site 2 (known)")
grp_bracket(7,7,"His103\nctrl")
grp_bracket(8,9,"ZMEG\nde novo")
panel_label(axB,"B",x=-0.13,y=1.04)

save(fig,os.path.join(OUT,"figS_calpro_control"))
print("wrote figS_calpro_control.png/.pdf")
for l,g,d,p,ty in pts: print(f"  {l:5s} {g:5s} dist={d:.2f} pLDDT={p:.2f} {ty}")
