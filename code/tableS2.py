#!/usr/bin/env python3
"""Supplementary Table S2 -- TSV + austere booktabs-style PDF/PNG.

  S2  ZMEG (CA)n tract length -- within D. melanogaster strains (task3_tractlength.tsv)
      and across syntenic Drosophila orthologues (zmeg_synteny_summary.tsv):
      tract units + whether the H-x-H motif is preserved.

(Table S1 is the machine-readable per-replicate expression CSV of the 267
mappability-clean neORFs -- generated separately by tableS1_expression.py, written
to ../data/Dataset_S1_neORF_expression.csv.)"""
import os, csv
import matplotlib.pyplot as plt
from hs_mpl import apply_style, save, INK, GREY1, GREY2, GREY3, GREYbg, WHITE

ROOT=os.environ.get("ZROOT","../data/analysis")   # analysis data tree (sibling of submission/)
OUT =os.environ.get("ZOUT","figures")                # release figures/ dir
TRACT  =os.path.join(ROOT,"avenueA","tables","task3_tractlength.tsv")
SYNTENY=os.path.join(ROOT,"zmeg_history","tables","zmeg_synteny_summary.tsv")
os.makedirs(OUT,exist_ok=True)
apply_style()

# ============================================================ booktabs renderer
def render_table(path_stem, title, headers, rows, colx, colalign, note,
                 section_rows=None, winch=9.5, rowh=0.34, fs=9.2, hfs=9.4):
    """One vertical unit = one body row; y increases downward (inverted axis)."""
    nsec=len(section_rows) if section_rows else 0
    nbody=len(rows)+nsec
    nnote=len(note) if note else 0
    HEADBLK=1.7                                        # units between top rule and header rule (fits 2-line headers)
    top_rule=1.15
    header_rule=top_rule+HEADBLK
    body_last=header_rule+0.6+(nbody-1)
    bottom_rule=body_last+0.6
    note_end=bottom_rule+0.55+nnote*0.9
    total=note_end+0.4
    fig_h=rowh*total
    fig=plt.figure(figsize=(winch,fig_h)); ax=fig.add_axes([0.015,0.0,0.97,1.0]); ax.axis("off")
    ax.set_xlim(0,1); ax.set_ylim(0,total); ax.invert_yaxis()
    ax.text(0.0,0.15,title,fontsize=11.5,fontweight="bold",va="top",ha="left")
    ax.plot([0,1],[top_rule,top_rule],color=INK,lw=1.6)
    for h,x,al in zip(headers,colx,colalign):
        ax.text(x,top_rule+HEADBLK*0.52,h,fontsize=hfs,fontweight="bold",va="center",ha=al)
    ax.plot([0,1],[header_rule,header_rule],color=INK,lw=1.0)
    yr=header_rule+0.6
    for idx,cells in enumerate(rows):
        if section_rows and idx in section_rows:
            ax.text(0.0,yr,section_rows[idx],fontsize=fs+0.3,fontweight="bold",fontstyle="italic",
                    va="center",ha="left",color=GREY1)
            ax.plot([0,1],[yr+0.5,yr+0.5],color=GREY3,lw=0.7)
            yr+=1.0
        for c,x,al in zip(cells,colx,colalign):
            ax.text(x,yr,c,fontsize=fs,va="center",ha=al,color=INK)
        yr+=1.0
    ax.plot([0,1],[bottom_rule,bottom_rule],color=INK,lw=1.6)
    yn=bottom_rule+0.7
    for nline in (note or []):
        ax.text(0.0,yn,nline,fontsize=7.7,va="top",ha="left",color=GREY1); yn+=0.9
    save(fig,path_stem); plt.close(fig)

# ============================================================ Table S2 : (CA)n tract length
# ---- within D. melanogaster (ZMEG (CA)n microsat) ----
mel_lab={"dm6":"D. melanogaster dm6 (ISO1, ref.)","ORw1118":"OreR / w1118 (lab)",
    "GC_UA_denovo":"Ukraine (de novo transcript)","gDNA_DK":"Denmark (DK)","gDNA_ES":"Spain (ES)",
    "gDNA_FI":"Finland (FI)","gDNA_SE":"Sweden (SE)","gDNA_TR":"Turkey (TR)","gDNA_UA":"Ukraine (UA)",
    "gDNA_ZI":"Zambia (ZI)"}
mel_order=["dm6","ORw1118","gDNA_DK","gDNA_ES","gDNA_FI","gDNA_SE","gDNA_TR","gDNA_UA","gDNA_ZI"]
VARMAP={
 "reference":"reference assembly",
 "identical_to_dm6(9 units)":"identical to reference",
 "identical":"identical to reference",
 "del6@107(CA-microsat/inframe)":"6-bp in-frame deletion (tract shortened)",
 "del2@111(CA-microsat/FS);del1@119(unique/FS)":"2- & 1-bp indels (frameshift)",
 "del2@111(CA-microsat/FS)":"2-bp deletion (frameshift)"}
def clean_var(v): return VARMAP.get(v.strip(), v.strip())
def preserved_mark(cleaned):
    return "Y‡" if "frameshift" in cleaned else "Y"
mel={}
with open(TRACT) as f:
    for r in csv.DictReader((l for l in f if not l.startswith("#")),delimiter="\t"):
        if r["symbol"]=="ZMEG" and r["metric_type"]=="(CA)n_microsat":
            mel[r["reference"]]=r
# ---- across Drosophila (synteny) ----
sp_rows=[]
with open(SYNTENY) as f:
    for r in csv.DictReader(f,delimiter="\t"): sp_rows.append(r)
def sp_note(r):
    q=r["ZMEG_region_aln_quality"].split("(")[0].strip(); extra=""
    pm=r["poly-His_motif_region"]
    if "collapsed" in pm: extra="; tract collapsed"
    elif "truncated" in pm: extra="; ORF truncated pre-tract"
    mya=str(r["approx_Mya_from_mel"]).lstrip("~").strip()
    mya = "0" if mya=="0" else mya
    return f"~{mya} Mya; {q}{extra}"

# build combined rows: [strain/species, (CA)n units, motif preserved, note]
rowsS2=[]; sect={}
# section A
sect[0]="Within D. melanogaster  (population strains; ZMEG (CA)n microsatellite)"
for key in mel_order:
    r=mel[key]; cv=clean_var(r["notable_variants"])
    rowsS2.append([mel_lab[key], r["CAn_units"], preserved_mark(cv), cv])
# section B
idxB=len(rowsS2)
sect[idxB]="Across Drosophila  (syntenic orthologues; His-x-H triple)"
for r in sp_rows:
    sp=r["species"].split("(")[0].strip()
    units=r["CAn_units(raw,single-asm)"]
    mp=r["HxHxH_triple"]
    rowsS2.append([sp, units, mp, sp_note(r)])
# TSV
with open(os.path.join(OUT,"tableS2.tsv"),"w") as f:
    f.write("group\tstrain_or_species\tCAn_tract_units\tHxH_motif_preserved\tnote\n")
    grp_lbl="within_D.melanogaster_strain"
    for i,cells in enumerate(rowsS2):
        if i==idxB: grp_lbl="across_Drosophila_species"
        f.write(grp_lbl+"\t"+"\t".join(cells)+"\n")
# formatted
headersS2=["Strain / species","(CA)$_n$\nunits","H-x-H\npreserved","Note"]
colxS2=[0.0,0.42,0.55,0.70]
alignS2=["left","center","center","left"]
render_table(os.path.join(OUT,"tableS2"),
    "Table S2.  ZMEG (CA)$_n$ tract length and H-x-H motif preservation, within and across species.",
    headersS2,rowsS2,colxS2,alignS2,section_rows=sect,
    note=["(CA)$_n$ units = number of CA dinucleotide repeats at the syntenic ZMEG locus "
          "(single-assembly counts; short tracts are indel-assembly-sensitive).  ‡ = tract present but frameshifted "
          "(H-x-H encoded in DNA, out of the ZMEG reading frame).",
          "Within-species: every population strain retains the (CA)$_n$ tract (6–10 units).  Across species: the "
          "(CA)$_n$→His-x-H triple is conserved throughout the melanogaster subgroup (5–9 units) and lost deeper "
          "(ananassae tract collapsed; pseudoobscura / virilis ORF truncated before the tract)."],
    winch=10.2,rowh=0.335,fs=8.8,hfs=9.2)
print("wrote tableS2.tsv/.png/.pdf")
