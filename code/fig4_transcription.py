#!/usr/bin/env python3
"""Manuscript Figure 4 -- unified, single-script, austere.

  A  Expression density in S2R+ (de novo neORFs vs annotated genes; STAR/featureCounts)
  B  Number of neORFs transcribed (counts at TPM thresholds; all / H-x-H / poly-H-x-H)
  C  De novo origin & activation trajectory of ZMEG across Drosophila
       (cladogram, all 9 tips; polarised events as asterisks (*) ON the preceding
        branches; per-tip character-state column: (CA)n units + ORF length, colour-coded)
  D  ZMEG ORF alignment (ATG -> first stop); His boxed; (CA)n-derived tract shaded

Layout: A,B across the top; C (tree + character states + branch asterisks) on the
LEFT and D (amino-acid alignment) on the RIGHT of the lower two-thirds.

Transcription panels (A,B) are recomputed here directly from star_counts/counts.txt
+ neORF_clean_mappable.tsv + neORF_clean_dm6.tsv (STAR + featureCounts unique reads,
length-normalized to TPM).  Trajectory panels (C,D) use the authoritative tables
zmeg_history/tables/{zmeg_synteny_summary,zmeg_orf_raw}.tsv.

Data root defaults to ../data/analysis (override: argv[1] or run from submission/code/).
"""
import sys, csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle

# ----------------------------------------------------------------------------- style
# Unified font across Figs 3,4,5,S3: sans-serif, Arial/Helvetica preferred (to match
# the R "sans" figures); DejaVu Sans is the fallback actually rendered where Arial is absent.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.default": "regular",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.9,
})
# austere greys (grey45/grey10/grey25/grey55/grey82)
G_ANN, G_DEN = "#737373", "#1a1a1a"
G25, G55, G82 = "#404040", "#8c8c8c", "#d1d1d1"
# meaningful state code (the (CA)n / ORF-quality gradient) -- the ONLY hues in C/D
COL = {"full": "#08519c", "poly": "#6baed6", "lost": "#e6550d", "trunc": "#bdbdbd"}
HIS_FILL, HIS_EDGE, HIS_TXT = "#d5e5f2", "#08519c", "#08519c"   # His accent (Panel D)
STOP = "#b2182b"
EVT  = "#111111"   # polarised-event asterisk colour (Panel C branches)

ROOT = sys.argv[1] if len(sys.argv) > 1 else "../data/analysis"  # analysis data tree (sibling of submission/)
OUT  = sys.argv[2] if len(sys.argv) > 2 else "figures"             # release figures/ dir
DPI  = int(sys.argv[3]) if len(sys.argv) > 3 else 300

def plab(ax, s, x=-0.085, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=18, fontweight="bold",
            ha="left", va="top")

# ============================================================ transcription data (A,B)
d = pd.read_csv(f"{ROOT}/star_counts/counts.txt", sep="\t", comment="#")
L = d["Length"].values.astype(float)
cnt = d.iloc[:, 6:].values.astype(float)
rpk = cnt / (L[:, None] / 1000.0)
tpm = rpk / rpk.sum(axis=0, keepdims=True) * 1e6
mt = tpm.mean(axis=1)
gid = d["Geneid"].values
neo = pd.read_csv(f"{ROOT}/neORF_clean_mappable.tsv", sep="\t")
allneo = pd.read_csv(f"{ROOT}/neORF_clean_dm6.tsv", sep="\t")
isneo = np.isin(gid, neo["name"].values)
den = mt[isneo]
ann = mt[~np.isin(gid, allneo["name"].values)]
mp = dict(zip(neo["name"], zip(neo["has_hxh"], neo["poly_hxh"])))
hx = np.array([mp[g][0] for g in gid[isneo]])
pp = np.array([mp[g][1] for g in gid[isneo]])
denH, denP = den[hx == 1], den[pp == 1]
THS = [0.1, 0.5, 1, 5]

def r_density(x, frm=-1.0, npts=512, cut=3):
    """Reproduce R stats::density() default (nrd0 bandwidth, gaussian) with from=frm."""
    from scipy.stats import gaussian_kde
    x = np.asarray(x, float); n = len(x)
    sd = np.std(x, ddof=1)
    q75, q25 = np.percentile(x, 75), np.percentile(x, 25)
    lo = min(sd, (q75 - q25) / 1.349)
    if lo == 0: lo = sd if sd > 0 else (abs(x[0]) if x[0] != 0 else 1.0)
    bw = 0.9 * lo * n ** (-0.2)
    kde = gaussian_kde(x, bw_method=bw / sd)
    grid = np.linspace(frm, x.max() + cut * bw, npts)
    return grid, kde(grid)

# ============================================================ trajectory data (C,D)
orf = {}
for r in csv.DictReader(open(f"{ROOT}/zmeg_history/tables/zmeg_orf_raw.tsv"), delimiter="\t"):
    orf[r["species"]] = r
def O(sp, k): return orf[sp][k]

# tips top->bottom ; ladderised mel-group topology
# state: full / poly / lost / trunc
TAX = [  # (short, italic label, orf_raw key, state, status text)
    ("mel", "D. melanogaster", "D. melanogaster (dm6)", "full",  "full ORF + His-triple"),
    ("sim", "D. simulans",     "D. simulans",           "poly",  "poly-His"),
    ("sec", "D. sechellia",    "D. sechellia",          "poly",  "poly-His"),
    ("mau", "D. mauritiana",   "D. mauritiana",         "full",  "full ORF + His-triple"),
    ("yak", "D. yakuba",       "D. yakuba",             "poly",  "poly-His"),
    ("ere", "D. erecta",       "D. erecta",             "poly",  "poly-His"),
    ("ana", "D. ananassae",    "D. ananassae",          "lost",  "ORF open, His-triple lost"),
    ("pse", "D. pseudoobscura","D. pseudoobscura",      "trunc", "ORF truncated (pre-tract)"),
    ("vir", "D. virilis",      "D. virilis",            "trunc", "ORF truncated (pre-tract)"),
]
YT = {t[0]: 9 - i for i, t in enumerate(TAX)}  # mel=9 ... vir=1
TIPX = 6.0
NODES = {  # (x, y)  -- schematic cladogram (branch lengths not to scale)
    "secmau": (5.4, 6.5), "simclade": (4.7, 7.25), "melcomplex": (4.0, 8.125),
    "yakcomplex": (4.0, 4.5), "subgroup": (3.0, 6.3125), "melgroup": (2.0, 4.65625),
    "soph": (1.0, 3.328125), "root": (0.0, 2.1640625),
}
def P(k):
    if k in NODES: return NODES[k]
    return (TIPX, YT[k])
CHILDREN = {
    "secmau": ["sec", "mau"], "simclade": ["sim", "secmau"],
    "melcomplex": ["mel", "simclade"], "yakcomplex": ["yak", "ere"],
    "subgroup": ["melcomplex", "yakcomplex"], "melgroup": ["subgroup", "ana"],
    "soph": ["melgroup", "pse"], "root": ["soph", "vir"],
}

# ============================================================ figure scaffold
#  A , B  across the top;  C (tree) left  |  D (alignment) right  below.
fig = plt.figure(figsize=(13.0, 11.0))

axA = fig.add_axes([0.070, 0.740, 0.505, 0.205])
axB = fig.add_axes([0.680, 0.740, 0.290, 0.205])
axC = fig.add_axes([0.010, 0.105, 0.470, 0.575]); axC.axis("off")
axD = fig.add_axes([0.500, 0.105, 0.490, 0.575]); axD.axis("off")

# ------------------------------------------------------------ PANEL A: density
aD_x, aD_y = r_density(np.log10(ann[ann >= 0.1]))
dD_x, dD_y = r_density(np.log10(den[den >= 0.1]))
ymax = max(aD_y.max(), dD_y.max()) * 1.05
fUN, nUN = (den < 0.1).mean(), int((den < 0.1).sum())
fUA, nUA = (ann < 0.1).mean(), int((ann < 0.1).sum())

axA.plot(aD_x, aD_y, color=G_ANN, lw=2.3, zorder=3)
axA.fill_between(aD_x, aD_y, color=G_ANN, alpha=0.16, zorder=1)
axA.plot(dD_x, dD_y, color=G_DEN, lw=2.3, zorder=4)
axA.fill_between(dD_x, dD_y, color=G_DEN, alpha=0.26, zorder=2)
# undetected strip (left of the dotted separator).  The two bars are EVENLY
# distributed in the region [left axis (-2.75) -> dotted line (-1)] so the three
# horizontal gaps are equal: axis->bar1, bar1->bar2, bar2->dotted line.
# strip width = 1.75, bar width w = 0.35 -> each gap = (1.75 - 2*0.35)/3 = 0.35.
# de novo bar -2.40..-2.05, annotated bar -1.70..-1.35.  Labels centred on bars.
axA.add_patch(Rectangle((-2.40, 0), 0.35, fUN * ymax, facecolor=G_DEN, edgecolor="black", lw=0.8))
axA.add_patch(Rectangle((-1.70, 0), 0.35, fUA * ymax, facecolor=G_ANN, edgecolor="black", lw=0.8))
axA.text(-2.225, fUN * ymax + 0.11 * ymax, f"{100*nUN/len(den):.0f}%", ha="center", fontsize=10, fontweight="bold")
axA.text(-2.225, fUN * ymax + 0.045 * ymax, f"{nUN}/{len(den)}", ha="center", fontsize=8, fontweight="bold")
axA.text(-1.525, fUA * ymax + 0.11 * ymax, f"{100*nUA/len(ann):.0f}%", ha="center", fontsize=9)
axA.text(-1.525, fUA * ymax + 0.045 * ymax, f"{nUA}/{len(ann)}", ha="center", fontsize=7.5)
axA.text(-1.875, ymax * 0.99, "<0.1 TPM", ha="center", fontsize=8.5)
axA.axvline(-1, ls=(0, (2, 2)), color="grey", lw=1.2)
axA.set_xlim(-2.75, 4.2); axA.set_ylim(0, ymax)
axA.set_xticks([-1, 0, 1, 2, 3, 4])
axA.set_xlabel(r"$\log_{10}$ TPM", fontsize=11.5)
axA.set_ylabel("density", fontsize=11.5)
axA.set_title("Expression in S2R+", fontsize=12.5, fontweight="bold", pad=8)
axA.tick_params(labelsize=10)
axA.spines[["top", "right"]].set_visible(False)
axA.plot([], [], color=G_ANN, lw=2.6, label="annotated genes")
axA.plot([], [], color=G_DEN, lw=2.6, label="de novo ORFs (non-overlapping)")
axA.legend(loc="upper right", frameon=False, fontsize=9.2, handlelength=1.5)
plab(axA, "A", x=-0.11)

# ------------------------------------------------------------ PANEL B: counts
def cnt_ge(x, t): return int((x >= t).sum())
mat = np.array([[cnt_ge(den, t) for t in THS],
                [cnt_ge(denH, t) for t in THS],
                [cnt_ge(denP, t) for t in THS]], float)
greys = [G25, G55, G82]
grp = np.arange(len(THS)); bw = 0.26
for i in range(3):
    xs = grp + (i - 1) * bw
    axB.bar(xs, mat[i], width=bw, color=greys[i], edgecolor="black", lw=0.8, zorder=3)
    for x, v in zip(xs, mat[i]):
        axB.text(x, v + mat.max() * 0.025, f"{int(v)}", ha="center", va="bottom", fontsize=8.6)
axB.set_xticks(grp); axB.set_xticklabels([f"≥{t:g}" for t in THS], fontsize=10)
axB.set_ylim(0, mat.max() * 1.22)
axB.set_ylabel("number of neORFs detected", fontsize=11)
axB.set_xlabel("mean TPM threshold", fontsize=11.5)
axB.tick_params(labelsize=10)
axB.spines[["top", "right"]].set_visible(False)
hlab = [f"all non-overlapping neORFs (n={len(den)})", f"H-x-H (n={len(denH)})", f"poly-H-x-H (n={len(denP)})"]
handles = [Rectangle((0, 0), 1, 1, facecolor=greys[i], edgecolor="black", lw=0.8) for i in range(3)]
axB.legend(handles, hlab, loc="upper right", frameon=False, fontsize=8.6, handlelength=1.2, handleheight=1.2)
plab(axB, "B", x=-0.27, y=1.10)

# ------------------------------------------------------------ PANEL C: cladogram + states
# tightened y-range: dead top/bottom margins removed so C/D content hugs A/B above.
CY_LO, CY_HI = 0.4, 10.0
axC.set_xlim(-0.75, 13.2); axC.set_ylim(CY_LO, CY_HI)
def seg(x1, y1, x2, y2, lw=2.0, c="#2b2b2b"):
    axC.plot([x1, x2], [y1, y2], color=c, lw=lw, zorder=2, solid_capstyle="round")
for node, kids in CHILDREN.items():
    nx, ny = P(node)
    (c1x, c1y), (c2x, c2y) = P(kids[0]), P(kids[1])
    seg(nx, c1y, nx, c2y)                 # vertical connector
    seg(nx, c1y, c1x, c1y)                # horiz to child 1
    seg(nx, c2y, c2x, c2y)                # horiz to child 2
seg(-0.60, P("root")[1], 0.0, P("root")[1])   # root stub

# polarised events -- asterisks placed ON the branch PRECEDING each event's node
# (character-state changes are inferred along branches, not at speciation nodes).
AST_DY = -0.16   # nudge so the "*" ink sits centred on the branch line
EVENT_AST = [
    (-0.30, P("root")[1]),        # (1) ancestral: root stub
    ( 1.50, P("melgroup")[1]),    # (2) mel-group stem  (branch soph->melgroup)
    ( 2.50, P("subgroup")[1]),    # (3) mel-subgroup stem (branch melgroup->subgroup)
    ( 3.50, P("melcomplex")[1]),  # (4) mel-complex stem (branch subgroup->melcomplex)
    ( 4.30, YT["ana"]),           # (5) ananassae loss: terminal ananassae branch
]
for ex, ey in EVENT_AST:
    axC.text(ex, ey + AST_DY, "*", color=EVT, fontsize=25, fontweight="bold",
             ha="center", va="center", zorder=7,
             path_effects=[pe.withStroke(linewidth=3.2, foreground="white")])

# tips: name + state colour square + compact character-state text.
# box+label column nudged right (SQX 9.35->9.90) so long species names (e.g.
# D. pseudoobscura) clear the colour box with air; xlim widened to hold the text.
NAMEX, SQX, TXTX = 6.25, 9.90, 10.45
for short, ital, key, state, status in TAX:
    y = YT[short]
    seg(P(short)[0] - 0.02, y, TIPX + 0.18, y, lw=1.6)   # tidy tip extension
    axC.text(NAMEX, y, ital, fontsize=10, style="italic", va="center", fontweight="bold")
    axC.add_patch(Rectangle((SQX, y - 0.24), 0.30, 0.48, facecolor=COL[state], edgecolor="black", lw=0.5))
    can = O(key, "CAn_units_raw"); aa = O(key, "ORF_aa_to_stop")
    axC.text(TXTX, y, f"(CA)$_{{{can}}}$ · {aa} aa", fontsize=12.5, va="center")

# (state-colour legend is drawn figure-level, full-width, at the very bottom -- see below)
plab(axC, "C", x=-0.01, y=1.02)

# ------------------------------------------------------------ PANEL D: ORF alignment (right side)
PROT = [(t[1].replace("D. ", "D."), O(t[2], "translation") + "*") for t in TAX]
maxlen = max(len(s) for _, s in PROT)          # 59 incl. trailing stop
axD.set_xlim(0, 1); axD.set_ylim(0, 1)
# per-row species labels removed (rows align to Panel C's tree tips), so the left
# margin can shrink -- this widens the per-residue advance and lets the AA text grow.
x0 = 0.015                                        # thin left margin (no species labels)
dx = (0.997 - x0) / maxlen                        # per-residue advance
# rows are aligned to Panel C's tree tips (axC and axD share the same vertical
# rectangle): tip y=9..1 in axC ylim -> the same fractional heights here.
ytop = (9 - CY_LO) / (CY_HI - CY_LO)              # D. melanogaster row == its tree tip
dy   = 1.0 / (CY_HI - CY_LO)                      # one tip step
TRACT0, TRACT1 = 25, 42                           # (CA)n-derived poly-His tract cols
RES_FS = 9.0                                       # enlarged AA text
BOX_H = 0.30 * dy                                 # His-box half-heights kept modest
for r, (name, seq) in enumerate(PROT):
    y = ytop - r * dy
    for i, ch in enumerate(seq):
        x = x0 + (i + 0.5) * dx
        if TRACT0 <= i < TRACT1:
            axD.add_patch(Rectangle((x - dx / 2, y - dy * 0.44), dx, dy * 0.88,
                          facecolor="#efefef", edgecolor="none", zorder=0))
        if ch == "H":
            axD.add_patch(Rectangle((x - dx * 0.45, y - BOX_H), dx * 0.90, 2 * BOX_H,
                          facecolor=HIS_FILL, edgecolor=HIS_EDGE, lw=0.8, zorder=1))
            axD.text(x, y, "H", fontsize=RES_FS, color=HIS_TXT,
                     fontweight="bold", ha="center", va="center", zorder=2)
        elif ch == "*":
            axD.text(x, y, "*", fontsize=RES_FS + 2, color=STOP,
                     fontweight="bold", ha="center", va="center", zorder=2)
        else:
            axD.text(x, y, ch, fontsize=RES_FS, color="#222",
                     ha="center", va="center", zorder=2)
# tract bracket label
xb0, xb1 = x0 + TRACT0 * dx, x0 + TRACT1 * dx
ybk = ytop + dy * 0.60
axD.plot([xb0, xb1], [ybk, ybk], color="#666", lw=1.0)
axD.text((xb0 + xb1) / 2, ybk + 0.008, "(CA)$_n$-derived poly-His tract",
         fontsize=11.0, color="#444", ha="center", va="bottom")
plab(axD, "D", x=-0.005, y=1.02)

# ---------------------------------------------------------- figure-level state legend
# One evenly-spaced, full-width row across A-D at the very bottom.  Each swatch+label
# unit is centred at 1/8, 3/8, 5/8, 7/8 of the figure width; text widths are measured
# from the renderer so the centring is exact and nothing overlaps or runs off-edge.
COLABS = [("full",  "full ORF + His-triple (mature)"),
          ("poly",  "poly-His, ORF frame labile"),
          ("lost",  "ORF open, His-triple lost"),
          ("trunc", "ORF truncated (no protein)")]
axL = fig.add_axes([0.0, 0.045, 1.0, 0.055]); axL.axis("off")
axL.set_xlim(0, 1); axL.set_ylim(0, 1)
fig.canvas.draw(); _rend = fig.canvas.get_renderer()
LEG_FS = 11.5
sq_w = 0.013
sq_h = sq_w * 13.0 / (0.055 * 11.0)               # aspect-corrected -> visually square
gap  = 0.009
for (k, lab), cx in zip(COLABS, [0.125, 0.375, 0.625, 0.875]):
    _t = axL.text(0, 0, lab, fontsize=LEG_FS)
    tw = _t.get_window_extent(renderer=_rend).width / fig.bbox.width
    _t.remove()
    x_sq = cx - (sq_w + gap + tw) / 2.0           # centre the whole swatch+label unit
    axL.add_patch(Rectangle((x_sq, 0.5 - sq_h / 2), sq_w, sq_h,
                  facecolor=COL[k], edgecolor="black", lw=0.6))
    axL.text(x_sq + sq_w + gap, 0.5, lab, fontsize=LEG_FS, va="center", ha="left")

fig.savefig(f"{OUT}/fig4.png", dpi=DPI)
fig.savefig(f"{OUT}/fig4.pdf")
print(f"wrote {OUT}/fig4.png + fig4.pdf  (dpi={DPI})")
