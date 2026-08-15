#!/usr/bin/env python3
"""Compose Figure 5 -- ZMEG's candidate bis-histidine metal site (3 panels).

Panel A  ZMEG predicted fold (pLDDT) + His32-His36 site with modeled ion   [PyMOL]
Panel B  Metal3D max p(zinc): ZMEG vs calprotectin positive control        [matplotlib bar]
Panel C  Structural His-cage comparison: ZMEG His32/36 vs calprotectin site 2 [PyMOL]

Austere/monochrome base; accent RED = histidine, grey sphere = metal, pLDDT scale in A.
Writes 300-dpi PNG + PDF.

This script COMPOSES the ray-traced PyMOL panels (produced by fig5_render_pymol.py)
into the final figure; it needs the panel PNGs under $FIG5_SCR (default ./fig5_render),
not PyMOL itself.  Render the panels first (see fig5_render_pymol.py, needs PyMOL +
the ZMEG / calprotectin structures), then run this script.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

# unified font (matches Figs 3,4,S3): sans-serif, Arial/Helvetica/DejaVu Sans
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

SCR = os.environ.get("FIG5_SCR", "fig5_render")   # dir with the PyMOL-rendered panels (see fig5_render_pymol.py)
OUTDIR = os.environ.get("FIG5_OUTDIR", "figures") # release figures/ dir
os.makedirs(OUTDIR, exist_ok=True)
DPI = 300

HIS_RED = "#C0392B"
CTRL_GREY = "#6f6f6f"
BG_GREY = "#c9c9c9"

# ---- Metal3D values (ZMEG re-derived from JSON this session; see report) ----
ZMEG_MAXP = 0.283          # result_cand3_ZMEG_mb.json max_p = 0.2829
# background "other de novo ORFs" bars (grey); labels are internal only -- the
# figure renders these as an unnamed group, so only the values are drawn.
OTHER = [("other1", 0.053), ("other2", 0.091), ("other3", 0.112)]  # cand1/2/4 mb max_p
# calprotectin positive control (Metal3D). Independently re-run this session:
# metalbinding mode, weights metal_0.5A_v3_d0.2_16Abox.pth -> max p = 0.899,
# top site (1.14,12.54,9.47) = 2.46 A from the site-2 His-N centroid (<1 voxel).
CALPRO_MAXP = float(os.environ.get("CALPRO_MAXP", "0.899"))

MARKER_RGB255 = [(255, 0, 0), (0, 178, 0), (0, 0, 255),
                 (255, 0, 255), (0, 216, 216), (255, 153, 0)]


# ----------------------------- image helpers -----------------------------
def autocrop(path, pad_frac=0.02, bg=245):
    im = Image.open(path).convert("RGB")
    arr = np.asarray(im)
    mask = np.any(arr < bg, axis=2)
    if not mask.any():
        return im
    ys, xs = np.where(mask)
    y0, y1, x0, x1 = ys.min(), ys.max(), xs.min(), xs.max()
    ph = int((y1 - y0) * pad_frac) + 10
    pw = int((x1 - x0) * pad_frac) + 10
    y0 = max(0, y0 - ph); x0 = max(0, x0 - pw)
    y1 = min(arr.shape[0], y1 + ph); x1 = min(arr.shape[1], x1 + pw)
    return im.crop((x0, y0, x1, y1))


def crop_and_markers(prefix, res_list):
    """Return (cropped display image, {idx:(x,y)}) for callout leaders."""
    disp = Image.open(f"{SCR}/{prefix}.png").convert("RGB")
    arr = np.asarray(disp)
    mask = np.any(arr < 245, axis=2)
    ys, xs = np.where(mask)
    pad = 22
    x0 = max(0, xs.min() - pad); y0 = max(0, ys.min() - pad)
    x1 = min(arr.shape[1], xs.max() + pad); y1 = min(arr.shape[0], ys.max() + pad)
    disp_c = disp.crop((x0, y0, x1, y1))
    mk = np.asarray(Image.open(f"{SCR}/{prefix}_markers.png").convert("RGB")).astype(int)
    cen = {}
    for i in range(len(res_list)):
        tr, tg, tb = MARKER_RGB255[i]
        d = np.sqrt((mk[:, :, 0]-tr)**2 + (mk[:, :, 1]-tg)**2 + (mk[:, :, 2]-tb)**2)
        m = d < 80
        if not m.any():
            continue
        yy, xx = np.where(m)
        cen[i] = (xx.mean() - x0, yy.mean() - y0)
    return disp_c, cen


# ----------------------------- Panel A -----------------------------
def panelA_fold(ax):
    his = [27, 29, 32, 34, 36]
    disp_c, cen = crop_and_markers("zmeg_overview", [("A", h) for h in his])
    W, Himg = disp_c.size
    ax.imshow(disp_c, aspect="equal")
    ax.set_xlim(-W*0.02, W*1.02)
    ax.set_ylim(Himg*1.02, -Himg*0.22)
    ax.axis("off")
    ax.text(0.5, 1.005, "Predicted fold (ESMFold, apo), coloured by confidence",
            transform=ax.transAxes, fontsize=11.5, ha="center", va="bottom",
            fontweight="semibold")
    # Explicit per-His placement (declutter). His29/His36 keep their TOP-row slots;
    # His27 and His34 move to the opposite (lower) face, into the open dead space
    # below the ribbon; His32 moves right into the top slot His34 vacated. Targets are
    # read from the baked sidechain pixels (cen = CB markers; ring positions measured
    # from the render) so no leader crosses another residue's sidechain, and each
    # leader ends just short of (His32, His34) or right at (His27, His29, His36) its
    # own sidechain. His27 sidechain points DOWN (ring tip ~(835,348)); His34 points
    # DOWN (ring bottom ~y320); His32 points UP-left (ring x850-915, y40-105).
    ly = -Himg * 0.12        # top label row (His29, His32, His36)
    lyb = Himg * 0.795       # lower label row in the dead space (His27, His34)
    placements = {
        27: dict(xytext=(810, lyb),   xy=(826, 310)),   # below; leader reaches its own (down) sidechain tip
        29: dict(xytext=(W*0.30, ly), xy=cen[1]),        # unchanged (top)
        32: dict(xytext=(W*0.70, ly), xy=(930, 75)),     # moved right; leader stops short, just right of ring
        34: dict(xytext=(1025, lyb),  xy=(983, 297)),    # below; leader extended slightly closer to the (down) sidechain
        36: dict(xytext=(W*0.90, ly), xy=cen[4]),        # unchanged (top)
    }
    for h in his:
        p = placements[h]
        ax.annotate(f"His{h}", xy=p["xy"], xytext=p["xytext"],
                    fontsize=11.5, fontweight="bold", ha="center", va="center",
                    annotation_clip=False,
                    bbox=dict(boxstyle="round,pad=0.22", fc="white",
                              ec="#555555", lw=0.8),
                    arrowprops=dict(arrowstyle="-", color="#444444", lw=1.0,
                                    shrinkA=5, shrinkB=3))


def panelA_site(ax):
    # cleaned render: the baked His32/His36 labels were inpainted out (they overlapped
    # the ribbon); relabelled here off-structure with short leader lines.
    ax.imshow(autocrop(f"{SCR}/zmeg_site_ion_nolabel.png"))
    ax.axis("off")
    ax.text(0.5, 1.005, "His32-His36 site with modelled ion",
            transform=ax.transAxes, fontsize=11.5, ha="center", va="bottom",
            fontweight="semibold")
    lab_kw = dict(fontsize=11.5, fontweight="bold", ha="center", va="center",
                  annotation_clip=False,
                  bbox=dict(boxstyle="round,pad=0.22", fc="white",
                            ec="#555555", lw=0.8),
                  arrowprops=dict(arrowstyle="-", color="#444444", lw=1.1,
                                  shrinkA=4, shrinkB=3))
    ax.annotate("His36", xy=(700, 376), xytext=(320, 215), **lab_kw)
    # His32 label moved DOWN + RIGHT so its leader approaches from below the ring
    # (near-horizontal, from the right) and ends at the edge of the CB down-stick,
    # instead of the old upper-right position whose leader clipped the imidazole.
    ax.annotate("His32", xy=(1455, 460), xytext=(1760, 555), **lab_kw)


def plddt_legend(ax):
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 1.02, "Model confidence (pLDDT)", fontsize=9.6,
            fontweight="bold", va="top", ha="center")
    bands = [("#0053D6", "Very high", ">90"), ("#65CBF3", "Confident", "70-90"),
             ("#FFDB13", "Low", "50-70"), ("#FF7D45", "Very low", "<50")]
    w, gap, x0 = 0.205, 0.03, 0.06
    ytop, h = 0.42, 0.28
    for i, (c, lab, rng) in enumerate(bands):
        xx = x0 + i * (w + gap)
        ax.add_patch(Rectangle((xx, ytop), w, h, facecolor=c,
                     edgecolor="black", lw=0.6))
        ax.text(xx + w/2, ytop - 0.05, lab, fontsize=8.2, ha="center",
                va="top", fontweight="semibold")
        ax.text(xx + w/2, ytop - 0.42, rng, fontsize=7.6, ha="center",
                va="top", color="#333333")


# ----------------------------- Panel B -----------------------------
def panelB_bar(ax):
    # bar positions: calprotectin (control), ZMEG (candidate), then 3 background greys
    xc, xz = 0.0, 1.0
    xo = [2.35, 2.95, 3.55]
    # threshold shading (grey, austere)
    ax.axhspan(0.75, 1.0, color="#000000", alpha=0.05, zorder=0)
    ax.axhspan(0.5, 0.75, color="#000000", alpha=0.025, zorder=0)
    ax.axhline(0.75, ls="--", lw=1.2, color="#444444", zorder=2)
    ax.axhline(0.5, ls=":", lw=1.2, color="#777777", zorder=2)

    ax.bar(xc, CALPRO_MAXP, width=0.62, color=CTRL_GREY, edgecolor="black",
           lw=0.8, zorder=3)
    ax.bar(xz, ZMEG_MAXP, width=0.62, color=HIS_RED, edgecolor="black",
           lw=0.8, zorder=3)
    for (nm, v), x in zip(OTHER, xo):
        ax.bar(x, v, width=0.42, color=BG_GREY, edgecolor="#8a8a8a",
               lw=0.6, zorder=3)

    # value labels
    ax.text(xc, CALPRO_MAXP + 0.02, f"{CALPRO_MAXP:.2f}", ha="center",
            va="bottom", fontsize=11, fontweight="bold")
    ax.text(xz, ZMEG_MAXP + 0.02, f"{ZMEG_MAXP:.2f}", ha="center",
            va="bottom", fontsize=11, fontweight="bold", color=HIS_RED)
    for (nm, v), x in zip(OTHER, xo):
        ax.text(x, v + 0.02, f"{v:.2f}", ha="center", va="bottom",
                fontsize=9.5, color="#555555")

    # threshold labels
    ax.text(3.95, 0.755, "0.75  high-confidence", ha="right", va="bottom",
            fontsize=10, color="#444444")
    ax.text(3.95, 0.505, "0.50  site call", ha="right", va="bottom",
            fontsize=10, color="#777777")

    ax.set_ylim(0, 1.02)
    ax.set_xlim(-0.6, 4.0)
    ax.set_ylabel("Metal3D  max  p(zinc)", fontsize=10.5, fontweight="semibold")
    ax.set_xticks([xc, xz])
    ax.set_xticklabels(["Calprotectin", "ZMEG"], fontsize=9.6)
    ax.tick_params(axis="x", length=0)
    # bracket for background group
    ax.annotate("", xy=(xo[0]-0.28, -0.075), xytext=(xo[-1]+0.28, -0.075),
                xycoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="-", color="#8a8a8a", lw=1.0),
                annotation_clip=False)
    ax.text((xo[0]+xo[-1])/2, -0.11, "other de novo ORFs", ha="center",
            va="top", fontsize=10, color="#666666",
            transform=ax.get_xaxis_transform())
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ----------------------------- Panel C -----------------------------
def panelC(axL, axR):
    axL.imshow(autocrop(f"{SCR}/zmeg_cluster.png")); axL.axis("off")
    axL.text(0.5, 1.005, "ZMEG  His32 / His36", transform=axL.transAxes,
             fontsize=10.8, ha="center", va="bottom", fontweight="semibold")

    axR.imshow(autocrop(f"{SCR}/calpro_cluster.png")); axR.axis("off")
    axR.text(0.5, 1.005, "Calprotectin  site 2  (His$_6$)",
             transform=axR.transAxes, fontsize=10.8, ha="center", va="bottom",
             fontweight="semibold")


# NB: the on-figure honesty caption was removed (moved to the manuscript caption);
# see main() below. Its former helper/text lived here and has been dropped as dead code.
def main():
    # The on-figure honesty caption has been removed (-> manuscript caption); crop the
    # freed bottom band so no empty gap remains. Panel rects/labels are kept at their
    # original absolute inch sizes/positions via a fraction remap (H_OLD -> H_NEW).
    H_OLD, DY = 12.2, 1.45
    H_NEW = H_OLD - DY
    ry = lambda y: (y * H_OLD - DY) / H_NEW                 # remap a y fraction
    R  = lambda x, y, w, h: [x, (y * H_OLD - DY) / H_NEW, w, h * H_OLD / H_NEW]

    fig = plt.figure(figsize=(13.5, H_NEW))
    fig.patch.set_facecolor("white")
    fig.text(0.5, ry(0.982),
             "ZMEG: a candidate bis-histidine metal site in a de novo protein",
             fontsize=16, fontweight="bold", ha="center", va="top")

    # ---- Panel A (top) ----
    axAf = fig.add_axes(R(0.030, 0.585, 0.44, 0.345))
    panelA_fold(axAf)
    axAs = fig.add_axes(R(0.505, 0.585, 0.45, 0.345))
    panelA_site(axAs)
    axLeg = fig.add_axes(R(0.055, 0.525, 0.39, 0.045))
    plddt_legend(axLeg)
    fig.text(0.012, ry(0.945), "A", fontsize=24, fontweight="bold", va="top", ha="left")

    # ---- Panel B (bottom-left) ----
    axB = fig.add_axes(R(0.075, 0.185, 0.345, 0.285))
    panelB_bar(axB)
    fig.text(0.012, ry(0.485), "B", fontsize=24, fontweight="bold", va="top", ha="left")

    # ---- Panel C (bottom-right) ----
    axCL = fig.add_axes(R(0.475, 0.175, 0.245, 0.305))
    axCR = fig.add_axes(R(0.725, 0.175, 0.255, 0.305))
    panelC(axCL, axCR)
    fig.text(0.452, ry(0.485), "C", fontsize=24, fontweight="bold", va="top", ha="left")

    png = os.path.join(OUTDIR, "fig5.png")
    pdf = os.path.join(OUTDIR, "fig5.pdf")
    fig.savefig(png, dpi=DPI, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    print("wrote", png)
    print("wrote", pdf)


if __name__ == "__main__":
    main()
