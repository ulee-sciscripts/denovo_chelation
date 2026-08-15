"""Shared austere greyscale style for the de novo zinc figures (matplotlib).
Monochrome/greyscale, clean sans-serif, minimal chartjunk. 300 dpi PNG + PDF.
Companion to submission/code/hs_style.R (which governs the R base-graphics figures)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# greyscale ramp used across all panels (dark -> light)
INK   = "#111111"   # near-black: hero / de novo / His
GREY1 = "#3f3f3f"
GREY2 = "#6e6e6e"
GREY3 = "#9a9a9a"
GREY4 = "#c4c4c4"
GREYbg= "#ececec"   # light fill / spacer residues
WHITE = "#ffffff"

def apply_style():
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 10.5,
        "axes.linewidth": 0.9,
        "axes.edgecolor": INK,
        "axes.titlesize": 11.5,
        "axes.labelsize": 11,
        "xtick.color": INK, "ytick.color": INK,
        "text.color": INK, "axes.labelcolor": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,   # editable text in PDF
    })

def panel_label(ax, letter, x=-0.02, y=1.06, fs=15):
    """Bold uppercase floating panel label."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", va="bottom", ha="right")

def save(fig, stem):
    fig.savefig(stem + ".png", bbox_inches="tight", facecolor="white")
    fig.savefig(stem + ".pdf", bbox_inches="tight", facecolor="white")
