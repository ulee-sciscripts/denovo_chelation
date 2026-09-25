#!/usr/bin/env Rscript
## Fig 2 -- the (CA)n reading-frame mechanism, stripped to the minimum.  Two
## panels, side by side (swapped order):
##   A = His-codon usage: H-x-H His is CAC-biased relative to non-H-x-H His.
##   B = middle residue of de novo H-x-H motifs (Thr enriched; top middle codon ACA).
## Base plot() + points only; series told apart by pch (filled vs open); ONE
## comparison marker per panel (A: CAC at two x -> hs_bracket; B: Thr mid-vs-bg at
## one x -> hs_topstar).  No colour / fill / grid / bar.
source("hs_style.R")                # run from the code/ directory (see run_all.sh)
DATA <- file.path("..", "data", "hisspacing")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

mid   <- read.csv(file.path(DATA, "fig2_middle.csv"), stringsAsFactors = FALSE)
codon <- read.csv(file.path(DATA, "fig2_codon.csv"),  stringsAsFactors = FALSE)
st    <- read.csv(file.path(DATA, "fig2_stats.csv"),  stringsAsFactors = FALSE)
getv  <- function(k) st$value[st$key == k]
star_thr <- getv("star_thr")   # ***  (Thr middle-residue enrichment vs background)
star_cac <- getv("star_cac")   # ***  (H-x-H His CAC vs all His CAC)

hs_pdf(file.path(FIGS, "fig2.pdf"), width = 12.0, height = 5.9)
par(oma = c(1.8, 0, 0, 0))
layout(matrix(c(1, 2), nrow = 1), widths = c(1.0, 1.3))   # A codon (narrow), B middle (wide)

## ============================ Panel A ============================
## His-codon usage; x = 1 H-x-H His, x = 2 all His.  Filled = CAC, open = CAT.
hxh <- codon[codon$group == "H-x-H His", ]
nhx <- codon[codon$group == "non-H-x-H His", ]
cac  <- c(hxh$CAC_pct, nhx$CAC_pct)   # 56.2, 47.2
cat_ <- c(hxh$CAT_pct, nhx$CAT_pct)   # 43.8, 52.8

par(mar = c(4.6, 5.2, 3.0, 1.2))
plot(NA, xlim = c(0.5, 2.5), ylim = c(36, 66), xaxt = "n", cex.lab = 1.18,
     xlab = "", ylab = "His-codon usage (%)")
axis(1, at = c(1, 2), labels = c("H-x-H His", "non-H-x-H His"))
points(c(1, 2), cac,  pch = 19, cex = 1.5)   # filled = CAC
points(c(1, 2), cat_, pch = 1,  cex = 1.5)   # open   = CAT

## the ONE comparison: CAC at two different x -> horizontal bracket
hs_bracket(1, 2, 61, star_cac)

legend(x = 1.5, y = 38.5, xjust = 0.5, horiz = TRUE, bty = "n", cex = 1.0,
       x.intersp = 0.85, pch = c(19, 1), legend = c("CAC", "CAT"))
mtext("A   His codon CAC vs CAT", side = 3, line = 1.0, adj = 0, font = 2, cex = 1.15)

## ============================ Panel B ============================
## middle residue of de novo H-x-H; aa ordered by mid-frequency descending.
mid  <- mid[order(mid$mid_freq_pct, decreasing = TRUE), ]
xa   <- seq_len(nrow(mid))            # 1..20, Thr at x = 1
thrx <- which(mid$aa == "T")

par(mar = c(4.6, 5.2, 3.0, 1.0))
plot(NA, xlim = c(0.3, 20.7), ylim = c(0, 26), xaxt = "n", cex.lab = 1.18,
     xlab = "amino acid (single-letter)",
     ylab = "frequency as middle residue of H-x-H (%)")
axis(1, at = xa, labels = FALSE)                            # mtext draws every
mtext(mid$aa, side = 1, at = xa, line = 0.7, cex = 1.0)     # letter (no thinning)
segments(xa, mid$bg_freq_pct, xa, mid$mid_freq_pct, lwd = 1.0)
points(xa, mid$bg_freq_pct,  pch = 1,  cex = 1.1)   # open   = genomic background
points(xa, mid$mid_freq_pct, pch = 19, cex = 1.1)   # filled = middle residue

## the ONE comparison: Thr middle-vs-background at a shared x -> topstar
star_y <- mid$mid_freq_pct[thrx] + 1.6
hs_topstar(thrx, star_y, star_thr)
text(thrx + 0.5, star_y, "4.2x background", pos = 4, cex = 0.9)   # thr_fold = 4.18

legend("topright", bty = "n", cex = 1.0, pch = c(19, 1), y.intersp = 1.15,
       legend = c("middle residue", "background"))
mtext("B   middle residue of de novo H-x-H", side = 3, line = 1.0,
      adj = 0, font = 2, cex = 1.15)

mtext("*** p<0.001   ** p<0.01   * p<0.05   ns not significant",
      side = 1, outer = TRUE, line = 0.5, cex = 0.82)
invisible(dev.off())
cat("wrote fig2.pdf\n")
