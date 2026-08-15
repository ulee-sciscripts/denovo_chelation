#!/usr/bin/env Rscript
## Fig S1 -- ELM / short-linear-motif (SLiM) classes most ENRICHED in neORFs.
## Horizontal dot plot; 15 classes ordered so the MOST enriched (largest
## neORF/established fold) sits at the TOP.  Three series per row at the same y
## (no jitter): neORF, established proteome, random ORF.
## The comparison is de novo vs established (all BH-FDR q<0.05).  A bracket per
## row (15 of them) would be unreadable, so the shared result is stated once as a
## caption and each row's `star` sits at the right margin instead.
## Austere spec (hs_style.R): base plot() only; black only; series told apart by
## pch alone; legend in its own bottom strip.  Text sized to match the enlarged
## Fig 1 exemplar (class labels inherit cex.axis 1.18; notes/subtitle ~0.9).
source("hs_style.R")                # run from the code/ directory (see run_all.sh)
DATA <- file.path("..", "data", "hisspacing")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
s   <- read.csv(file.path(DATA, "figS1_slim.csv"), stringsAsFactors = FALSE)
hxh <- read.csv(file.path(DATA, "figS1_hxh.csv"),  stringsAsFactors = FALSE)

## ELM classes ordered largest-fold (most enriched) at the top of the ELM block;
## the H-x-H motif is pinned as a separate reference row above them.
s  <- s[order(s$fold_dn_over_est, decreasing = FALSE), ]
n  <- nrow(s)                          # 15 ELM classes at y = 1..n
yy <- seq_len(n)
yh <- n + 1                            # H-x-H reference row at the top

xmax <- max(s$denovo, s$established, s$randorf,
            hxh$denovo, hxh$established, hxh$randorf) * 1.10
hs_pdf(file.path(FIGS, "figS1.pdf"), width = 9.6, height = 8.2)
par(oma = c(0, 0, 2.2, 0))             # top outer margin holds the title
layout(matrix(c(1, 2), nrow = 2), heights = c(6.8, 0.7))

## ---- main dot-plot panel ----
par(mar = c(4.6, 12.6, 1.0, 1.2))
plot(NA, xlim = c(0, xmax), ylim = c(0.4, yh + 0.6), yaxt = "n",
     xlab = "motif density (regex hits per 100 aa)", ylab = "")
axis(2, at = yy, labels = s$class, las = 1, tick = FALSE)             # ELM classes
axis(2, at = yh, labels = "H-x-H", las = 1, tick = FALSE, font = 2)   # reference motif (bold)
axis(1)
abline(h = n + 0.5, lwd = 1.1)         # divider: H-x-H (above) vs ELM classes (below)

## points at shared y (established -> random -> neORF last, hero on top)
ally <- c(yy, yh)
points(c(s$established, hxh$established), ally, pch = 2,                 cex = 1.2)  # established
points(c(s$randorf,     hxh$randorf),     ally, pch = hs_pch("randorf"), cex = 1.2)  # random ORF
points(c(s$denovo,      hxh$denovo),      ally, pch = hs_pch("denovo"),  cex = 1.4)  # neORF (hero)

mtext("SLiM classes most enriched in neORFs (H-x-H shown for reference)",
      outer = TRUE, side = 3, line = 0.5, font = 2, cex = 1.15)

## ---- bottom legend strip (own panel: cannot overlap data) ----
par(mar = c(0.2, 5.0, 0.2, 1.2)); plot.new()
legend("center", horiz = TRUE, bty = "n", cex = 1.0, pt.cex = 1.35, x.intersp = 0.8,
       legend = c("neORF", "established proteome", "random ORF"),
       pch = c(hs_pch("denovo"), 2, hs_pch("randorf")))
invisible(dev.off())
cat("wrote figS1.pdf\n")
