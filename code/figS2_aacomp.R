#!/usr/bin/env Rscript
## Fig S2 -- amino-acid composition across FIVE sets, austere (no color/fill/grid).
## TWO panels in one figure + shared legend strip; points told apart by pch only.
##   A  grouped points for the H,C,D,E class, sets ordered LEFT->RIGHT
##      intergenic, random ORF, de novo, proteome, metal-binders (rule 6); each
##      residue carries a de-novo-vs-proteome bracket (two-proportion z-test, ***).
##   B  full-composition overview: all 20 aa as a horizontal dot plot (H,C,D,E bold).
source("hs_style.R")                # run from the code/ directory (see run_all.sh)
DATA <- file.path("..", "data", "hisspacing")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
d   <- read.csv(file.path(DATA, "figS2_comp.csv"), check.names = FALSE, stringsAsFactors = FALSE)
sig <- read.csv(file.path(DATA, "figS2_sig.csv"),  stringsAsFactors = FALSE)

## per-set pch (task override: proteome = 2, not the matched-key 5); de novo largest cex
P_IN <- 1; P_RO <- 4; P_DN <- 19; P_PR <- 2; P_MB <- 17
SETCOL <- c("intergenic", "random ORF", "de novo", "proteome", "metal-binders")
SETOFF <- c(-0.28, -0.14, 0.00, 0.14, 0.28)     # rule-6 left->right ordering
SETPCH <- c(P_IN, P_RO, P_DN, P_PR, P_MB)
SETCEX <- c(1.00, 1.00, 1.30, 1.00, 1.00)

hs_pdf(file.path(FIGS, "figS2.pdf"), width = 10.5, height = 10.8)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(3.0, 6.2, 1.0))

## ================= PANEL A -- H,C,D,E grouped (rule-6 order) =================
resid  <- c("H", "C", "D", "E")
rlab   <- c("His (H)", "Cys (C)", "Asp (D)", "Glu (E)")
ridx   <- match(resid, d$aa)
gx     <- seq_along(resid)                       # group centres 1..4

par(mar = c(3.4, 4.8, 4.8, 1.2))
plot(NA, xlim = c(0.6, 4.4), ylim = c(0, 8.0), xaxt = "n", xlab = "",
     ylab = "% of residues")
axis(1, at = gx, labels = rlab)

for (s in seq_along(SETCOL))
  points(gx + SETOFF[s], d[[SETCOL[s]]][ridx], pch = SETPCH[s], cex = SETCEX[s])

## de novo (offset 0) vs proteome (offset +0.14) bracket, above the group's tallest point
for (j in seq_along(resid)) {
  top <- max(as.numeric(d[ridx[j], SETCOL]))
  st  <- sig$star[match(resid[j], sig$aa)]
  hs_bracket(gx[j], gx[j] + 0.14, top + 0.45, st, tick = 0.15)
}

mtext("A   His / Cys / Asp / Glu residues", side = 3, line = 0.5, cex = 1.15, adj = 0, font = 2)
mtext("*** p<0.001", side = 3, line = 0.5, cex = 0.82, adj = 1)

## ================= PANEL B -- full composition, all 20 aa ====================
n    <- nrow(d)                       # 20 amino acids
ypos <- n:1                           # first row (A) at top, last (Y) at bottom
aa   <- d$aa
bold <- aa %in% resid                 # H,C,D,E bold

par(mar = c(4.6, 3.2, 2.4, 1.2))
plot(NA, xlim = c(0, 12), ylim = c(0.5, n + 0.5), xaxt = "n", yaxt = "n",
     xlab = "residue frequency (% of residues)", ylab = "")
axis(1, at = seq(0, 12, 2))
axis(2, at = ypos, labels = FALSE)                                              # tick marks
axis(2, at = ypos[bold],  labels = aa[bold],  font = 2, tick = FALSE, las = 1)  # H,C,D,E bold
axis(2, at = ypos[!bold], labels = aa[!bold], font = 1, tick = FALSE, las = 1)

## points, one set at a time; de novo (hero) drawn last so it sits on top
points(d[["intergenic"]],    ypos, pch = P_IN, cex = 1.00)
points(d[["random ORF"]],    ypos, pch = P_RO, cex = 1.00)
points(d[["proteome"]],      ypos, pch = P_PR, cex = 1.00)
points(d[["metal-binders"]], ypos, pch = P_MB, cex = 1.00)
points(d[["de novo"]],       ypos, pch = P_DN, cex = 1.30)

## (the de novo vs proteome comparison is bracketed in Panel A; Panel B is a clean
## full-composition overview -- H, C, D, E are highlighted by bold row labels.)
mtext("B   full composition (all 20 aa)", side = 3, line = 0.6, cex = 1.15, adj = 0, font = 2)

## ================= shared legend strip (own panel, full width) ===============
par(mar = c(0.1, 1.2, 0.1, 1.2)); plot.new()
legend("center", horiz = TRUE, bty = "n", cex = 1.0, x.intersp = 0.6,
       legend = c("intergenic", "random ORF", "neORF", "proteome", "metal-binders"),
       pch    = c(P_IN, P_RO, P_DN, P_PR, P_MB),
       pt.cex = c(1.00, 1.00, 1.30, 1.00, 1.00))
invisible(dev.off())
cat("wrote figS2.pdf\n")
