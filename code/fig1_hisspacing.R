#!/usr/bin/env Rscript
## Manuscript Figure 1 -- the de novo His-spacing panels (A profile; B H-x-H,
## C C-x-x-C fold; D incidence; E density).
##   A  His-pair order-fold vs spacing (profile; full width)
##   B  H-x-H order-fold spectrum      C  C-x-x-C order-fold spectrum
##   D  cluster incidence vs matched   E  per-residue density vs matched
## Two separate legends: a LINE legend for the profile (A) directly under it, and
## a POINT legend for the categorical panels (B-E) at the bottom.  Panel letters
## are floating labels to the left of the centred titles.  Austere spec (hs_style.R).
source("hs_style.R")                # run from the code/ directory (see run_all.sh)
DATA <- file.path("..", "data", "hisspacing")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
p1  <- read.csv(file.path(DATA, "fig1_profile.csv"))
hx  <- read.csv(file.path(DATA, "fig1_hxh5.csv"),   stringsAsFactors = FALSE)
cx  <- read.csv(file.path(DATA, "fig1_cxxc5.csv"),  stringsAsFactors = FALSE)
sp  <- read.csv(file.path(DATA, "fig1_stats.csv"), stringsAsFactors = FALSE)
inc <- read.csv(file.path(DATA, "fig1_incidence.csv"), stringsAsFactors = FALSE)
den <- read.csv(file.path(DATA, "fig1_density.csv"),   stringsAsFactors = FALSE)
inc <- inc[match(c("HxH", "CxxC"), inc$tag), ]
den <- den[match(c("HxH", "CxxC"), den$tag), ]
sp_h <- sp$star[sp$panel == "HxH"]; sp_c <- sp$star[sp$panel == "CxxC"]
bs  <- read.csv(file.path(DATA, "fig1_hxh5_stats.csv"), stringsAsFactors = FALSE)
bstar <- function(k) bs$star[bs$comp == k]

CEX_LET <- 2.05      # floating panel letter (A, B, ...)
CEX_TIT <- 1.58     # centred panel title
ptitle <- function(letter, title, cex_tit = CEX_TIT, line = 1.0) {
  mtext(letter, side = 3, line = line, adj = 0,   font = 2, cex = CEX_LET)
  mtext(title,  side = 3, line = line, adj = 0.5, font = 2, cex = cex_tit)
}

## centred custom legend: line (2 segments w/ gap) only where `lines` is TRUE, plus
## the pch marker and label; optional italic prefix. y in the strip's 0..1 coords.
draw_legend <- function(keys, labs, lines, prefix, cex, symw, gapw, sym_cex, y = 0.6) {
  g <- 0.008
  pw <- if (nzchar(prefix)) strwidth(prefix, cex = cex, font = 3) + 0.016 else 0
  lw <- strwidth(labs, cex = cex) + 0.008
  total <- pw + sum(symw + lw) + gapw * (length(labs) - 1)
  cur <- (1 - total) / 2
  if (nzchar(prefix)) { text(cur, y, prefix, adj = c(0, 0.5), cex = cex, font = 3); cur <- cur + pw }
  for (k in seq_along(keys)) {
    xc <- cur + symw / 2
    if (lines[k]) {
      segments(cur, y, xc - g, y, lty = hs_lty(keys[k]), lwd = hs_lwd(keys[k]))
      segments(xc + g, y, cur + symw, y, lty = hs_lty(keys[k]), lwd = hs_lwd(keys[k]))
    }
    points(xc, y, pch = hs_pch(keys[k]), cex = sym_cex)
    text(cur + symw + 0.006, y, labs[k], adj = c(0, 0.5), cex = cex)
    cur <- cur + symw + lw[k] + gapw
  }
}

hs_pdf(file.path(FIGS, "fig1.pdf"), width = 15.8, height = 11.4)
## rows: A | line-legend(6) | B C D E | point-legend + p-key(7)
layout(matrix(c(1, 1, 1, 1,
                6, 6, 6, 6,
                2, 3, 4, 5,
                7, 7, 7, 7), nrow = 4, byrow = TRUE), heights = c(4.3, 0.8, 4.1, 1.5))

## ================= Panel A : His-pair spacing profile =================
par(mar = c(5.6, 6.4, 3.2, 1.2), cex.axis = 1.62, cex.lab = 1.68)
plot(NA, xlim = c(0.7, 6.3), ylim = c(0.72, 2.74), xaxt = "n", xlab = "",
     ylab = "His-pair order-fold  (obs / scramble)")
axis(1, at = 1:6)
mtext(c("H-H", "H-x-H", "H-x-x-H", "H-x-x-x-H", "H-x-x-x-x-H", "H-x-x-x-x-x-H"),
      side = 1, at = 1:6, line = 2.0, cex = 1.12)
mtext("His-pair spacing  d  (residues apart)", side = 1, line = 4.5, cex = 1.52)
hs_hline1()
lines(p1$d, p1$metal, lty = hs_lty("metal"), lwd = hs_lwd("metal"))
points(p1$d, p1$metal, pch = hs_pch("metal"), cex = 1.45)
lines(p1$d, p1$clean, lty = hs_lty("intergenic"), lwd = hs_lwd("intergenic"))
hs_ci(p1$d, p1$clean, p1$cl_lo, p1$cl_hi, pch = hs_pch("intergenic"), cex = 1.4)
lines(p1$d, p1$denovo, lty = hs_lty("denovo"), lwd = hs_lwd("denovo"))
hs_ci(p1$d, p1$denovo, p1$dn_lo, p1$dn_hi, pch = hs_pch("denovo"), cex = 1.7)
for (i in seq_len(nrow(p1))) hs_topstar(p1$d[i], 2.34, p1$star[i], cex = 1.55)
ptitle("A", "His-pair spacing profile")

## ===== Panels B/C : 5 series, order random -> proteome -> intergenic -> de novo -> metal =====
## pch: random x(4), proteome diamond(5), intergenic circle(1), de novo dot(19), metal tri(17)
pch5 <- c(hs_pch("randorf"), hs_pch("matched"), hs_pch("intergenic"), hs_pch("denovo"), hs_pch("metal"))
cex5 <- c(1.95, 1.95, 1.95, 2.4, 2.2)
par(mar = c(2.8, 6.2, 3.6, 0.8), cex.axis = 1.48, cex.lab = 1.52)
plot(NA, xlim = c(0.5, 5.5), ylim = c(0.9, 2.4), xaxt = "n", xlab = "",
     ylab = "order-fold  (obs / scramble)")
axis(1, at = 1:5, labels = FALSE); hs_hline1()
for (i in 1:5) hs_ci(hx$order[i], hx$fold[i], hx$lo[i], hx$hi[i], pch = pch5[i], cex = cex5[i])
hs_bracket(1, 3, max(hx$hi[1], hx$hi[2], hx$hi[3]) + 0.12, bstar("random_vs_intergenic"), cex = 1.55)
hs_bracket(4, 5, max(hx$hi[4], hx$hi[5]) + 0.12, bstar("denovo_vs_metal"), cex = 1.55)
ptitle("B", "H-x-H order-fold", cex_tit = 1.5)

par(mar = c(2.8, 6.2, 3.6, 0.8), cex.axis = 1.48, cex.lab = 1.52)
plot(NA, xlim = c(0.5, 5.5), ylim = c(0.9, 5.9), xaxt = "n", xlab = "",
     ylab = "order-fold  (obs / scramble)")
axis(1, at = 1:5, labels = FALSE); hs_hline1()
for (i in 1:5) hs_ci(cx$order[i], cx$fold[i], cx$lo[i], cx$hi[i], pch = pch5[i], cex = cex5[i])
hs_bracket(4, 5, max(cx$hi[4], cx$hi[5]) + 0.20, sp_c, cex = 1.55)   # de novo vs metal
ptitle("C", "C-x-x-C order-fold", cex_tit = 1.5)

par(mar = c(4.2, 6.2, 3.6, 0.8), cex.axis = 1.48, cex.lab = 1.52)
plot(NA, xlim = c(0.55, 2.45), ylim = c(0, 0.58), xaxt = "n", xlab = "",
     ylab = "incidence (frac. >=1 cluster)")
axis(1, at = 1:2, labels = FALSE)
mtext(c("H-x-H", "C-x-x-C"), side = 1, at = 1:2, line = 1.6, cex = 1.38)
offD <- 0.22
for (i in 1:2) {
  x <- i
  hs_ci(x - offD, inc$matched[i], inc$matched[i] - inc$matched_sd[i],
        inc$matched[i] + inc$matched_sd[i], pch = hs_pch("matched"), cex = 1.95)
  points(x, inc$denovo[i], pch = hs_pch("denovo"), cex = 2.4)
  points(x + offD, inc$metal[i], pch = hs_pch("metal"), cex = 2.2)
  hs_bracket(x - offD, x, max(inc$denovo[i], inc$matched[i] + inc$matched_sd[i]) + 0.055, inc$star[i], cex = 1.55)
}
ptitle("D", "cluster incidence", cex_tit = 1.5)

par(mar = c(4.2, 6.2, 3.6, 0.8), cex.axis = 1.48, cex.lab = 1.52)
plot(NA, xlim = c(0.55, 2.45), ylim = c(0.05, 0.52), xaxt = "n", xlab = "",
     ylab = "density (per 100 aa)")
axis(1, at = 1:2, labels = FALSE)
mtext(c("H-x-H", "C-x-x-C"), side = 1, at = 1:2, line = 1.6, cex = 1.38)
offE <- 0.22
for (i in 1:2) {
  x <- i
  hs_ci(x - offE, den$mt[i], den$mt_lo[i], den$mt_hi[i], pch = hs_pch("matched"), cex = 1.95)
  hs_ci(x,        den$dn[i], den$dn_lo[i], den$dn_hi[i], pch = hs_pch("denovo"),  cex = 2.4)
  hs_ci(x + offE, den$me[i], den$me_lo[i], den$me_hi[i], pch = hs_pch("metal"),   cex = 2.2)
  hs_bracket(x - offE, x, max(den$dn_hi[i], den$mt_hi[i]) + 0.03, den$star[i], cex = 1.55)
}
ptitle("E", "density (per 100 aa)", cex_tit = 1.5)

## ================= line legend for the profile (Panel A) =================
par(mar = c(0.1, 1.0, 0.1, 1.0)); plot.new()
draw_legend(c("intergenic", "denovo", "metal"),
            c("intergenic (repeat-masked)", "neORF", "metal-binders"),
            c(TRUE, TRUE, TRUE), "", cex = 1.62, symw = 0.062,
            gapw = 0.045, sym_cex = 2.0, y = 0.5)

## ================= point legend for panels B-E + p-value key =================
par(mar = c(0.1, 1.0, 0.1, 1.0)); plot.new()
draw_legend(c("randorf", "matched", "intergenic", "denovo", "metal"),
            c("random ORF", "proteome (length-matched)", "intergenic (repeat-masked)",
              "neORF", "metal-binders"),
            rep(FALSE, 5), "", cex = 1.5, symw = 0.05,
            gapw = 0.022, sym_cex = 2.35, y = 0.72)
text(0.5, 0.24, "*** p<0.001    ** p<0.01    * p<0.05    ns not significant",
     cex = 1.5, adj = c(0.5, 0.5))
invisible(dev.off())
cat("wrote fig1.pdf\n")
