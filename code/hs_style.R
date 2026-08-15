## hs_style.R -- shared austere style for the de novo His-spacing / chelation figures.
##
## SPEC (docs/results/FIGURE_WORK.md §4), enforced here so every figure matches:
##   * base plot()/plotCI() only -- NO barplot(), image(), heatmap, etc.
##   * the ONLY channels distinguishing elements are pch, lwd, lty
##   * NO colour, NO fills, NO shading, NO background grids
##   * legends never overlap data (place in a verified-empty region / outer margin)
##   * every categorical comparison gets a horizontal bracket + significance marker
##   * left-to-right order (unless a figure overrides): background -> de novo -> metal
##
## Everything is black. Series are told apart by shape (pch) and line (lwd/lty).
suppressMessages(library(gplots))   # plotCI

## --- stable per-series visual keys, identical across ALL nine figures ---------
## de novo is always the solid dot + solid heavy line (the hero series).
HS_KEY <- list(
  denovo     = list(pch = 19, lwd = 2.6, lty = 1, lab = "neORF"),
  intergenic = list(pch = 1,  lwd = 1.8, lty = 2, lab = "intergenic (repeat-masked)"),
  metal      = list(pch = 17, lwd = 1.8, lty = 3, lab = "metal-binders"),
  randorf    = list(pch = 4,  lwd = 1.5, lty = 4, lab = "random ORF"),
  matched    = list(pch = 5,  lwd = 1.8, lty = 5, lab = "proteome (length-matched)"),
  repeatd    = list(pch = 0,  lwd = 1.8, lty = 3, lab = "repeat-derived")
)
hs_pch <- function(k) HS_KEY[[k]]$pch
hs_lwd <- function(k) HS_KEY[[k]]$lwd
hs_lty <- function(k) HS_KEY[[k]]$lty
hs_lab <- function(k) HS_KEY[[k]]$lab

## --- device with consistent, austere par (black only, no grid) ---------------
hs_pdf <- function(path, width = 7.2, height = 4.8) {
  pdf(path, width = width, height = height)
  par(mar = c(4.8, 5.0, 2.8, 1.0), mgp = c(3.0, 0.8, 0), las = 1,
      cex.axis = 1.18, cex.lab = 1.30, tcl = -0.38, xpd = FALSE)
}

## --- horizontal comparison bracket with significance marker above ------------
## Draws  |___|  from x1..x2 at height y (small downward ticks) and the marker
## (ns / * / ** / ***) centred above the line.  Use for categorical (different-x)
## comparisons.  'ns' is rendered lighter/normal-weight; stars are bold.
hs_bracket <- function(x1, x2, y, star, tick = NULL, lwd = 1.2, cex = 1.35) {
  if (is.null(tick)) tick <- (par("usr")[4] - par("usr")[3]) * 0.018
  segments(x1, y, x2, y, lwd = lwd)
  segments(x1, y, x1, y - tick, lwd = lwd)
  segments(x2, y, x2, y - tick, lwd = lwd)
  ns <- identical(star, "ns")
  text((x1 + x2) / 2, y, star, pos = 3, offset = 0.12,
       font = if (ns) 1 else 2, cex = if (ns) cex * 0.8 else cex)
}

## --- same-x significance star (profiles: two curves compared at one x) --------
## A horizontal bracket is geometrically meaningless when the two compared points
## share an x (e.g. de novo vs intergenic at spacing d).  There we place the
## marker at a fixed top row instead; the compared pair is named in the caption.
hs_topstar <- function(x, y, star, cex = 1.28) {
  ns <- identical(star, "ns")
  text(x, y, star, font = if (ns) 1 else 2, cex = if (ns) cex * 0.72 else cex)
}

## --- point + vertical CI (adds to an existing plot canvas) -------------------
## Thin wrapper on gplots::plotCI, forced monochrome; call after plot(NA, ...).
hs_ci <- function(x, y, lo, hi, pch = 19, lwd = 1.6, cex = 1.2) {
  plotCI(x, y, ui = hi, li = lo, add = TRUE, pch = pch, cex = cex,
         lwd = lwd, col = "black", barcol = "black", sfrac = 0.006, gap = 0)
}

## --- reference line at fold = 1 (no order); dashed, thin, grey-free ----------
hs_hline1 <- function() abline(h = 1, lty = 2, lwd = 1.0)

cat("hs_style.R loaded\n")
