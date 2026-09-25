#!/usr/bin/env Rscript
## Combined summary figure in the manuscript's austere style (hs_style.R: base graphics, black only,
## series distinguished by pch, point + 95% bootstrap CI, no grid, no legend -- labels under each tick).
## n per category, significance stars and test descriptions live in the figure caption.
## Input:  results/fig_motif_enrichment.csv, results/fig_prevalence.csv (from fig_stats.py)
## Output: figures/fig_summary.{pdf,png}
##   A-C  order-fold (obs / within-sequence scramble) for H-x-H, C-x-C, C-x-x-C
##        stars: Fisher's exact, 2x2 of position pairs at the motif spacing (motif vs not; observed vs expected)
##   D-E  Cys and His prevalence (% of residues)
##        stars: binomial exact, residue count vs intergenic residue frequency ("ref")
source("../hs_style.R")
D <- "."
me <- read.csv(file.path(D, "results", "fig_motif_enrichment.csv"), stringsAsFactors = FALSE)
pv <- read.csv(file.path(D, "results", "fig_prevalence.csv"), stringsAsFactors = FALSE)

## order: background -> neORF -> young -> older -> metal-binders (n per category lives in the caption)
KEYS <- c("intergenic", "neORF", "young_all", "older", "metal")
X    <- 1:5
PCH  <- c(1, 19, 0, 5, 17)
CEX  <- c(1.8, 2.1, 1.9, 1.9, 2.0)
LABS <- c("intergenic", "neORFs", "young", "older", "metal-binders")

CEX_LET <- 2.0; CEX_TIT <- 1.6
ptitle <- function(letter, title) {
  mtext(letter, side = 3, line = 1.0, adj = 0, font = 2, cex = CEX_LET)
  mtext(title,  side = 3, line = 1.0, adj = 0.5, font = 2, cex = CEX_TIT)
}
## rotated, right-justified multi-line tick labels: R centres the lines of a rotated multi-line string
## on each other, so draw each line separately; line ends share one y and are centred on the tick.
rot_labels <- function(x, y, labs, srt = 45, cex = 1.35, nudge = 0.08) {   # nudge (in): shift multi-line labels right
  th <- srt * pi / 180
  lh <- strheight("M", units = "inches", cex = cex) * 1.55            # line pitch (inches)
  u  <- par("usr"); pin <- par("pin")
  for (i in seq_along(x)) {
    ln <- strsplit(labs[i], "\n")[[1]]
    n <- length(ln)
    for (k in seq_along(ln)) {
      ## all line ends share one y; successive lines shift sideways by lh / sin(th) (keeps the
      ## perpendicular line spacing = lh), and the set of line ends is centred on the tick
      dx <- (((k - 1) - (n - 1) / 2) * lh / sin(th) + if (n > 1) nudge else 0) * diff(u[1:2]) / pin[1]
      dy <- 0
      text(x[i] + dx, y + dy, ln[k], srt = srt, adj = c(1, 0.5), xpd = NA, cex = cex)   # anchor = end of text midline
    }
  }
}
panel <- function(vals, lo, hi, stars, ylim, ylab, letter, title, hline = TRUE, no_ci = integer(0), labs = LABS) {
  par(cex = 1, mar = c(7.2, 6.6, 3.4, 0.8), mgp = c(4.0, 0.8, 0), cex.axis = 1.3, cex.lab = 1.45)   # cex=1: undo layout() auto-shrink
  plot(NA, xlim = c(0.5, 5.5), ylim = ylim, xaxt = "n", xlab = "", ylab = ylab)
  axis(1, at = X, labels = FALSE)
  ## labels end just below each tick, rotated
  yl <- par("usr")[3] - (par("usr")[4] - par("usr")[3]) * 0.045
  rot_labels(X, yl, labs, srt = 45, cex = 1.35)
  if (hline) hs_hline1()
  for (i in seq_along(X)) {
    if (i %in% no_ci) points(X[i], vals[i], pch = PCH[i], cex = CEX[i])      # reference: no error bar
    else hs_ci(X[i], vals[i], lo[i], hi[i], pch = PCH[i], cex = CEX[i])
  }
  pad <- (ylim[2] - ylim[1]) * 0.045
  for (i in seq_along(X)) if (stars[i] != "ref") hs_topstar(X[i], hi[i] + pad, stars[i], cex = 1.45)
  ptitle(letter, title)
}
get_m <- function(m) { d <- me[me$motif == m, ]; d[match(KEYS, d$key), ] }
get_p <- function(r) { d <- pv[pv$residue == r, ]; d[match(KEYS, d$key), ] }

draw <- function() {
  par(oma = c(1.0, 1.5, 0, 0))
  ## top: A B C (2 cols each) ; bottom: D E (3 cols each) ; star key
  layout(matrix(c(1, 1, 2, 2, 3, 3,
                  4, 4, 4, 5, 5, 5,
                  6, 6, 6, 6, 6, 6), nrow = 3, byrow = TRUE), heights = c(5.8, 5.8, 0.45))
  mm <- list(c("HxH", "H-x-H"), c("CxC", "C-x-C"), c("CxxC", "C-x-x-C"))
  for (j in 1:3) {
    d <- get_m(mm[[j]][1])
    panel(d$fold, d$lo, d$hi, d$fisher_star,
          if (mm[[j]][1] == "CxxC") c(0.5, 5.8) else c(0.5, 2.3),
          "order-fold  (obs / scramble)", LETTERS[j], paste(mm[[j]][2], "order-fold"))
  }
  rr <- list(c("Cys", "Cysteine prevalence"), c("His", "Histidine prevalence"))
  for (j in 1:2) {
    d <- get_p(rr[[j]][1])
    panel(d$pct, d$lo, d$hi, d$binom_star, c(1.5, 4.4), "% of residues", LETTERS[3 + j], rr[[j]][2],
          hline = FALSE, no_ci = 1,   # intergenic = binomial reference
          labs = replace(LABS, 1, "intergenic\n(reference)"))
  }
  par(mar = c(0.1, 1.0, 0.1, 1.0)); plot.new()
  text(0.5, 0.5, "*** p<0.001    ** p<0.01    * p<0.05    ns not significant", cex = 1.5)
}

dir.create(file.path("..", "figures"), showWarnings = FALSE)
f <- file.path("..", "figures", "fig_summary")
hs_pdf(paste0(f, ".pdf"), width = 17.0, height = 13.6); draw(); invisible(dev.off())
png(paste0(f, ".png"), width = 17.0, height = 13.6, units = "in", res = 300)
par(las = 1, mgp = c(3.0, 0.8, 0), tcl = -0.38); draw(); invisible(dev.off())
cat("wrote figures/fig_summary.{pdf,png}\n")
