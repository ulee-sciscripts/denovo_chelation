#!/usr/bin/env bash
# Rasterize the vector figure PDFs to 300-dpi PNG + LZW-compressed TIFF (submission-grade).
# Uses ghostscript only (correct 300-dpi resolution metadata). Output: figures/raster_300dpi/.
set -euo pipefail
cd "$(dirname "$0")"
OUT=figures/raster_300dpi; mkdir -p "$OUT"
# R vector figures always; the Python figs (fig3-5, tables, calpro) already emit
# 300-dpi PNGs but are re-rastered here too when present (all guarded by -f).
for f in figures/fig1.pdf figures/fig2.pdf figures/figS1.pdf figures/figS2.pdf \
         figures/fig3.pdf figures/fig4.pdf figures/fig5.pdf \
         figures/tableS2.pdf figures/figS_calpro_control.pdf; do
  [ -f "$f" ] || continue
  b=$(basename "$f" .pdf); o="$OUT/${b}_300dpi"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -o "${o}.png" "$f"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=tifflzw -r300 -o "${o}.tif" "$f"
done
echo "300-dpi rasters (PNG + LZW TIFF) written to $OUT/"
