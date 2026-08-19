#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_all.sh -- reproduce the manuscript figures.
#
# Manuscript figure numbering (see README.md):
#   Fig 1  H-x-H / His-spacing enrichment ............ fig1_hisspacing.R   (R)
#   Fig 2  (CA)n codon-usage mechanism ............... fig2_camechanism.R  (R)
#   Fig 3  (CA)n -> H-x-H mechanism schematic ........ fig3_mechanism.py   (py)
#   Fig 4  de novo ORF transcription + ZMEG trajectory fig4_transcription.py (py)
#   Fig 5  ZMEG candidate metal site ................. fig5_metalsite.py   (py, +PyMOL)
#   Fig S1 enriched ELM/SLiM classes ................. figS1_slim.R        (R)
#   Fig S2 amino-acid composition .................... figS2_aacomp.R      (R)
#   Fig S3 calprotectin apo-model control ............ figS_calpro.py      (py)
#   Table S1 neORF expression (472 neORFs, CSV) ...... tableS1_expression.py (-> ../data/)
#   Table S2 ZMEG (CA)n tract across strains/species . tableS2.py          (py)
#
# Execution modes:
#   (default)      figures-from-CSVs: run ONLY the R scripts for Figs 1, 2, S1, S2,
#                  which read the shipped intermediate CSVs in ../data/hisspacing/ .
#                  No raw data, no Python, no network.  This is the tested path.
#
#   --from-raw     full pipeline for the R figures: download raw inputs
#                  (fetch_data.sh), run the Python drivers to regenerate the
#                  intermediate CSVs, then run the R scripts.  Requires network +
#                  numpy + scipy.  Path-wired but NOT tested in this release.
#
#   --zinc-figs    generate the zinc/ZMEG figures (Figs 3, 4, 5 + Table S1/S2 +
#                  the calprotectin control).  These are Python and read the
#                  analysis data tree at ../data/analysis/ (neORF.fa, star_counts/,
#                  zmeg_history/, structures/, ...), which is now BUNDLED in this
#                  figure release.  Requires matplotlib (+ numpy/scipy/pandas/PIL;
#                  Fig 5 also needs PyMOL-rendered panels, see fig5_render_pymol.py).
#                  Override the data root with ZROOT=/path ; best-effort (a missing
#                  input/dep warns and continues rather than aborting the batch).
#
# Interpreters can be overridden:  PYTHON=... RSCRIPT=... ./run_all.sh
# Figures are written to ./figures/ .
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
RSCRIPT="${RSCRIPT:-Rscript}"
MODE="${1:-figures}"

run_r_figures () {
  echo "== drawing R figures (Figs 1, 2, S1, S2) =="
  "$RSCRIPT" fig1_hisspacing.R           # Figure 1  (panels A-E)
  "$RSCRIPT" fig2_camechanism.R          # Figure 2  (A codon usage, B middle residue)
  "$RSCRIPT" figS1_slim.R                # Figure S1 (top-15 enriched ELM SLiM classes)
  "$RSCRIPT" figS2_aacomp.R              # Figure S2 (amino-acid composition by class)
}

if [ "$MODE" = "--zinc-figs" ]; then
  echo "== zinc/ZMEG figures (Figs 3, 4, 5, S3) + Tables S1/S2 =="
  echo "   data root ZROOT=${ZROOT:-../data/analysis} (bundled; standalone)"
  # best-effort: a missing input/dependency warns but does not abort the batch
  "$PYTHON" fig3_mechanism.py       || echo "WARN fig3_mechanism.py failed (missing data/deps?)"
  "$PYTHON" fig4_transcription.py   || echo "WARN fig4_transcription.py failed (missing data/deps?)"
  "$PYTHON" figS_calpro.py          || echo "WARN figS_calpro.py failed (missing data/deps?)"
  "$PYTHON" tableS1_expression.py   || echo "WARN tableS1_expression.py failed (missing data/deps?)"   # Table S1 -> ../data/TableS1_neORF_expression_472.csv
  "$PYTHON" tableS2.py              || echo "WARN tableS2.py failed (missing data/deps?)"                # Table S2 -> figures/tableS2.{tsv,png,pdf}
  echo "   Fig 5 (fig5_metalsite.py) needs PyMOL-rendered panels first:"
  echo "     pymol -cq fig5_render_pymol.py   # writes panels to \$FIG5_SCR (default ./fig5_render)"
  echo "     python3 fig5_metalsite.py        # composes figures/fig5.pdf"
  "$PYTHON" fig5_metalsite.py       || echo "WARN fig5_metalsite.py skipped (render panels first)"
  echo "Done. Zinc figures are in ./figures/ ; the Table S1 CSV is in ../data/ ."
  exit 0
fi

if [ "$MODE" = "--from-raw" ]; then
  echo "== [1/3] fetching raw inputs =="
  bash fetch_data.sh
  echo "== [2/3] running Python drivers (regenerate ../data/hisspacing/*.csv) =="
  # Dependency order: the main driver writes fig1_profile.csv, which fig1_5series
  # and lenmatch_regen read back; the wholeprotein/enriched drivers read shipped intermediates.
  "$PYTHON" denovo_hisspacing_analysis.py   # fig1_profile, fig2_codon/middle/stats, figS2_comp/sig
  "$PYTHON" denovo_fig1_5series.py          # fig1_hxh5, fig1_cxxc5
  "$PYTHON" denovo_fig1_bstats.py           # fig1_hxh5_stats
  "$PYTHON" denovo_lenmatch_regen.py        # fig1_stats
  "$PYTHON" denovo_fig1_wholeprotein.py     # fig1_incidence, fig1_density
  "$PYTHON" denovo_figS1_enriched.py        # figS1_slim, figS1_hxh (needs slim scan JSON, see script)
  echo "== [3/3] =="
else
  echo "== figures-from-CSVs mode (default); pass --from-raw for the full R pipeline, --zinc-figs for Figs 3-5 =="
fi

run_r_figures

echo "Done. Figures are in ./figures/ .  (Figs 3, 4, 5 + Tables: run './run_all.sh --zinc-figs')"
