# denovo_chelation

Code and data release for the manuscript on *de novo* ORFs as candidate metal (zinc/manganese) chelators in *Drosophila*.

- `code/` — analysis and figure-generation scripts. `code/run_all.sh` reproduces the figures; see `code/README.md` for details and dependencies.
- `data/` — machine-readable supplementary data (`TableS1_neORF_expression_267.csv`: per-replicate expression of the 267 mappability-clean, non-overlapping neORFs, with H-x-H content and dm6 coordinates).
- `manuscript_denovo_chelation.docx` — manuscript (Science-style numbered references).

Notes: figure scripts reproduce the paper's figures (pixel-verified). The zinc/ZMEG figure scripts read an unbundled analysis tree (`../../zinc_analysis/`); Figure 5 additionally requires PyMOL. See `code/README.md`.
