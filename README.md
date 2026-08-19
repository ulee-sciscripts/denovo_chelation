# denovo_chelation

Code and data release for the manuscript on *de novo* ORFs as candidate metal (zinc/manganese) chelators in *Drosophila*.

- `code/` — analysis and figure-generation scripts. `code/run_all.sh` reproduces the figures; see `code/README.md` for details and dependencies.
- `data/` — machine-readable supplementary data:
  - `TableS1_neORF_expression_472.csv` — per-replicate S2R+ expression of the clean-intergenic, uniquely-mappable neORF set (**472 identified**; **466 quantified** in this featureCounts run — 6 X-linked, near-zero loci fell outside the run's annotation, listed in `data/analysis/clean_set_unquantified.tsv`) (H-x-H content, dm6 coordinates, ZMEG labelled).
  - `analysis/` — the bundled analysis inputs the zinc/ZMEG figures read (neORF sequences, featureCounts counts, the clean-set table, cross-species tables, control structures). This makes **all** figures reproducible without any external tree.
- `manuscript_denovo_chelation.docx` — manuscript (Science-style numbered references).

## Reproducibility

Every figure and table regenerates from committed inputs — there is **no dependency on an external `zinc_analysis/` tree** (a previous release read one; it is now bundled under `data/analysis/`).

- Figs 1, 2, S1, S2 render from `data/hisspacing/*.csv` (`run_all.sh`, R only).
- Figs 3, 4, S3 + Tables S1, S2 render from `data/analysis/` (`run_all.sh --zinc-figs`, Python). Fig 5 additionally requires PyMOL-rendered panels (`fig5_render_pymol.py`; input structures are bundled).

### The neORF expression set (472) is now reproducible from primary inputs

The clean-intergenic set consumed by Fig 4 / Table S1 is the output of a documented,
auditable pipeline (previously it shipped as an opaque list):

- `code/rebuild_clean_intergenic.py` — the **"mappability-clean, non-overlapping"** selector: BLAT neORFs → dm6 (unique best hit), require zero FlyBase gene overlap, then genmap K=101/E=2 unique-mapping screen. Emits `neORF_clean_mappable.tsv` (472 loci).
- `code/build_genmap.sh` — regenerates the dm6 mappability track.
- `code/build_augmented_gtf.sh` — builds the featureCounts annotation (FlyBase genes + the clean neORFs).
- `code/run_star_featurecounts.sh` — aligns the 18 S2R+ RNA-seq libraries (**SRR5608193–SRR5608210**, single-end) with STAR and runs `featureCounts -s 0 -t exon -g gene_id` → `star_counts/counts.txt`.

### Coordinates

dm6 coordinates in the shipped tables are the **BLAT-verified** genomic coordinates
(`data/analysis/ref/neorf_dm6.bed`), i.e. neORF sequences re-mapped to dm6. An earlier
metadata layer carried stale pre-liftover coordinates (and gene-overlap labels computed
on them); those are not used here — the clean-set producer recomputes coordinates and
gene overlap from the genome directly.
