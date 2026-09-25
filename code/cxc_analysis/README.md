# cxc_analysis — translation of the (CA)ₙ complement, and gene-age stratification

Analysis code added in the GBE revision for the Results section *"Translation of the
complement of (CA)ₙ"* and the new supplementary figure (C-x-C / (TG)ₙ and gene-age
stratification of motif order-fold and residue prevalence), plus the revision's
intergenic-vs-proteome and bis-histidine-spacing exact-test scripts.

## Run convention
Run every script **from inside this directory** (`code/cxc_analysis/`). Paths resolve
relative to the script's own location (`_HERE`): shared modules in `../` (`hs_common.py`,
`hs_style.R`, `denovo_hisspacing_analysis.py`), bundled analysis inputs in
`../../data/analysis/` and `../../data/intermediate/`, and external inputs (not committed)
under `./data/` (see below).

## Reproducing the supplementary figure (self-contained)
The figure is reproducible from committed derived tables — no external data needed:

```
cd code/cxc_analysis
Rscript fig_summary.R      # reads results/fig_motif_enrichment.csv + results/fig_prevalence.csv
                           # writes ../figures/fig_summary.{pdf,png}
```

`fig_stats.py` regenerates those two CSVs from the raw inputs below; the CSVs are committed
so the figure builds without them.

## External inputs (NOT committed — populate `code/cxc_analysis/data/`)
These large/third-party inputs are git-ignored. Fetch them into `code/cxc_analysis/data/`
(mirroring the sub-paths the scripts expect) to re-run the upstream analyses:

| input | expected path | source |
|---|---|---|
| FlyBase proteome | `data/md_hp6_scratch/data/dmel_pep.fa.gz` | FlyBase (dmel r6.61) |
| FlyBase metal-binders | `data/md_hp6_scratch/data/dmel_metalbinding.fasta` | UniProt KW-0479 |
| Grandchamp neORF nucleotide FASTAs | `data/md_hp6_scratch/gc2023/*_DeNovoGenesNuc.fa` | Grandchamp et al. 2023 (Zenodo 7322757) |
| NCBI orthologs + taxonomy | `data/` (build) → `data/gene_age_ncbi.json` | `gene_age_ncbi.py` from NCBI `gene_orthologs.gz` + Taxonomy |
| Heames 2020 de novo set | `data/heames2020/heames2020_dmel_all_orphans.csv` | Heames et al. 2020 |
| intergenic ORF intervals | `data/interg_intervals/` (gene + RepeatMasker interval dir) | built from dm6 GTF + UCSC rmsk |
| dm6 genome / CDS | `../../data/analysis/ref/` (`dm6_genome.fa`, `dmel_cds.fa.gz`) | UCSC/FlyBase; see `../fetch_data.sh`, `../build_*.sh` |

The upstream scripts (`fig_stats.py`, `cxc_analysis.py`, `cxc_registers.py`,
`ca_tg_tracts_4sets.py`, `cxxc_*.py`, `age_motif_foldenrich.py`, `gene_age_ncbi.py`,
`heames_denovo_motifs.py`, `young_motif_scan.py`, `composition_table.py`,
`val_middle_cxc.py`, `origin_composition.py`) were validated in the original analysis
environment; only `fig_summary.R` is verified to run from the committed data here.

## Files
- `fig_summary.R` — the supplementary figure (5 panels: H-x-H / C-x-C / C-x-x-C order-fold; Cys/His prevalence).
- `fig_stats.py` — computes `results/fig_motif_enrichment.csv` and `results/fig_prevalence.csv` (Fisher's exact on spacing position-pairs; binomial exact on residue prevalence vs intergenic; bootstrap CIs).
- `length_quartile_hxh.py` — H-x-H incidence by length quartile for intergenic ORFs and neORFs, each with an amino-acid-shuffle expected incidence (P(≥1 motif) under composition- and length-preserving residue permutation), obs/exp, a Fisher's-exact p (observed vs shuffle-expected with/without ≥1 motif, per quartile), and His prevalence. Shows the incidence gradient is length-driven (flat His prevalence) with a modest ordering excess in neORFs/intergenic (obs/exp>1) but not the proteome (obs/exp≤1). Writes `results/length_quartile_hxh.csv`; caches the intergenic scan to `data/interg_LNH.npy`.
- `cxc_analysis.py`, `cxc_registers.py` — C-x-C / (CA)ₙ-(TG)ₙ reading-frame register order-fold.
- `ca_tg_tracts_4sets.py` — (CA)ₙ vs (TG)ₙ coding-strand tract counts (≥4 units; binomial vs 50:50).
- `cxxc_prevalence.py`, `cxxc_shuffle_incidence.py`, `age_motif_foldenrich.py` — C-x-x-C and gene-age-stratified motif enrichment.
- `gene_age_ncbi.py` — builds `data/gene_age_ncbi.json` from NCBI orthologs + taxonomy.
- `heames_denovo_motifs.py`, `young_motif_scan.py`, `composition_table.py`, `val_middle_cxc.py`, `origin_composition.py` — supporting analyses.
- `spacing_fisher_tests.py` — bis-histidine spacing incidence Fisher's exact against two composition-aware nulls (random-ORF + dinucleotide Altschul–Erikson shuffle).
- `intergenic_vs_proteome.py`, `intergenic_stratified.py`, `interg_canonical.py`, `interg_neORFmatch.py`, `can_fingerprint.py` — intergenic-vs-proteome H-x-H enrichment and the (CA)ₙ codon-usage fingerprint (revision R2 lines 180–182).
- `results/` — committed derived tables for the figure.
