# Code & data release — de novo neORF histidine/cysteine spacing & the ZMEG zinc story

Reproduces the manuscript figures for the *Drosophila melanogaster* de novo ORF
(neORF) metal-site analysis: His/Cys spacing enrichment, the (CA)n codon-usage
mechanism, de novo ORF transcription, and the ZMEG candidate metal site.

All paths in the code are **relative**. Run everything through `run_all.sh`
(it `cd`s into this `code/` directory first), or run the individual scripts
from within `code/`.

## Manuscript figure numbering

| Fig | content | script | style | reproducible from |
|---|---|---|---|---|
| **1** | H-x-H / His-spacing enrichment (A profile; B H-x-H, C C-x-x-C fold; D incidence; E density) | `fig1_hisspacing.R` | `hs_style.R` (mono) | shipped CSVs ✅ |
| **2** | (CA)n codon usage (A His codon CAC bias; B H-x-H middle residue) | `fig2_camechanism.R` | `hs_style.R` (mono) | shipped CSVs ✅ |
| **3** | (CA)n → H-x-H **mechanism schematic** + poly-HxH multiplicity | `fig3_mechanism.py` | `hs_mpl.py` (mono) | bundled ✓ |
| **4** | de novo ORF **transcription** (STAR) + ZMEG cross-species **trajectory** | `fig4_transcription.py` | inline (color)† | bundled ✓ |
| **5** | ZMEG candidate **bis-His metal site** (fold, Metal3D, His-cage) | `fig5_metalsite.py` (+ `fig5_render_pymol.py`) | inline (color)† | PyMOL + structures ⚠️ |
| **S1** | top-15 enriched ELM/SLiM classes | `figS1_slim.R` | `hs_style.R` (mono) | shipped CSVs ✅ |
| **S2** | amino-acid composition by class | `figS2_aacomp.R` | `hs_style.R` (mono) | shipped CSVs ✅ |
| **S3** | apo-model calibration control (calprotectin vs ZMEG) | `figS_calpro.py` | `hs_mpl.py` (mono) | structures ⚠️ |
| **Dataset S1** | per-replicate expression of the clean-intergenic mappable neORFs (472 identified; 466 quantified — 6 X-linked near-zero loci outside the featureCounts run, see `data/analysis/clean_set_unquantified.tsv`; ZMEG labelled) | `tableS1_expression.py` | — (CSV → `../data/`) | bundled ✓ |
| **Table S1** | ZMEG (CA)n tract length within/across species | `tableS2.py` | `hs_mpl.py` (mono) | bundled ✓ |
| **Table 2** | perfect dinucleotide microsatellite tracts by class (2L/2R/3L/3R, ≥5 units) | `dinuc_census.py` | — (TSV → `../data/hisspacing/table2_dinuc_census.tsv`) | dm6 arms via `fetch_data.sh` ✓ |

`hs_style.R` (R base graphics) and `hs_mpl.py` (matplotlib) are the two shared
**austere/monochrome** style modules; every figure that can be monochrome sources one of them.

† **Style deviation (intentional).** Figs 4 and 5 carry *one* meaningful colour channel
each and do **not** source `hs_mpl.py`: Fig 4 (C/D) colour-codes the ZMEG character
state ((CA)n / ORF / His status) across the tree; Fig 5 uses the pLDDT confidence bands
and a single His-red accent. These are semantic encodings, kept deliberately; the rest of
each figure is greyscale. All other figures are strictly monochrome (shape/line/greyscale only).

## Two reproducibility tracks

### Track 1 — R figures from shipped CSVs (default, tested)

Figs 1, 2, S1, S2 render with **no raw data, no Python, no network**. The R scripts
read `../data/hisspacing/*.csv`:

```
./run_all.sh
```

Full-from-raw variant (regenerates those CSVs with the Python drivers first, then
draws the R figures — path-wired, not tested here; needs network + numpy + scipy):

```
./run_all.sh --from-raw
```

### Track 2 — zinc/ZMEG figures (Figs 3, 4, 5, S3 + Tables S1/S2)

```
./run_all.sh --zinc-figs
```

These are Python and read the **analysis data tree** at `../data/analysis/`
(sibling of `submission/`), which is now **bundled** under `data/analysis/`. Override
the root with `ZROOT=/path/to/analysis`. The mode is best-effort: a missing
input or dependency warns and continues. Python deps: `matplotlib` (+ `numpy`,
`scipy`, `pandas` for Fig 4; `PIL` for Fig 5). Fig 5 additionally needs **PyMOL** and
structure files — render its panels first, then compose:

```
pymol -cq fig5_render_pymol.py     # writes ray-traced panels to $FIG5_SCR (default ./fig5_render)
python3 fig5_metalsite.py          # composes figures/fig5.pdf
```

Inputs each zinc figure expects under `$ZROOT` (`../data/analysis/`):

| script | reads |
|---|---|
| `fig3_mechanism.py` | `neORF.fa` |
| `fig4_transcription.py` | `star_counts/counts.txt`, `neORF_clean_mappable.tsv`, `neORF_clean_dm6.tsv`, `zmeg_history/tables/zmeg_orf_raw.tsv` |
| `tableS1_expression.py` | `star_counts/counts.txt`, `sample_sheet.tsv`, `neORF_clean_mappable.tsv` (→ CSV in `../data/`) |
| `tableS2.py` | `avenueA/tables/task3_tractlength.tsv`, `zmeg_history/tables/zmeg_synteny_summary.tsv` |
| `figS_calpro.py` | `structures/cand3.pdb`, `structures/calprotectin/calprotectin_heterodimer.pdb` |
| `fig5_render_pymol.py` | `$FIG5_SCR/cand3.pdb`, `$FIG5_SCR/calprotectin_heterodimer.pdb` (place these, then render) |

Output PDFs/PNGs land in `figures/` (Fig 5 also honours `FIG5_OUTDIR`).

## Layout

```
submission/
  code/                     <- all scripts (this directory)
    run_all.sh              <- entry point (Track 1 default; --from-raw; --zinc-figs)
    fetch_data.sh           <- downloads raw inputs for the R --from-raw path into ../data/raw/
    fig{1,2}*.R  figS{1,2}*.R          <- R figures (mono; hs_style.R)
    fig{3,4,5}*.py tableS2.py figS_calpro.py    <- zinc/ZMEG figures (hs_mpl.py / inline)
    tableS1_expression.py   <- Dataset S1 generator (machine-readable CSV -> ../data/)
    hs_style.R  hs_mpl.py   <- shared austere style modules (R / Python)
    figures/                <- output PDFs/PNGs (created on run)
  data/
    hisspacing/             <- shipped intermediate CSVs (Figs 1/2/S1/S2 inputs)
    intermediate/           <- precomputed chelation tables the Python drivers read
    raw/                    <- raw inputs, populated by fetch_data.sh (empty here)
    Dataset_S1_neORF_expression.csv  <- Dataset S1: per-replicate STAR TPM (466 of the 472-locus clean set quantified; 6 X-linked near-zero loci unquantified)
```

## R figure inputs (Track 1, `../data/hisspacing/`)

| Figure | R script | reads |
|---|---|---|
| **Figure 1** | `fig1_hisspacing.R` | `fig1_profile.csv`, `fig1_hxh5.csv`, `fig1_cxxc5.csv`, `fig1_stats.csv`, `fig1_hxh5_stats.csv`, `fig1_incidence.csv`, `fig1_density.csv` |
| **Figure 2** | `fig2_camechanism.R` | `fig2_codon.csv`, `fig2_middle.csv`, `fig2_stats.csv` |
| **Figure S1** | `figS1_slim.R` | `figS1_slim.csv`, `figS1_hxh.csv` |
| **Figure S2** | `figS2_aacomp.R` | `figS2_comp.csv`, `figS2_sig.csv` |

## Python drivers (Track 1 `--from-raw` only)

Regenerate the `../data/hisspacing/*.csv` from raw inputs. Run in this order
(encoded in `run_all.sh --from-raw`):

| Script | writes | notes |
|---|---|---|
| `denovo_hisspacing_analysis.py` | `fig1_profile`, `fig2_codon`, `fig2_middle`, `fig2_stats`, `figS2_comp`, `figS2_sig` | main driver; genome/GTF/rmsk intergenic scan + neORF/proteome/metal loading |
| `denovo_fig1_5series.py` | `fig1_hxh5`, `fig1_cxxc5` | 5-series H-x-H / C-x-x-C order-fold; reads back `fig1_profile.csv` |
| `denovo_fig1_bstats.py` | `fig1_hxh5_stats` | bootstrap brackets for Fig 1B/C |
| `denovo_lenmatch_regen.py` | `fig1_stats` | length-matched proteome comparison |
| `denovo_fig1_wholeprotein.py` | `fig1_incidence`, `fig1_density` | reads `../data/intermediate/chelation_derived.json` |
| `denovo_figS1_enriched.py` | `figS1_slim`, `figS1_hxh` | see note below |
| `dinuc_census.py` | `table2_dinuc_census.tsv` | Table 2; reads the dm6 arms `dna_{2L,2R,3L,3R}.fa.gz` from `fetch_data.sh` (or any genome FASTA given as arguments) |

Shared helpers: `hs_common.py` (FASTA/translation/spacing), `dinuc_shuffle.py`
(Altschul–Erikson dinucleotide shuffle null).

**Note on `denovo_figS1_enriched.py`:** it also needs
`../data/intermediate/slim_chelation_results.json`, the output of a separate
pre-registered ELM/SLiM regex scan (`slim_chelation_scan.py`) that is *not* part
of this figure release and is *not* a raw download. It is therefore not fetched
by `fetch_data.sh`. The selected classes are already captured in the shipped
`figS1_slim.csv`; to re-derive them, place that scan's JSON at that path.

## RNA-seq expression (Figure 4, STAR)

Figure 4 panels A/B quantify de novo ORF transcription in *D. melanogaster* S2R+
cells with **STAR + featureCounts** (unique reads, length-normalized to TPM), on a
**mappability-honest, non-overlapping** clean set (`neORF_clean_mappable.tsv`:
472 loci). That set is reproducible from primary inputs via
`rebuild_clean_intergenic.py` (BLAT neORFs → dm6 unique best hit; drop any FlyBase
gene overlap; genmap K=101/E=2 unique-mapping screen). Dataset S1
(`tableS1_expression.py` → `../data/Dataset_S1_neORF_expression.csv`) carries the
per-replicate STAR TPMs (WT / ZnT63-KO / IA2-KO × control/Zn/Mn) for the 466 of 472 quantified in this run (6 X-linked near-zero loci unquantified; see `data/analysis/clean_set_unquantified.tsv`),
with per-condition means and H-x-H content. The featureCounts matrix
(`data/analysis/star_counts/counts.txt`) is **bundled**; regenerate it from raw FASTQ
(SRR5608193–SRR5608210) with `run_star_featurecounts.sh` (+ `build_augmented_gtf.sh`).

## Raw inputs (downloaded by `fetch_data.sh` into `../data/raw/`)

`fetch_data.sh` covers only the **R (`--from-raw`) inputs** — the neORF/proteome/
genome/metal data behind Figs 1/2/S1/S2. The zinc-figure inputs (Track 2) live in
the analysis tree and are not fetched here.

- **Grandchamp neORF nucleotide FASTAs** — Zenodo record 7322757
  (`gc2023/*_DeNovoGenesNuc.fa`).
- **Proteome + GTF** — FlyBase current *D. melanogaster* release
  (`dmel_pep.fa.gz`, `dmel.gtf.gz`).
- **Metal-binders** — UniProt REST, reviewed, organism 7227, keyword KW-0479
  (`dmel_metalbinding.fasta`).
- **dm6 genome arms + RepeatMasker** — UCSC goldenPath dm6
  (`dna_{2L,2R,3L,3R}.fa.gz`, `dm6_rmsk.txt.gz`).

`fetch_data.sh` uses stable, documented URLs; where a filename/version cannot be
pinned without querying the server (the FlyBase release directory and the Zenodo
archive contents) it auto-detects at run time and prints a clearly-marked
`TODO:` fallback rather than guessing.

## Dependencies / tested versions

- **R 4.2** with base graphics + `gplots` (for `plotCI`). No `plotrix` needed.
- **Python 3** (tested 3.9) with `numpy` (and `scipy` for
  `denovo_fig1_wholeprotein.py`), in a conda env.
- **Zinc figures (Track 2):** `matplotlib` (+ `numpy`/`scipy`/`pandas`/`PIL`);
  Fig 5 also needs **PyMOL**.

The R figures-from-CSVs path (Track 1 default) was tested end to end with
`Rscript` (R 4.2) and `gplots`; Figs 1, 2, S1, S2 render. The full-from-raw Python
path is path-wired and imports cleanly but was not executed end to end here (raw
inputs not bundled). The Track-2 zinc figures now read the bundled `data/analysis/` tree
and are reproducible standalone; regenerated in this
release.
