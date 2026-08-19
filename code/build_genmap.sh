#!/usr/bin/env bash
# build_genmap.sh -- regenerate the dm6 mappability track used by the
# clean-intergenic mappability screen (rebuild_clean_intergenic.py step 3).
# Produces ../data/analysis/ref/dm6_k101_e2.bedgraph (genmap, K=101, E=2).
#
# Requires: genmap (conda install -c bioconda genmap), and the dm6 genome FASTA
# (../data/analysis/ref/dm6_genome.fa; fetch via fetch_data.sh / UCSC or FlyBase).
set -euo pipefail
cd "$(dirname "$0")"
BUNDLE="${ZROOT:-../data/analysis}"
GENOME="${DM6_GENOME:-$BUNDLE/ref/dm6_genome.fa}"
IDX="${GENMAP_IDX:-$BUNDLE/ref/genmap_index}"
OUTPFX="${GENMAP_OUT:-$BUNDLE/ref/dm6_k101_e2}"
GENMAP="${GENMAP:-genmap}"
[ -f "$GENOME" ] || { echo "ERROR: genome $GENOME missing (see fetch_data.sh)"; exit 1; }
# genmap requires the index dir to NOT pre-exist
rm -rf "$IDX"
"$GENMAP" index -F "$GENOME" -I "$IDX"
"$GENMAP" map -K 101 -E 2 -I "$IDX" -O "$OUTPFX" -bg
echo "wrote ${OUTPFX}.bedgraph"
