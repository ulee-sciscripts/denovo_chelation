#!/usr/bin/env bash
# build_augmented_gtf.sh -- build the featureCounts annotation:
#   FlyBase genes  +  the clean-intergenic neORF loci (as single-exon features).
# Produces $OUT (default /tmp/augmented_v2.gtf; large, not committed).
# Consumed by run_star_featurecounts.sh. neORF gene_ids match the clean set.
#
# Requires: a FlyBase dm6 GTF ($FLY_GTF) and the clean set
# (../data/analysis/neORF_clean_mappable.tsv).
set -euo pipefail
cd "$(dirname "$0")"
BUNDLE="${ZROOT:-../data/analysis}"
FLY_GTF="${FLY_GTF:?set FLY_GTF to a FlyBase dm6 .gtf (see fetch_data.sh)}"
CLEAN="${CLEAN_SET:-$BUNDLE/neORF_clean_mappable.tsv}"
OUT="${OUT:-/tmp/augmented_v2.gtf}"
# FlyBase gene/exon/etc lines
awk -F'\t' '$2=="FlyBase"' "$FLY_GTF" > "$OUT"
# append neORF features (gene+transcript+exon+CDS) from clean-set dm6 coords
python3 - "$CLEAN" >> "$OUT" <<'PY'
import csv,sys
for r in csv.DictReader(open(sys.argv[1]),delimiter="\t"):
    c,s,e,st,g=r["chrom"],r["start_dm6"],r["end_dm6"],r["strand"],r["name"]
    attr=f'gene_id "{g}"; transcript_id "{g}.t"; gene_biotype "denovo_neORF";'
    for feat in ("gene","transcript","exon","CDS"):
        fr="0" if feat=="CDS" else "."
        a=attr+(' exon_number "1";' if feat in("exon","CDS") else "")
        print(f"{c}\tneORF\t{feat}\t{s}\t{e}\t.\t{st}\t{fr}\t{a}")
PY
echo "wrote $OUT ($(grep -c neORF "$OUT") neORF feature lines)"
