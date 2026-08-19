#!/usr/bin/env bash
# run_star_featurecounts.sh -- regenerate the neORF expression counts that Fig 4
# and Table S1 read (../data/analysis/star_counts/counts.txt).
#
# Aligns the 18 S2R+ RNA-seq libraries (Wang et al. GEO; SRR5608193..SRR5608210,
# single-end) to dm6 with STAR, then featureCounts against the augmented GTF
# (FlyBase genes + clean-intergenic neORFs) with unique reads only.
# This reproduces the exact featureCounts invocation used for the release:
#   featureCounts -s 0 -t exon -g gene_id -a augmented.gtf <18 sorted BAMs>
#
# Requires: STAR, subread (featureCounts), samtools, aria2c (or wget); a dm6
# genome FASTA; and augmented.gtf from build_augmented_gtf.sh. FASTQs are pulled
# from ENA (faster/steadier than NCBI SRA). ulimit -n is raised (STAR BAM-sort
# opens many temp files; the default 1024 causes a fatal sort error).
set -uo pipefail
cd "$(dirname "$0")"
BUNDLE="${ZROOT:-../data/analysis}"
WORK="${WORK:-/tmp/neorf_requant}"
GENOME="${DM6_GENOME:-$BUNDLE/ref/dm6_genome.fa}"
GTF="${AUG_GTF:?set AUG_GTF to augmented.gtf from build_augmented_gtf.sh}"
IDX="${STAR_INDEX:-$WORK/star_index}"
THREADS="${THREADS:-$(nproc)}"
OUT="${OUT:-$BUNDLE/star_counts/counts.txt}"
SRRS="SRR5608193 SRR5608194 SRR5608195 SRR5608196 SRR5608197 SRR5608198 SRR5608199 SRR5608200 SRR5608201 SRR5608202 SRR5608203 SRR5608204 SRR5608205 SRR5608206 SRR5608207 SRR5608208 SRR5608209 SRR5608210"
mkdir -p "$WORK/fastq" "$WORK/bam" "$(dirname "$OUT")"
ulimit -n 65535 || true

# STAR index (small dm6 genome -> --genomeSAindexNbases 12)
if [ ! -f "$IDX/SAindex" ]; then
  mkdir -p "$IDX"
  STAR --runMode genomeGenerate --genomeDir "$IDX" --genomeFastaFiles "$GENOME" \
       --sjdbGTFfile "$GTF" --sjdbOverhang 100 --genomeSAindexNbases 12 --runThreadN "$THREADS"
fi

# download (ENA direct FASTQ) + align
BAMS=""
for SRR in $SRRS; do
  fq="$WORK/fastq/$SRR.fastq.gz"
  if [ ! -s "$fq" ]; then
    url=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR&result=read_run&fields=fastq_ftp&format=tsv" | tail -1 | cut -f2 | tr ';' '\n' | head -1)
    case "$url" in http://*|https://*|ftp://*) : ;; *) url="https://$url" ;; esac
    if command -v aria2c >/dev/null; then aria2c -x8 -s8 -c -d "$WORK/fastq" -o "$SRR.fastq.gz" "$url"; else wget -c -O "$fq" "$url"; fi
  fi
  bam="$WORK/bam/$SRR.Aligned.sortedByCoord.out.bam"
  if [ ! -s "$bam" ]; then
    STAR --runThreadN "$THREADS" --genomeDir "$IDX" --readFilesIn "$fq" --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$WORK/bam/$SRR." \
         --limitBAMsortRAM 24000000000 --outBAMsortingBinsN 50
  fi
  BAMS="$BAMS $bam"
done

# featureCounts (unique reads, unstranded, exon/gene_id) -- matches the release
featureCounts -a "$GTF" -o "$OUT" -T "$THREADS" -s 0 -t exon -g gene_id $BAMS
echo "wrote $OUT"
