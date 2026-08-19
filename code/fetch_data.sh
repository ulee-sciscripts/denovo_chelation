#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# fetch_data.sh -- download the heavy RAW inputs for the full-from-raw pipeline
# into ../data/raw/ .  These are NOT needed to redraw the figures from the
# shipped intermediate CSVs (see run_all.sh, default mode); they are only needed
# to regenerate those CSVs with the Python drivers (run_all.sh --from-raw).
#
# SCOPE: this covers the R figures only (Figs 1, 2, S1, S2).  The zinc/ZMEG
# figures (Figs 3, 4, 5 + Tables) read the bundled analysis data at ../data/analysis/
# (star_counts/, zmeg_history/, structures/, neORF.fa, ...), which is NOT bundled
# and NOT fetched here; see README.md, "Track 2".
#
# Sources (all public, documented, stable):
#   Grandchamp neORF nucleotide FASTAs .. Zenodo record 7322757
#   D. melanogaster proteome + GTF ....... FlyBase (current dmel release)
#   metal-binders (UniProt KW-0479) ...... UniProt REST (reviewed, organism 7227)
#   dm6 genome (2L,2R,3L,3R) + RepeatMasker  UCSC goldenPath dm6
#
# The drivers expect these exact filenames under ../data/raw/ :
#   gc2023/*_DeNovoGenesNuc.fa   dmel_pep.fa.gz   dmel.gtf.gz
#   dmel_metalbinding.fasta      dm6_rmsk.txt.gz  dna_{2L,2R,3L,3R}.fa.gz
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"
RAW="../data/raw"
mkdir -p "$RAW/gc2023"

# small download helper (prefers curl, falls back to wget)
dl () {  # dl <url> <output-path>
  local url="$1" out="$2"
  if [ -s "$out" ]; then echo "  exists, skip: $out"; return 0; fi
  echo "  downloading: $url"
  if command -v curl >/dev/null 2>&1; then curl -fSL "$url" -o "$out"
  else wget -O "$out" "$url"; fi
}

# --- 1. UCSC dm6 genome (per chromosome arm) + RepeatMasker -------------------
# UCSC per-chromosome FASTA headers (>chr2L) are not used by the code; only the
# sequence is read, so chr2L.fa.gz is simply saved as dna_2L.fa.gz.
UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/dm6"
echo "[UCSC] dm6 genome arms + rmsk"
for c in 2L 2R 3L 3R; do
  dl "$UCSC/chromosomes/chr${c}.fa.gz" "$RAW/dna_${c}.fa.gz"
done
# RepeatMasker table; code reads genoName (chr2L) and strips the 'chr' prefix.
dl "$UCSC/database/rmsk.txt.gz" "$RAW/dm6_rmsk.txt.gz"

# --- 2. FlyBase proteome + GTF (current release) -----------------------------
# The release directory is version-stamped (e.g. dmel_r6.61); auto-detect it
# from the 'current' symlink.  dm6 == FlyBase R6, so FlyBase gene coordinates and
# the UCSC dm6 repeat coordinates share one coordinate system.
echo "[FlyBase] proteome + GTF (current release)"
FB_BASE="https://ftp.flybase.net/releases/current"
FB_REL="$(curl -fsSL "$FB_BASE/" 2>/dev/null | grep -oE 'dmel_r[0-9.]+' | sort -u | head -1 || true)"
if [ -z "${FB_REL:-}" ]; then
  echo "  TODO: could not auto-detect the current FlyBase release directory."
  echo "        Browse $FB_BASE/ , pick the dmel_rX.XX dir, and set FB_REL below."
  # FB_REL="dmel_r6.61"   # <-- uncomment & set to a real release if auto-detect fails
else
  echo "  detected FlyBase release: $FB_REL"
  # translated proteome (all isoforms) and the gene-model GTF
  dl "$FB_BASE/$FB_REL/fasta/dmel-all-translation-${FB_REL#dmel_}.fasta.gz" "$RAW/dmel_pep.fa.gz"
  dl "$FB_BASE/$FB_REL/gtf/dmel-all-${FB_REL#dmel_}.gtf.gz"                  "$RAW/dmel.gtf.gz"
fi

# --- 3. UniProt metal-binders (reviewed, D. melanogaster, KW-0479) -----------
echo "[UniProt] metal-binders (KW-0479, reviewed, organism 7227)"
UNIPROT="https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%29+AND+%28organism_id%3A7227%29+AND+%28keyword%3AKW-0479%29"
dl "$UNIPROT" "$RAW/dmel_metalbinding.fasta"

# --- 4. Grandchamp 2023 neORF nucleotide FASTAs (Zenodo 7322757) -------------
# The archive layout is not hard-coded here: query the Zenodo API for the file
# list, download every file, then move any *_DeNovoGenesNuc.fa into gc2023/.
echo "[Zenodo] Grandchamp neORF FASTAs (record 7322757)"
ZEN_API="https://zenodo.org/api/records/7322757"
if command -v curl >/dev/null 2>&1 && command -v python3 >/dev/null 2>&1; then
  urls="$(curl -fsSL "$ZEN_API" | python3 -c 'import sys,json; [print(f["links"]["self"]) for f in json.load(sys.stdin).get("files",[])]' 2>/dev/null || true)"
  if [ -n "$urls" ]; then
    tmp="$RAW/gc2023/_download"; mkdir -p "$tmp"
    while read -r u; do [ -n "$u" ] && dl "$u" "$tmp/$(basename "${u%%\?*}")"; done <<< "$urls"
    # unpack any archives, then collect the neORF nucleotide FASTAs
    for a in "$tmp"/*.zip;    do [ -e "$a" ] && (cd "$tmp" && unzip -o "$a" >/dev/null) || true; done
    for a in "$tmp"/*.tar.gz; do [ -e "$a" ] && tar -xzf "$a" -C "$tmp" || true; done
    find "$tmp" -name '*_DeNovoGenesNuc.fa' -exec mv -t "$RAW/gc2023/" {} + 2>/dev/null || true
  fi
fi
if ! ls "$RAW"/gc2023/*_DeNovoGenesNuc.fa >/dev/null 2>&1; then
  echo "  TODO: could not auto-place the neORF FASTAs."
  echo "        Download record 7322757 from https://zenodo.org/record/7322757 ,"
  echo "        and put the seven *_DeNovoGenesNuc.fa files in $RAW/gc2023/ ."
fi

echo "fetch_data.sh done. Raw inputs are in $RAW/ ."
