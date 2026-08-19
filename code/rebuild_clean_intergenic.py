#!/usr/bin/env python3
"""
rebuild_clean_intergenic.py
==============================================================================
AUDITABLE producer of the "clean intergenic, uniquely-mappable" neORF set that
Fig 4 and Table S1 consume (data/analysis/neORF_clean_mappable.tsv, 472 loci).

This script exists so the neORF expression set is REPRODUCIBLE from primary
inputs rather than shipped as an opaque list. Given the neORF nucleotide FASTA
(Grandchamp 2023, Zenodo 7322757), the dm6 genome, a FlyBase gene annotation,
and a genmap mappability track, it regenerates the set end-to-end with every
threshold explicit:

  1. MAP  : BLAT each neORF -> dm6; keep the single best hit by score, requiring
            identity >= MIN_IDENT, query-coverage >= MIN_QCOV, and a UNIQUE
            placement (2nd-best score < UNIQ_FRAC * best) -> uniquely dm6-mapped.
  2. GENE : bedtools-intersect the best-hit spans against FlyBase genes;
            keep genehits == 0  -> clean intergenic (non-overlapping).
  3. MAP2 : genmap K=101 E=2 frac_unique >= MAP_THRESH over the locus span
            -> mappability-clean.

This is the "mappability-clean, non-overlapping" definition the manuscript
states; the shipped set is its output. (An earlier release shipped a curated
list without this producer; this closes that reproducibility gap.)

Paths resolve via env vars, else bundle-relative defaults; tools via $BLAT /
$BEDTOOLS or PATH (shutil.which). Build the genmap track with build_genmap.sh.

Usage:  python3 rebuild_clean_intergenic.py
Outputs (next to the bundle): neORF_clean_mappable.tsv  (drop-in for Fig4/TableS1)
"""
import csv, os, shutil, subprocess, collections

HERE   = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.environ.get("ZROOT", os.path.join(HERE, "..", "data", "analysis"))

GENOME    = os.environ.get("DM6_GENOME",  os.path.join(BUNDLE, "ref", "dm6_genome.fa"))
NEORF_FA  = os.environ.get("NEORF_FA",    os.path.join(BUNDLE, "neORF.fa"))
FLY_GENES = os.environ.get("FLY_GENES_BED", os.path.join(BUNDLE, "ref", "flybase_genes.bed"))
GENMAP_BG = os.environ.get("GENMAP_BG",   os.path.join(BUNDLE, "ref", "dm6_k101_e2.bedgraph"))
META      = os.environ.get("NEORF_META",  os.path.join(BUNDLE, "neORF_metadata.tsv"))
OUT       = os.environ.get("CLEAN_OUT",   os.path.join(BUNDLE, "neORF_clean_mappable.tsv"))
SCR       = os.environ.get("SCRATCH",     os.path.join(HERE, "_rebuild_scratch"))

BLAT     = os.environ.get("BLAT")     or shutil.which("blat")
BEDTOOLS = os.environ.get("BEDTOOLS") or shutil.which("bedtools")

# ---- explicit, auditable thresholds ----------------------------------------
MIN_IDENT, MIN_QCOV, UNIQ_FRAC = 0.90, 0.80, 0.90
MAP_THRESH, UNIQ_MAP           = 0.50, 0.9999

def need(tool, name):
    if not tool:
        raise SystemExit(f"ERROR: {name} not found (set ${name.upper()} or add to PATH)")
    return tool

def main():
    need(BLAT, "blat"); need(BEDTOOLS, "bedtools")
    for p in (GENOME, NEORF_FA, FLY_GENES, GENMAP_BG):
        if not os.path.exists(p):
            raise SystemExit(f"ERROR: missing input {p} (see build_genmap.sh / build_augmented_gtf.sh / fetch_data.sh)")
    os.makedirs(SCR, exist_ok=True)

    # 1. BLAT all neORFs -> dm6
    psl = os.path.join(SCR, "neorf_vs_dm6.psl")
    if not (os.path.exists(psl) and os.path.getsize(psl) > 1000):
        print(f"[1] BLAT {NEORF_FA} vs dm6 (minIdentity={int(MIN_IDENT*100)}) ...", flush=True)
        subprocess.run([BLAT, "-t=dna", "-q=dna", f"-minIdentity={int(MIN_IDENT*100)}",
                        "-noHead", GENOME, NEORF_FA, psl], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hits = collections.defaultdict(list)
    with open(psl) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 17: continue
            match, mism = int(f[0]), int(f[1]); q, qsize = f[9], int(f[10])
            hits[q].append(dict(score=match - mism, ident=match / max(1, match + mism),
                                qcov=(int(f[12]) - int(f[11])) / qsize, strand=f[8],
                                chrom=f[13], tstart=int(f[15]), tend=int(f[16])))
    mapped = {}
    for q, hs in hits.items():
        hs.sort(key=lambda h: h["score"], reverse=True); best = hs[0]
        if best["ident"] < MIN_IDENT or best["qcov"] < MIN_QCOV: continue
        if len(hs) > 1 and hs[1]["score"] >= UNIQ_FRAC * best["score"]: continue
        mapped[q] = best
    print(f"[1] uniquely dm6-mapped: {len(mapped)}", flush=True)

    mapped_bed = os.path.join(SCR, "neorf_dm6.bed")
    with open(mapped_bed, "w") as out:
        for name, h in sorted(mapped.items()):
            out.write(f'{h["chrom"]}\t{h["tstart"]}\t{h["tend"]}\t{name}\t0\t{h["strand"]}\n')

    # 2. gene overlap
    res = subprocess.run([BEDTOOLS, "intersect", "-a", mapped_bed, "-b", FLY_GENES, "-c"],
                         capture_output=True, text=True, check=True)
    genehits = {f.split("\t")[3]: int(f.split("\t")[6]) for f in res.stdout.splitlines()}
    clean = {n: mapped[n] for n in mapped if genehits.get(n, 0) == 0}
    print(f"[2] clean intergenic (genehits==0): {len(clean)}", flush=True)

    # 3. genmap mappability screen
    clean_bed = os.path.join(SCR, "clean.bed")
    with open(clean_bed, "w") as out:
        for name, h in sorted(clean.items()):
            out.write(f'{h["chrom"]}\t{h["tstart"]}\t{h["tend"]}\t{name}\t0\t{h["strand"]}\n')
    res = subprocess.run([BEDTOOLS, "intersect", "-a", clean_bed, "-b", GENMAP_BG, "-wo"],
                         capture_output=True, text=True, check=True)
    agg = collections.defaultdict(lambda: [0, 0])
    for line in res.stdout.splitlines():
        f = line.split("\t"); name = f[3]
        agg[name][0] += int(f[10]); agg[name][1] += int(f[10]) if float(f[9]) >= UNIQ_MAP else 0
    mappable = {}
    for name, h in clean.items():
        span = h["tend"] - h["tstart"]; cov, uniq = agg.get(name, [0, 0])
        if span and uniq / span >= MAP_THRESH:
            mappable[name] = dict(h, frac_unique=round(uniq / span, 4), span_nt=span)
    print(f"[3] clean intergenic AND mappable (frac_unique>={MAP_THRESH}): {len(mappable)}", flush=True)

    # 4. join seq-derived His metadata, write
    meta = {r["name"]: r for r in csv.DictReader(open(META), delimiter="\t")} if os.path.exists(META) else {}
    cols = ["name", "chrom", "start_dm6", "end_dm6", "strand", "len_aa",
            "has_hxh", "poly_hxh", "n_hxh", "frac_unique"]
    with open(OUT, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t"); w.writeheader()
        for name, h in sorted(mappable.items()):
            m = meta.get(name, {})
            w.writerow(dict(name=name, chrom=h["chrom"], start_dm6=h["tstart"] + 1, end_dm6=h["tend"],
                            strand=h["strand"], len_aa=m.get("len_aa", ""), has_hxh=m.get("has_hxh", ""),
                            poly_hxh=m.get("poly_hxh", ""), n_hxh=m.get("n_hxh", ""),
                            frac_unique=h["frac_unique"]))
    print(f"[out] {OUT}  ({len(mappable)} loci)")

if __name__ == "__main__":
    main()
