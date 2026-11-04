#!/usr/bin/env python3
"""
Genome-wide dinucleotide microsatellite census for D. melanogaster (dm6).

This is the reproducible, in-house basis for the manuscript's Table 2 (abundance of
perfect dinucleotide microsatellite tracts by class) and the claim that (CA)n / AC is
the most abundant dinucleotide repeat. The literature (e.g. Ross et al. 2003) reports
the same qualitative dominance; here we compute our OWN number directly from dm6.

What it does
------------
Scans each chromosome's (+)-strand sequence for PERFECT dinucleotide microsatellite
tracts and tabulates them into the 4 canonical strand-and-rotation-collapsed classes:

    AC class = {AC, CA, GT, TG}   (all the same repeat up to rotation + revcomp)
    AG class = {AG, GA, CT, TC}
    AT class = {AT, TA}
    CG class = {CG, GC}

For each class and each minimum-unit threshold n_min we report:
    (1) tract count      -- number of maximal perfect tracts
    (2) total repeat bp  -- summed tract lengths
and the AC-class share of all dinucleotide microsatellite content, by BOTH metrics.

Counting logic (the double-counting trap)
-----------------------------------------
We reduce the problem to an unambiguous, overlap-free primitive: a "dinucleotide
repeat tract" is a MAXIMAL alternating run of two DISTINCT bases (period-2 periodicity),
e.g. ...ACACACAC... or ...GTGTGTG...  We detect these on the (+) strand only.

We work in "bonds" (adjacent base pairs). Bond i is the unordered pair {s[i], s[i+1]}.
A maximal run of consecutive bonds that all carry the SAME unordered heterodinucleotide
pair is, provably, a single period-2 alternating tract: two consecutive same-pair bonds
share their middle base, which forces strict A/C/A/C alternation (a same-pair bond cannot
be a homopolymer, so no A/A ambiguity). Thus one maximal same-pair bond-run == one tract,
with no phase ambiguity: (AC)n and the embedded (CA)(n-1) are the SAME run, counted once.

Because we scan only the (+) strand, the AC class already captures both a (CA)n tract and,
at a different locus, a (GT)n tract (GT is the reverse complement of AC). We do NOT also
reverse-complement / separately scan GT, so there is no strand double-counting.

The 6 unordered heterodinucleotide pairs collapse onto the 4 classes:
    {A,C} -> AC     {G,T} -> AC   (revcomp of AC)
    {A,G} -> AG     {C,T} -> AG   (revcomp of AG)
    {A,T} -> AT
    {C,G} -> CG
Homopolymer bonds ({X,X}) and any bond touching an N are treated as tract breaks.

Tract geometry: a bond-run of k bonds spans k+1 bases -> bp = k+1, units = (k+1)//2.
A tract qualifies at threshold n_min iff units >= n_min. We additionally clip the rare
case of two qualifying tracts that abut and share a single boundary base (assigning the
shared base to the earlier tract) so that tracts are strictly disjoint; this is enforced
by asserting that summed tract bp per class never exceeds the chromosome length.

Scope: the paper (Table 2) uses the major euchromatic arms 2L/2R/3L/3R, n_min = 5 units
(>= 10 bp). If the input also contains other sequences (X/4/Y/scaffolds, e.g. a full
FlyBase dm6 FASTA), the whole assembly is reported as a secondary scope so the choice is
explicit. The census is run on the RAW (unmasked) genome.

Inputs
------
    python3 dinuc_census.py                    # default: the per-arm UCSC files that
                                               # fetch_data.sh writes to ../data/raw/
                                               # (dna_2L.fa.gz ... dna_3R.fa.gz)
    python3 dinuc_census.py genome.fa [...]    # any FASTA(s), plain or .gz
UCSC-style headers (>chr2L) and FlyBase-style headers (>2L) are both accepted.

Writes: ../data/hisspacing/table2_dinuc_census.tsv  and prints a summary.

Std-lib + numpy only, deterministic.
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hs_common import read_fasta  # yields (header, seq) from plain or .gz FASTA

# ---- config ----------------------------------------------------------------
ARMS = ["2L", "2R", "3L", "3R"]          # major euchromatic arms (primary)
NMINS = [5, 7]                            # 5 units (>=10 bp) primary; 7 matches text
CLASSES = ["AC", "AG", "AT", "CG"]

RAW = os.path.join(HERE, os.pardir, "data", "raw")
DEFAULT_FASTAS = [os.path.join(RAW, f"dna_{a}.fa.gz") for a in ARMS]   # from fetch_data.sh
OUT = os.path.join(HERE, os.pardir, "data", "hisspacing")
os.makedirs(OUT, exist_ok=True)
OUT_TSV = os.path.join(OUT, "table2_dinuc_census.tsv")

t0 = time.time()
def log(*a):
    print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)

# ---- base coding -----------------------------------------------------------
BASE2INT = np.full(256, 4, np.uint8)      # A/C/G/T -> 0..3 ; everything else -> 4 (N)
for i, b in enumerate("ACGT"):
    BASE2INT[ord(b)] = i

# unordered heterodinucleotide pair id (sorted lo*4+hi) -> class name
PAIR2CLASS = {
    0 * 4 + 1: "AC",   # {A,C}
    2 * 4 + 3: "AC",   # {G,T} revcomp of AC
    0 * 4 + 2: "AG",   # {A,G}
    1 * 4 + 3: "AG",   # {C,T} revcomp of AG
    0 * 4 + 3: "AT",   # {A,T}
    1 * 4 + 2: "CG",   # {C,G}
}


def find_tracts(seq):
    """Scan a (+)-strand sequence; return {class: list of unit-counts}.

    Detects maximal period-2 alternating dinucleotide runs (see module docstring).
    Every base belongs to at most one tract (shared boundary base -> earlier tract).
    """
    n = len(seq)
    per_class_units = {c: [] for c in CLASSES}
    if n < 3:
        return per_class_units
    codes = BASE2INT[np.frombuffer(seq.encode("ascii", "replace"), np.uint8)]

    a = codes[:-1].astype(np.int16)       # bond first base
    b = codes[1:].astype(np.int16)        # bond second base
    valid = (a < 4) & (b < 4) & (a != b)  # heterodinuc, no N, no homopolymer
    lo = np.minimum(a, b)
    hi = np.maximum(a, b)
    pcode = np.where(valid, lo * 4 + hi, -1).astype(np.int16)  # per-bond pair id

    # maximal runs of identical pcode
    m = pcode.size
    change = np.empty(m, bool)
    change[0] = True
    np.not_equal(pcode[1:], pcode[:-1], out=change[1:])
    run_starts = np.flatnonzero(change)                 # bond indices where a run starts
    run_vals = pcode[run_starts]
    run_ends = np.empty_like(run_starts)                # exclusive bond index
    run_ends[:-1] = run_starts[1:]
    run_ends[-1] = m

    # NOTE on boundaries: consecutive bond-runs always share exactly one boundary base
    # (the last bond of one run and the first bond of the next share a base). For a tract
    # flanked by ordinary unique/non-repeat sequence the flanking bonds are sub-threshold
    # singletons and are ignored, so an isolated tract keeps its full k+1 bp. Only when two
    # RECORDED (qualifying) tracts directly abut do we clip the shared base, assigning it to
    # the earlier tract, so recorded tracts are strictly disjoint (asserted downstream).
    prev_base_end = 0                                   # base coord; only recorded tracts update it
    for rs, re, pc in zip(run_starts.tolist(), run_ends.tolist(), run_vals.tolist()):
        if pc < 0:
            continue
        cls = PAIR2CLASS[pc]
        base_start = rs
        base_end = re + 1                               # exclusive; bonds [rs,re) span bases [rs,re]
        if base_start < prev_base_end:                  # abuts a previously RECORDED tract
            base_start = prev_base_end
        bp = base_end - base_start
        units = bp // 2
        if units >= NMINS[0]:                           # keep anything >= smallest threshold
            per_class_units[cls].append((units, bp))
            prev_base_end = base_end                    # only recorded tracts reserve bases
    return per_class_units


def census_over(seqs):
    """seqs: list of (name, seq). Return dict keyed by n_min -> {class: (count, bp)}
    and the raw AC unit list (units) for the length distribution, plus total scanned bp."""
    # accumulate all tracts as (class, units, bp)
    all_units = {c: [] for c in CLASSES}
    total_len = 0
    for name, seq in seqs:
        seq = seq.upper()
        total_len += len(seq)
        pcu = find_tracts(seq)
        # disjointness sanity: summed tract bp per class <= chromosome length
        for c in CLASSES:
            s = sum(bp for (_u, bp) in pcu[c])
            assert s <= len(seq), f"{name}/{c}: tract bp {s} > chrom len {len(seq)}"
            all_units[c].extend(pcu[c])

    res = {}
    for nmin in NMINS:
        tab = {}
        for c in CLASSES:
            kept = [(u, bp) for (u, bp) in all_units[c] if u >= nmin]
            cnt = len(kept)
            bp = sum(x[1] for x in kept)
            tab[c] = (cnt, bp)
        res[nmin] = tab
    # AC unit-length distribution (at the primary threshold) is returned separately
    ac_units_primary = sorted(u for (u, _bp) in all_units["AC"] if u >= NMINS[0])
    return res, ac_units_primary, total_len


def pct(x, tot):
    return (100.0 * x / tot) if tot else float("nan")


def summarize(label, res, ac_units, total_len, rows_out):
    print()
    print("=" * 72)
    print(f"SCOPE: {label}   (scanned {total_len:,} bp)")
    print("=" * 72)
    for nmin in NMINS:
        tab = res[nmin]
        tot_cnt = sum(tab[c][0] for c in CLASSES)
        tot_bp = sum(tab[c][1] for c in CLASSES)
        print(f"\n  n_min = {nmin} units (>= {2*nmin} bp)")
        print(f"  {'class':<6} {'tracts':>10} {'%tracts':>9} {'repeat_bp':>12} {'%bp':>8}")
        for c in CLASSES:
            cnt, bp = tab[c]
            print(f"  {c:<6} {cnt:>10,} {pct(cnt,tot_cnt):>8.2f}% {bp:>12,} {pct(bp,tot_bp):>7.2f}%")
            rows_out.append([label, nmin, c, cnt, f"{pct(cnt,tot_cnt):.4f}",
                             bp, f"{pct(bp,tot_bp):.4f}"])
        print(f"  {'ALL':<6} {tot_cnt:>10,} {'100.00':>8}% {tot_bp:>12,} {'100.00':>7}%")
        rows_out.append([label, nmin, "ALL", tot_cnt, "100.0000", tot_bp, "100.0000"])
        ac_cnt, ac_bp = tab["AC"]
        print(f"    -> AC class share: {pct(ac_cnt,tot_cnt):.2f}% by tract count, "
              f"{pct(ac_bp,tot_bp):.2f}% by repeat bp   "
              f"(top-by-count={max(CLASSES,key=lambda c:tab[c][0])}, "
              f"top-by-bp={max(CLASSES,key=lambda c:tab[c][1])})")
    # AC length distribution (primary n_min)
    if ac_units:
        arr = np.array(ac_units)
        in79 = int(((arr >= 7) & (arr <= 9)).sum())
        print(f"\n  AC tract-length distribution (n_min={NMINS[0]} units): "
              f"n={arr.size:,}  min={arr.min()}  median={int(np.median(arr))}  "
              f"max={arr.max()} units")
        print(f"    fraction of AC tracts in 7-9 unit range: "
              f"{in79:,}/{arr.size:,} = {pct(in79, arr.size):.1f}%")


def seq_name(header):
    """'chr2L ...' (UCSC) or '2L ...' (FlyBase) -> '2L'."""
    name = header.split()[0]
    return name[3:] if name.startswith("chr") else name


def main():
    fastas = sys.argv[1:] or DEFAULT_FASTAS
    missing = [f for f in fastas if not os.path.exists(f)]
    if missing:
        sys.exit(f"missing input FASTA(s): {missing}\n"
                 f"run ./fetch_data.sh first, or pass genome FASTA path(s) as arguments.")
    seqs = []
    for fa in fastas:
        log(f"reading genome: {fa}")
        seqs.extend((seq_name(h), s) for h, s in read_fasta(fa))
    by_name = dict(seqs)
    arms_present = [a for a in ARMS if a in by_name]
    log(f"loaded {len(seqs)} sequences; arms present: {arms_present}")
    if arms_present != ARMS:
        sys.exit(f"expected all of {ARMS} in the input; found {arms_present}")

    rows_out = []

    # ---- primary scope (Table 2): major euchromatic arms 2L/2R/3L/3R -------
    arm_seqs = [(a, by_name[a]) for a in ARMS]
    log("scanning euchromatic arms 2L/2R/3L/3R ...")
    res_a, ac_a, len_a = census_over(arm_seqs)
    summarize("euchromatic arms (2L,2R,3L,3R)", res_a, ac_a, len_a, rows_out)

    # ---- secondary scope: whole assembly, only if the input has more than the arms
    if len(seqs) > len(ARMS):
        log("scanning whole assembly (all sequences) ...")
        res_g, ac_g, len_g = census_over(seqs)
        summarize("whole genome (all sequences, incl. X/4)", res_g, ac_g, len_g, rows_out)

    # ---- write TSV ---------------------------------------------------------
    with open(OUT_TSV, "w") as fh:
        fh.write("scope\tn_min_units\tclass\ttract_count\tpct_tracts\trepeat_bp\tpct_bp\n")
        for r in rows_out:
            fh.write("\t".join(str(x) for x in r) + "\n")
    log(f"wrote {OUT_TSV}")


if __name__ == "__main__":
    main()
