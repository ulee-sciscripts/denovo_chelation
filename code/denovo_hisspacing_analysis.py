#!/usr/bin/env python3
"""
De novo His-spacing / chelation analysis (drives the base-R figures).

Outputs JSON + CSVs under ../data/hisspacing/ for the R plotting layer.
All p-values are bootstrap or permutation (stated per figure); nothing fabricated.

Sets:
  - de novo neORFs (Grandchamp 2023, translated by us; 7,849 unique)
  - clean intergenic ORFs   = 6-frame stop-to-stop ORFs (>=10 aa) of
                              genome \ (gene UNION repeat)
  - repeat-derived intergenic ORFs = ORFs of (repeat \ gene)
  - established proteome (FlyBase)
  - known metal-binders (UniProt KW-0479, 856)
  - random ORFs at intergenic GC=0.413 (genetic-code null; for composition)
"""
import os, sys, glob, json, gzip, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hs_common import (SCRATCH, CODON2AA, AAS, read_fasta, translate, translate_full,
                       load_denovo, load_proteins, count_spacings,
                       observed_spacing_totals, analytic_expected_totals)
from dinuc_shuffle import altschul_erikson_shuffle

OUT = os.path.join(HERE, os.pardir, "data", "hisspacing")
os.makedirs(OUT, exist_ok=True)
SEED = 20260812
DMAX = 6
CHRS = ["2L", "2R", "3L", "3R"]
BIN_EDGES = [10, 40, 80, 160, 400]           # 4 bins: [10,40),[40,80),[80,160),[160,400)
BIN_LABELS = ["10-40", "40-80", "80-160", "160-400"]

t0 = time.time()
def log(*a):
    print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)

# ---- vectorised codon LUT (bases A=0 C=1 G=2 T=3) --------------------------
BASE2INT = np.full(256, 4, np.uint8)
for i, b in enumerate("ACGT"):
    BASE2INT[ord(b)] = i
AA_LUT = np.zeros(64, np.uint8)
for i0, a in enumerate("ACGT"):
    for i1, b in enumerate("ACGT"):
        for i2, c in enumerate("ACGT"):
            AA_LUT[i0 * 16 + i1 * 4 + i2] = ord(CODON2AA[a + b + c])
ORD_H, ORD_STAR = ord("H"), ord("*")

# =====================================================================
# Part 1: per-sequence His-spacing tables (store L,h,o1..o6 for h>=2 ORFs)
# =====================================================================
def spacing_row(prot):
    """Return (L, h, [o1..o6]) for a peptide (numpy)."""
    aa = np.frombuffer(prot.encode(), np.uint8)
    L = aa.size
    hmask = aa == ORD_H
    h = int(hmask.sum())
    o = np.zeros(DMAX, np.int64)
    if h >= 2:
        for d in range(1, DMAX + 1):
            if L > d:
                o[d - 1] = int(np.count_nonzero(hmask[:-d] & hmask[d:]))
    return L, h, o

def rows_from_prots(prots):
    """Build a table (n,8): L,h,o1..o6 for peptides with h>=2 (others add 0/0)."""
    rows = []
    for p in prots:
        L, h, o = spacing_row(p)
        if h >= 2 and L >= 10:
            rows.append((L, h, *o.tolist()))
    return np.array(rows, dtype=np.int64) if rows else np.zeros((0, 8), np.int64)

def analytic_exp_rows(tab, d):
    """Analytic scramble expectation per row for spacing d: (L-d) h(h-1)/(L(L-1))."""
    L = tab[:, 0].astype(float); h = tab[:, 1].astype(float)
    return np.where(L > d, (L - d) * h * (h - 1) / (L * (L - 1)), 0.0)

def pooled_fold(tab, d):
    o = tab[:, 1 + d].astype(float).sum()
    e = analytic_exp_rows(tab, d).sum()
    return (o / e) if e > 0 else np.nan, o, e

# =====================================================================
# Part 2: intergenic ORF scanning (6-frame stop-to-stop)
# =====================================================================
class Accum:
    def __init__(self):
        self.rows = []                     # L,h,o1..o6 for h>=2, L>=10
        self.comp = np.zeros(256, np.int64)  # aa char-code histogram (non-stop)
        self.nbp = 0
        self.norf = 0                      # number of ORFs >=10 aa (any h)

    def table(self):
        return np.array(self.rows, dtype=np.int64) if self.rows else np.zeros((0, 8), np.int64)

def _scan_frame(aa, acc):
    """aa: uint8 array of amino-acid codes ('*' = stop). Update accumulator."""
    is_stop = aa == ORD_STAR
    # composition over non-stop residues
    nonstop = aa[~is_stop]
    if nonstop.size:
        acc.comp += np.bincount(nonstop, minlength=256)
    # stop-to-stop runs
    idx = np.flatnonzero(is_stop)
    bounds = np.empty(idx.size + 2, np.int64)
    bounds[0] = -1; bounds[-1] = aa.size; bounds[1:-1] = idx
    for k in range(bounds.size - 1):
        s = bounds[k] + 1; e = bounds[k + 1]
        L = e - s
        if L < 10:
            continue
        acc.norf += 1
        seg = aa[s:e]
        hmask = seg == ORD_H
        h = int(hmask.sum())
        if h < 2:
            continue
        o = np.zeros(DMAX, np.int64)
        for d in range(1, DMAX + 1):
            if L > d:
                o[d - 1] = int(np.count_nonzero(hmask[:-d] & hmask[d:]))
        acc.rows.append((L, h, *o.tolist()))

def scan_region_seq(dna_bytes, acc):
    """6-frame translate one contiguous DNA region and scan for ORFs."""
    x = BASE2INT[dna_bytes]                        # 0..3, 4=non-ACGT
    comp = np.where(x == 4, 4, 3 - x).astype(np.uint8)   # complement (A<->T,C<->G)
    xr = comp[::-1]
    for src in (x, xr):
        n = src.size
        for frame in range(3):
            m = (n - frame) // 3
            if m < 10:
                continue
            cod = src[frame:frame + 3 * m].reshape(m, 3).astype(np.int64)
            bad = (cod == 4).any(axis=1)
            ci = cod[:, 0] * 16 + cod[:, 1] * 4 + cod[:, 2]
            ci[bad] = 0
            aa = AA_LUT[ci].copy()
            aa[bad] = ORD_STAR                     # non-ACGT codon breaks the ORF
            _scan_frame(aa, acc)

def load_intervals():
    """Return gene[chr] and repeat[chr] as lists of (start,end) 0-based half-open."""
    genes = {c: [] for c in CHRS}
    reps = {c: [] for c in CHRS}
    chrset = set(CHRS)
    with gzip.open(f"{SCRATCH}/dmel.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 5 or f[2] != "gene":
                continue
            c = f[0]
            if c in chrset:
                genes[c].append((int(f[3]) - 1, int(f[4])))     # GTF 1-based incl -> 0-based half-open
    with gzip.open(f"{SCRATCH}/dm6_rmsk.txt.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            gn = f[5]                                # genoName e.g. chr2L
            if gn.startswith("chr") and gn[3:] in chrset:
                reps[gn[3:]].append((int(f[6]), int(f[7])))     # UCSC 0-based start, end exclusive
    return genes, reps

def mask_from_intervals(intervals, L):
    m = np.zeros(L, bool)
    for s, e in intervals:
        if e <= 0 or s >= L:
            continue
        m[max(0, s):min(L, e)] = True
    return m

def runs_of_true(mask):
    """Yield (start,end) maximal True runs (0-based half-open)."""
    if mask.size == 0:
        return
    d = np.diff(mask.astype(np.int8))
    starts = list(np.flatnonzero(d == 1) + 1)
    ends = list(np.flatnonzero(d == -1) + 1)
    if mask[0]:
        starts = [0] + starts
    if mask[-1]:
        ends = ends + [mask.size]
    for s, e in zip(starts, ends):
        yield s, e

def build_intergenic():
    genes, reps = load_intervals()
    clean = Accum(); repint = Accum()
    for c in CHRS:
        seq = None
        for hdr, s in read_fasta(f"{SCRATCH}/dna_{c}.fa.gz"):
            seq = s.upper(); break
        Lc = len(seq)
        seqb = np.frombuffer(seq.encode(), np.uint8)
        gmask = mask_from_intervals(genes[c], Lc)
        rmask = mask_from_intervals(reps[c], Lc)
        clean_mask = (~gmask) & (~rmask)
        repint_mask = rmask & (~gmask)
        cbp = int(clean_mask.sum()); rbp = int(repint_mask.sum())
        clean.nbp += cbp; repint.nbp += rbp
        for s, e in runs_of_true(clean_mask):
            if e - s >= 30:
                scan_region_seq(seqb[s:e], clean)
        for s, e in runs_of_true(repint_mask):
            if e - s >= 30:
                scan_region_seq(seqb[s:e], repint)
        log(f"  chr{c}: len={Lc:,} clean_bp={cbp:,} repint_bp={rbp:,} "
            f"clean_orf>=2H={len(clean.rows):,} repint_orf>=2H={len(repint.rows):,}")
    return clean, repint

# =====================================================================
# bootstrap helpers for fold ratios
# =====================================================================
def boot_fold(tab, d, nboot=2000, rng=None):
    """Bootstrap distribution of pooled fold at spacing d by resampling rows."""
    if rng is None:
        rng = np.random.default_rng(SEED)
    n = tab.shape[0]
    if n == 0:
        return np.full(nboot, np.nan)
    o = tab[:, 1 + d].astype(float)
    e = analytic_exp_rows(tab, d)
    out = np.empty(nboot)
    for b in range(nboot):
        ix = rng.integers(0, n, n)
        se = e[ix].sum()
        out[b] = (o[ix].sum() / se) if se > 0 else np.nan
    return out

def boot_ci(tab, d, nboot=3000, seed=SEED):
    b = boot_fold(tab, d, nboot, np.random.default_rng(seed))
    b = b[~np.isnan(b)]
    if b.size == 0:
        return [np.nan, np.nan]
    return np.percentile(b, [2.5, 97.5]).tolist()

def stars(p):
    if p is None or np.isnan(p):
        return "ns"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"

def boot_diff_p(tabA, tabB, d, nboot=2000, seed=SEED):
    """Two-sided bootstrap p that foldA(d) != foldB(d) (independent resampling)."""
    rng = np.random.default_rng(seed)
    fa = boot_fold(tabA, d, nboot, rng)
    fb = boot_fold(tabB, d, nboot, rng)
    diff = fa - fb
    diff = diff[~np.isnan(diff)]
    if diff.size == 0:
        return np.nan, np.nan, np.nan
    p_gt = np.mean(diff <= 0); p_lt = np.mean(diff >= 0)
    p = 2 * min(p_gt, p_lt)
    p = min(1.0, max(p, 1.0 / (diff.size + 1)))   # floor at resolution
    return p, np.nanmean(fa), np.nanmean(fb)

# =====================================================================
# MAIN
# =====================================================================
def main():
    rng = np.random.default_rng(SEED)
    results = {"seed": SEED, "bin_labels": BIN_LABELS,
               "note": "folds are analytic scramble-null (== MC limit); "
                       "p-values are bootstrap over sequences (2-sided) unless stated"}

    # ---- load sets ----
    paths = sorted(glob.glob(f"{SCRATCH}/gc2023/*_DeNovoGenesNuc.fa"))
    den = load_denovo(paths, min_aa=10)
    prots_dn = [d["prot"] for d in den]
    log(f"de novo: {len(prots_dn)} proteins")
    meta = load_proteins(f"{SCRATCH}/dmel_metalbinding.fasta", min_aa=10)
    prot = load_proteins(f"{SCRATCH}/dmel_pep.fa.gz", min_aa=10)
    log(f"metal-binders {len(meta)}  proteome {len(prot)}")

    tab_dn = rows_from_prots(prots_dn)
    tab_meta = rows_from_prots(meta)
    tab_prot = rows_from_prots(prot)

    # ---- intergenic ORFs ----
    log("scanning intergenic ORFs ...")
    clean, repint = build_intergenic()
    tab_clean = clean.table(); tab_repint = repint.table()
    log(f"clean intergenic: {clean.nbp:,} bp, {clean.norf:,} ORFs>=10aa, "
        f"{tab_clean.shape[0]:,} with >=2 His")
    log(f"repeat intergenic: {repint.nbp:,} bp, {repint.norf:,} ORFs>=10aa, "
        f"{tab_repint.shape[0]:,} with >=2 His")

    # =====================================================================
    # FIGURE 1 : His-spacing profile d=1..6, 4 series + de novo vs clean stars
    # =====================================================================
    def profile(tab):
        return [pooled_fold(tab, d)[0] for d in range(1, DMAX + 1)]
    prof = {
        "denovo":   profile(tab_dn),
        "clean":    profile(tab_clean),
        "metal":    profile(tab_meta),
        "proteome": profile(tab_prot),
    }
    # bootstrap p: de novo vs clean per spacing
    p_prof = []
    for d in range(1, DMAX + 1):
        p, _, _ = boot_diff_p(tab_dn, tab_clean, d, nboot=3000, seed=SEED + d)
        p_prof.append(p)
    ci_dn_prof = [boot_ci(tab_dn, d, seed=SEED + 200 + d) for d in range(1, DMAX + 1)]
    ci_cl_prof = [boot_ci(tab_clean, d, seed=SEED + 300 + d) for d in range(1, DMAX + 1)]
    log("Fig1 profile de novo:", [round(x, 3) for x in prof["denovo"]])
    log("Fig1 profile clean  :", [round(x, 3) for x in prof["clean"]])
    log("Fig1 p(de novo vs clean):", [f"{x:.1e}" for x in p_prof])
    results["fig1"] = {"d": list(range(1, DMAX + 1)), "profiles": prof,
                       "ci_denovo": ci_dn_prof, "ci_clean": ci_cl_prof,
                       "p_denovo_vs_clean": p_prof,
                       "stars": [stars(p) for p in p_prof]}

    # =====================================================================
    # FIGURE 2 : length-binned HxH (d=2) fold, de novo / clean / repint
    # =====================================================================
    def binned_fold(tab, d=2, ci_seed=None):
        L = tab[:, 0]
        vals, ns, cis = [], [], []
        for i in range(len(BIN_EDGES) - 1):
            lo, hi = BIN_EDGES[i], BIN_EDGES[i + 1]
            sub = tab[(L >= lo) & (L < hi)]
            f, o, e = pooled_fold(sub, d)
            vals.append(f); ns.append(int(sub.shape[0]))
            cis.append(boot_ci(sub, d, seed=(ci_seed or SEED) + i) if sub.shape[0] else [np.nan, np.nan])
        return vals, ns, cis
    f_dn, n_dn, ci_dn2 = binned_fold(tab_dn, ci_seed=SEED + 400)
    f_cl, n_cl, ci_cl2 = binned_fold(tab_clean, ci_seed=SEED + 500)
    f_rp, n_rp, ci_rp2 = binned_fold(tab_repint, ci_seed=SEED + 600)
    # bootstrap p de novo vs clean per bin
    p_bin = []
    L_dn, L_cl = tab_dn[:, 0], tab_clean[:, 0]
    for i in range(len(BIN_EDGES) - 1):
        lo, hi = BIN_EDGES[i], BIN_EDGES[i + 1]
        A = tab_dn[(L_dn >= lo) & (L_dn < hi)]
        B = tab_clean[(L_cl >= lo) & (L_cl < hi)]
        p, _, _ = boot_diff_p(A, B, 2, nboot=3000, seed=SEED + 100 + i)
        p_bin.append(p)
    log("Fig2 de novo binned:", [round(x, 3) for x in f_dn], "target [1.31,1.76,1.79,1.97]")
    log("Fig2 clean  binned:", [round(x, 3) for x in f_cl], "target [1.41,1.40,1.60,1.65]")
    log("Fig2 repint binned:", [round(x, 3) for x in f_rp], "target [1.63,1.21,1.18,0.99]")
    log("Fig2 p(de novo vs clean):", [f"{x:.1e}" for x in p_bin])
    results["fig2"] = {"bins": BIN_LABELS,
                       "denovo": f_dn, "clean": f_cl, "repint": f_rp,
                       "ci_denovo": ci_dn2, "ci_clean": ci_cl2, "ci_repint": ci_rp2,
                       "n_denovo": n_dn, "n_clean": n_cl, "n_repint": n_rp,
                       "p_denovo_vs_clean": p_bin, "stars": [stars(p) for p in p_bin]}

    # =====================================================================
    # FIGURE 3 : CA-repeat mechanism (de novo only, needs coding DNA)
    # =====================================================================
    mid_counts = {}          # middle residue of H-x-H motifs
    mid_codon = {}           # middle-residue codon
    hxh_his_codon = {"CAC": 0, "CAT": 0}   # His codons at motif (H-x-H) positions (unique)
    non_his_codon = {"CAC": 0, "CAT": 0}   # His codons at non-motif positions (unique)
    all_his_codon = {"CAC": 0, "CAT": 0}   # His codons overall (de novo)
    n_hxh = 0
    aa_bg = np.zeros(256, np.int64)  # de novo aa background
    for rec in den:
        p = rec["prot"]; dna = rec["dna"]
        aa = np.frombuffer(p.encode(), np.uint8)
        aa_bg += np.bincount(aa, minlength=256)
        hmask = aa == ORD_H
        L = aa.size
        # His codons: overall, and classified by motif membership at UNIQUE positions.
        # A His at i is a "motif" His if it participates in >=1 H-x-H (an H at i-2 or i+2);
        # otherwise "non-motif". Each His is counted once -> removes the double-counting of
        # the central H of an H-x-H-x-H run that motivated the earlier permutation test.
        for i in np.flatnonzero(hmask):
            cod = dna[3 * i:3 * i + 3]
            if cod not in all_his_codon:
                continue
            all_his_codon[cod] += 1
            in_motif = (i - 2 >= 0 and hmask[i - 2]) or (i + 2 < L and hmask[i + 2])
            if in_motif:
                hxh_his_codon[cod] += 1
            else:
                non_his_codon[cod] += 1
        # H-x-H motifs (middle-residue tallies; motif instance count)
        if L > 2:
            starts = np.flatnonzero(hmask[:-2] & hmask[2:])
            for i in starts:
                n_hxh += 1
                mid = p[i + 1]
                mid_counts[mid] = mid_counts.get(mid, 0) + 1
                mc = dna[3 * (i + 1):3 * (i + 1) + 3]
                mid_codon[mc] = mid_codon.get(mc, 0) + 1
    tot_bg = aa_bg[[ord(a) for a in AAS]].sum()
    bg_freq = {a: aa_bg[ord(a)] / tot_bg for a in AAS}
    mid_freq = {a: mid_counts.get(a, 0) / n_hxh for a in AAS}
    # Thr enrichment as middle residue: permutation/bootstrap p
    thr_mid = mid_counts.get("T", 0)
    thr_bg = bg_freq["T"]
    # bootstrap: resample H-x-H motifs, distribution of Thr-middle fraction; p vs bg
    mids = np.array([1 if m == "T" else 0 for m in
                     sum([[k] * v for k, v in mid_counts.items()], [])], np.int8)
    rng2 = np.random.default_rng(SEED + 7)
    nb = 5000
    boot_thr = np.empty(nb)
    for b in range(nb):
        ix = rng2.integers(0, mids.size, mids.size)
        boot_thr[b] = mids[ix].mean()
    # one-sided p that Thr-middle fraction > background
    p_thr = np.mean(boot_thr <= thr_bg)
    p_thr = max(p_thr, 1.0 / (nb + 1)) * 2  # two-sided
    p_thr = min(1.0, p_thr)
    # CAC vs CAT bias: motif (H-x-H) His vs NON-motif His, unique positions -> Fisher exact 2x2.
    # (Reviewer 2 asked for the apples-to-apples motif-vs-non-motif contrast and a Fisher's/binomial
    #  test in place of the permutation shuffle; unique positions remove the double-count issue.)
    a1, b1 = hxh_his_codon["CAC"], hxh_his_codon["CAT"]   # motif His
    a2, b2 = non_his_codon["CAC"], non_his_codon["CAT"]   # non-motif His
    a0, b0 = all_his_codon["CAC"], all_his_codon["CAT"]   # all His (kept for the record / sanity)
    frac_hxh = a1 / (a1 + b1)
    frac_non = a2 / (a2 + b2)
    frac_all = a0 / (a0 + b0)
    from scipy.stats import fisher_exact
    _or_cac, p_cac = fisher_exact([[a1, b1], [a2, b2]], alternative="two-sided")
    top_mid_codons = sorted(mid_codon.items(), key=lambda kv: -kv[1])[:6]
    log(f"Fig3 n_HxH={n_hxh}  Thr-middle={thr_mid} ({100*mid_freq['T']:.1f}%) "
        f"bg Thr {100*thr_bg:.1f}%  fold {mid_freq['T']/thr_bg:.2f}x  p_thr={p_thr:.1e}")
    log(f"Fig3 CAC frac: H-x-H His {100*frac_hxh:.1f}% vs non-H-x-H His {100*frac_non:.1f}% "
        f"(all His {100*frac_all:.1f}%)  Fisher p={p_cac:.2e}")
    log(f"Fig3 top middle codons: {top_mid_codons}")
    results["fig3"] = {
        "n_hxh": n_hxh,
        "mid_freq": mid_freq, "bg_freq": bg_freq,
        "thr_mid_frac": mid_freq["T"], "thr_bg_frac": thr_bg,
        "thr_fold": mid_freq["T"] / thr_bg, "p_thr": float(p_thr), "stars_thr": stars(p_thr),
        "cac_frac_hxh": frac_hxh, "cac_frac_non": frac_non, "cac_frac_all": frac_all,
        "cat_frac_hxh": 1 - frac_hxh, "cat_frac_non": 1 - frac_non, "cat_frac_all": 1 - frac_all,
        "p_cac": float(p_cac), "stars_cac": stars(p_cac),
        "top_mid_codons": top_mid_codons,
    }

    # =====================================================================
    # FIGURE 4 : dinucleotide-preserving null (decisive)
    # =====================================================================
    log("dinucleotide-null: AE shuffles per de novo CDS ...")
    NDN = 25
    rng4 = np.random.default_rng(SEED + 11)
    o_i = np.zeros(len(den))          # observed HxH per seq
    e_aa = np.zeros(len(den))         # aa-scramble expected (analytic)
    e_dn = np.zeros(len(den))         # dinuc-null expected (mean over AE shuffles)
    skipped = 0
    for k, rec in enumerate(den):
        p = rec["prot"]; dna = rec["dna"]
        aa = np.frombuffer(p.encode(), np.uint8)
        hmask = aa == ORD_H
        L = aa.size; h = int(hmask.sum())
        if L > 2 and h >= 2:
            o_i[k] = int(np.count_nonzero(hmask[:-2] & hmask[2:]))
            e_aa[k] = (L - 2) * h * (h - 1) / (L * (L - 1))
        # dinuc null: AE-shuffle DNA, translate full (positions fixed), count HxH
        if len(dna) >= 6 and set(dna) <= set("ACGT"):
            acc = 0.0
            for _ in range(NDN):
                sh = altschul_erikson_shuffle(dna, rng4)
                pep = translate_full(sh)
                a2 = np.frombuffer(pep.encode(), np.uint8)
                hm2 = a2 == ORD_H
                if a2.size > 2:
                    acc += int(np.count_nonzero(hm2[:-2] & hm2[2:]))
            e_dn[k] = acc / NDN
        else:
            skipped += 1
            e_dn[k] = np.nan
        if (k + 1) % 2000 == 0:
            log(f"   dinuc {k+1}/{len(den)}")
    valid = ~np.isnan(e_dn)
    fold_aa = o_i.sum() / e_aa.sum()
    fold_dn = o_i[valid].sum() / e_dn[valid].sum()
    # bootstrap over sequences (valid set) for both folds + p vs 1 and vs each other
    idxv = np.flatnonzero(valid)
    ov, eav, edv = o_i[idxv], e_aa[idxv], e_dn[idxv]
    nb3 = 4000
    rng5 = np.random.default_rng(SEED + 12)
    bf_aa = np.empty(nb3); bf_dn = np.empty(nb3)
    n = idxv.size
    for b in range(nb3):
        ii = rng5.integers(0, n, n)
        bf_aa[b] = ov[ii].sum() / eav[ii].sum()
        bf_dn[b] = ov[ii].sum() / edv[ii].sum()
    ci_aa = np.percentile(bf_aa, [2.5, 97.5])
    ci_dn = np.percentile(bf_dn, [2.5, 97.5])
    p_aa1 = 2 * min(np.mean(bf_aa <= 1), np.mean(bf_aa >= 1)); p_aa1 = min(1, max(p_aa1, 1/(nb3+1)))
    p_dn1 = 2 * min(np.mean(bf_dn <= 1), np.mean(bf_dn >= 1)); p_dn1 = min(1, max(p_dn1, 1/(nb3+1)))
    dd = bf_aa - bf_dn
    p_aavdn = 2 * min(np.mean(dd <= 0), np.mean(dd >= 0)); p_aavdn = min(1, max(p_aavdn, 1/(nb3+1)))
    # fraction of excess explained by dinucleotide composition
    frac_explained = (fold_aa - fold_dn) / (fold_aa - 1.0)
    verdict = "SURVIVES" if (fold_dn > 1.0 and p_dn1 < 0.05) else "COLLAPSES"
    log(f"Fig4 fold_aa={fold_aa:.3f} CI{ci_aa.round(3)}  p_vs1={p_aa1:.1e}")
    log(f"Fig4 fold_dn={fold_dn:.3f} CI{ci_dn.round(3)}  p_vs1={p_dn1:.1e}  (skipped {skipped})")
    log(f"Fig4 p(aa vs dinuc)={p_aavdn:.1e}  frac of excess explained by dinuc={frac_explained:.2f}")
    log(f"Fig4 VERDICT: de novo H-x-H signal {verdict} the dinucleotide-preserving null")
    results["fig4"] = {
        "NDN": NDN, "n_valid": int(n), "skipped": int(skipped),
        "fold_aa": float(fold_aa), "ci_aa": ci_aa.tolist(), "p_aa_vs1": float(p_aa1),
        "fold_dn": float(fold_dn), "ci_dn": ci_dn.tolist(), "p_dn_vs1": float(p_dn1),
        "p_aa_vs_dn": float(p_aavdn), "frac_excess_explained_by_dinuc": float(frac_explained),
        "verdict": verdict,
        "stars_aa_vs1": stars(p_aa1), "stars_dn_vs1": stars(p_dn1),
        "stars_aa_vs_dn": stars(p_aavdn),
    }

    # =====================================================================
    # FIGURE 5 : amino-acid composition by class
    # =====================================================================
    # random ORFs at intergenic GC=0.413, length-matched to de novo
    def gen_random_orfs(lengths, gc=0.413, seed=SEED):
        rr = np.random.default_rng(seed)
        pG = pC = gc / 2; pA = pT = (1 - gc) / 2
        probs = np.array([pA, pC, pG, pT])   # A C G T
        bases = np.array(list("ACGT"))
        comp = np.zeros(256, np.int64)
        stops = {"TAA", "TAG", "TGA"}
        for L in lengths:
            got = []
            while len(got) < L:
                need = (L - len(got)) * 2 + 8
                draw = rr.choice(4, size=(need, 3), p=probs)
                for cod in draw:
                    c = bases[cod[0]] + bases[cod[1]] + bases[cod[2]]
                    if c in stops:
                        break
                    got.append(CODON2AA[c])
                    if len(got) >= L:
                        break
            aa = np.frombuffer("".join(got[:L]).encode(), np.uint8)
            comp += np.bincount(aa, minlength=256)
        return comp
    Ldn = [len(p) for p in prots_dn]
    comp_rand = gen_random_orfs(Ldn, gc=0.413, seed=SEED + 21)

    def comp_pct(hist):
        tot = sum(hist[ord(a)] for a in AAS)
        return {a: 100.0 * hist[ord(a)] / tot for a in AAS}
    comp_dn = np.zeros(256, np.int64)
    for p in prots_dn:
        comp_dn += np.bincount(np.frombuffer(p.encode(), np.uint8), minlength=256)
    comp_prot = np.zeros(256, np.int64)
    for p in prot:
        comp_prot += np.bincount(np.frombuffer(p.encode(), np.uint8), minlength=256)
    comp_meta = np.zeros(256, np.int64)
    for p in meta:
        comp_meta += np.bincount(np.frombuffer(p.encode(), np.uint8), minlength=256)
    comp_hist = {
        "de novo":       comp_dn,
        "intergenic":    clean.comp,
        "random ORF":    comp_rand,
        "proteome":      comp_prot,
        "metal-binders": comp_meta,
    }
    comp_sets = {k: comp_pct(v) for k, v in comp_hist.items()}
    comp_tot = {k: int(sum(v[ord(a)] for a in AAS)) for k, v in comp_hist.items()}
    log("Fig5 His%: " + "  ".join(f"{k}={comp_sets[k]['H']:.2f}" for k in comp_sets))
    log("Fig5 Cys%: " + "  ".join(f"{k}={comp_sets[k]['C']:.2f}" for k in comp_sets))
    # two-proportion z-test (real p) for highlighted residues: de novo vs proteome
    from math import sqrt, erfc
    def two_prop_p(x1, n1, x2, n2):
        p1, p2 = x1 / n1, x2 / n2
        p = (x1 + x2) / (n1 + n2)
        se = sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
        if se == 0:
            return 0.0, 0.0
        z = (p1 - p2) / se
        return z, erfc(abs(z) / sqrt(2))   # two-sided
    p_comp = {}
    for a in ["H", "C", "D", "E"]:
        z, pv = two_prop_p(int(comp_hist["de novo"][ord(a)]), comp_tot["de novo"],
                           int(comp_hist["proteome"][ord(a)]), comp_tot["proteome"])
        p_comp[a] = {"z": z, "p": pv, "star": stars(pv)}
    log("Fig5 de novo vs proteome (H,C,D,E) z/p: " +
        "  ".join(f"{a}:z={p_comp[a]['z']:.1f},p={p_comp[a]['p']:.1e}" for a in "HCDE"))
    results["fig5"] = {"aa_order": AAS, "sets": comp_sets, "totals": comp_tot,
                       "highlight": ["H", "C", "D", "E"],
                       "p_denovo_vs_proteome": p_comp}

    # ---- write outputs ----
    with open(f"{OUT}/hisspacing_results.json", "w") as fh:
        json.dump(results, fh, indent=2)
    log("wrote " + f"{OUT}/hisspacing_results.json")

    # tidy CSVs for R
    import csv
    with open(f"{OUT}/fig1_profile.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["d", "denovo", "clean", "metal", "proteome",
                                        "dn_lo", "dn_hi", "cl_lo", "cl_hi", "p_dn_clean", "star"])
        for i, d in enumerate(range(1, DMAX + 1)):
            w.writerow([d, prof["denovo"][i], prof["clean"][i], prof["metal"][i], prof["proteome"][i],
                        ci_dn_prof[i][0], ci_dn_prof[i][1], ci_cl_prof[i][0], ci_cl_prof[i][1],
                        p_prof[i], stars(p_prof[i])])
    with open(f"{OUT}/fig2_binned.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["bin", "denovo", "clean", "repint",
                                        "dn_lo", "dn_hi", "cl_lo", "cl_hi", "rp_lo", "rp_hi",
                                        "n_denovo", "n_clean", "n_repint", "p_dn_clean", "star"])
        for i, lab in enumerate(BIN_LABELS):
            w.writerow([lab, f_dn[i], f_cl[i], f_rp[i],
                        ci_dn2[i][0], ci_dn2[i][1], ci_cl2[i][0], ci_cl2[i][1], ci_rp2[i][0], ci_rp2[i][1],
                        n_dn[i], n_cl[i], n_rp[i], p_bin[i], stars(p_bin[i])])
    with open(f"{OUT}/fig2_middle.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["aa", "mid_freq_pct", "bg_freq_pct"])
        for a in AAS:
            w.writerow([a, 100 * mid_freq[a], 100 * bg_freq[a]])
    with open(f"{OUT}/fig2_codon.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["group", "CAC_pct", "CAT_pct"])
        w.writerow(["H-x-H His", 100 * frac_hxh, 100 * (1 - frac_hxh)])
        w.writerow(["non-H-x-H His", 100 * frac_non, 100 * (1 - frac_non)])
    with open(f"{OUT}/fig2_stats.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["key", "value"])
        for k, v in [("n_hxh", n_hxh),
                     ("thr_mid_pct", 100 * mid_freq["T"]), ("thr_bg_pct", 100 * thr_bg),
                     ("thr_fold", mid_freq["T"] / thr_bg),
                     ("p_thr", p_thr), ("star_thr", stars(p_thr)),
                     ("cac_hxh_pct", 100 * frac_hxh), ("cac_non_pct", 100 * frac_non),
                     ("p_cac", p_cac), ("star_cac", stars(p_cac)),
                     ("top_mid_codon", top_mid_codons[0][0]),
                     ("top_mid_codon_n", top_mid_codons[0][1])]:
            w.writerow([k, v])
    with open(f"{OUT}/fig4_dinuc.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["null", "fold", "lo", "hi", "p_vs1", "star_vs1"])
        w.writerow(["aa-scramble", fold_aa, ci_aa[0], ci_aa[1], p_aa1, stars(p_aa1)])
        w.writerow(["dinuc-preserving", fold_dn, ci_dn[0], ci_dn[1], p_dn1, stars(p_dn1)])
    with open(f"{OUT}/figS2_comp.csv", "w", newline="") as fh:
        setnames = list(comp_sets.keys())
        w = csv.writer(fh); w.writerow(["aa"] + setnames)
        for a in AAS:
            w.writerow([a] + [comp_sets[s][a] for s in setnames])
    with open(f"{OUT}/figS2_sig.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["aa", "denovo_pct", "proteome_pct", "z", "p", "star"])
        for a in ["H", "C", "D", "E"]:
            w.writerow([a, comp_sets["de novo"][a], comp_sets["proteome"][a],
                        p_comp[a]["z"], p_comp[a]["p"], p_comp[a]["star"]])
    log("all CSVs written. DONE.")

if __name__ == "__main__":
    main()
