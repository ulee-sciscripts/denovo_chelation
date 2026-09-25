#!/usr/bin/env python3
"""H-x-H incidence by length quartile, for intergenic ORFs, neORFs, and the
length-matched proteome, with an amino-acid-shuffle expectation alongside.

Motivation: raw H-x-H *incidence* (fraction of ORFs with >=1 motif) is strongly
length-dependent -- a 3-residue motif (His at i, i+2) is more likely present the more
positions an ORF has -- while His *prevalence* (% of residues that are His) is roughly
flat with length. This resolves the apparent 3.9%-vs-10.95% intergenic discrepancy
(pooled-raw vs neORF-length-matched) and isolates ordering from length+composition.

For each of the three sequence sets we split the ORFs into four equal-count quartiles
by length and report, per quartile:
  - observed H-x-H incidence (fraction with >=1 H-x-H)
  - amino-acid-shuffle expected incidence: each ORF's own residues placed at uniformly
    random positions (length- and composition-preserving), P(>=1 H-x-H); Monte-Carlo per
    unique (length, His-count), B=2000
  - obs / shuffle-expected (the incidence-level ordering excess)
  - His prevalence (% of residues that are His)

Run from code/cxc_analysis/. External inputs live under ./data/ and ../../data/analysis/
(see README.md); the intergenic scan is cached to ./data/interg_LNH.npy. Output:
results/length_quartile_hxh.csv (committed) + a printed table.
"""
import os, sys, glob, csv
import numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, ".."))
import denovo_hisspacing_analysis as A
from hs_common import load_denovo, load_proteins, read_fasta

SCRATCH  = f"{_HERE}/data/md_hp6_scratch"                 # external (git-ignored)
GENOME   = f"{_HERE}/../../data/analysis/ref/dm6_genome.fa"
INTERVAL = f"{_HERE}/data/interg_intervals"               # dmel.gtf.gz + dm6_rmsk.txt.gz
CACHE    = f"{_HERE}/data/interg_LNH.npy"
CHRS     = ["2L", "2R", "3L", "3R"]
SEED_MATCH, SEED_MC, B = 20260812, 20260925, 2000
ORD_H = ord("H")


def seq_arrays(seqs):
    """(L, nH, hxh) per sequence."""
    L  = np.array([len(s) for s in seqs], float)
    nH = np.array([s.count("H") for s in seqs], float)
    hx = np.array([sum(1 for i in range(len(s) - 2) if s[i] == "H" and s[i + 2] == "H")
                   for s in seqs], float)
    return L, nH, hx


def shuffle_expected_incidence(L, nH, b=B, seed=SEED_MC):
    """P(>=1 H-x-H) under a composition- and length-preserving amino-acid shuffle,
    Monte-Carlo per unique (length, His-count), mapped back to each ORF."""
    rng = np.random.default_rng(seed)
    Li, hi = L.astype(int), nH.astype(int)
    pairs, inv = np.unique(np.stack([Li, hi], 1), axis=0, return_inverse=True)
    p = np.zeros(len(pairs))
    for j, (Lp, hp) in enumerate(pairs):
        if hp < 2 or Lp < 3:
            continue
        keys = rng.random((b, Lp))
        thr = np.partition(keys, hp - 1, axis=1)[:, hp - 1][:, None]
        mask = keys <= thr                       # the hp smallest per row -> True
        p[j] = (mask[:, :-2] & mask[:, 2:]).any(axis=1).mean()
    return p[inv]


def quartile_rows(name, L, nH, hx):
    p = shuffle_expected_incidence(L, nH)
    n = L.size
    order = np.argsort(L, kind="stable")
    Ls, nHs, hxs, ps = L[order], nH[order], hx[order], p[order]
    rows = []
    for q, g in enumerate(np.array_split(np.arange(n), 4), 1):
        lg, hg, xg, pg = Ls[g], nHs[g], hxs[g], ps[g]
        obs = float(np.mean(xg >= 1)); exp = float(pg.mean())
        rows.append(dict(set=name, quartile=f"Q{q}",
                         len_min=int(lg.min()), len_max=int(lg.max()),
                         n=int(g.size), median_len=int(np.median(lg)),
                         obs_incidence_pct=round(100 * obs, 3),
                         shuffle_exp_pct=round(100 * exp, 3),
                         obs_over_exp=round(obs / exp, 3) if exp else float("nan"),
                         his_prev_pct=round(100 * hg.sum() / lg.sum(), 3)))
    return rows


def load_intergenic():
    if os.path.exists(CACHE):
        a = np.load(CACHE); return a[0].astype(float), a[1].astype(float), a[2].astype(float)
    A.SCRATCH = INTERVAL
    B2I, LUT, STAR = A.BASE2INT, A.AA_LUT, A.ORD_STAR
    seqs = {}
    for hdr, s in read_fasta(GENOME):
        c = hdr.split()[0]
        if c in CHRS:
            seqs[c] = s.upper()
    genes, reps = A.load_intervals()

    def scan(seqb):
        out = []; x = B2I[seqb]; xr = np.where(x == 4, 4, 3 - x).astype(np.uint8)[::-1]
        for src in (x, xr):
            n = src.size
            for fr in range(3):
                m = (n - fr) // 3
                if m < 10:
                    continue
                cod = src[fr:fr + 3 * m].reshape(m, 3).astype(np.int64)
                v = (cod < 4).all(axis=1)
                idx = np.where(v, cod[:, 0] * 16 + cod[:, 1] * 4 + cod[:, 2], 0)
                aa = np.where(v, LUT[idx], STAR).astype(np.uint8)
                st = np.flatnonzero(aa == STAR); bnd = np.concatenate(([-1], st, [aa.size]))
                for k in range(bnd.size - 1):
                    s0, e0 = bnd[k] + 1, bnd[k + 1]; Ln = e0 - s0
                    if Ln < 10:
                        continue
                    hm = aa[s0:e0] == ORD_H
                    out.append((Ln, int(hm.sum()),
                                int((hm[:-2] & hm[2:]).sum()) if Ln > 2 else 0))
        return out

    orfs = []
    for c in CHRS:
        sb = np.frombuffer(seqs[c].encode(), np.uint8); Lc = len(seqs[c])
        clean = (~A.mask_from_intervals(genes[c], Lc)) & (~A.mask_from_intervals(reps[c], Lc))
        for s, e in A.runs_of_true(clean):
            if e - s >= 30:
                orfs.extend(scan(sb[s:e]))
    L  = np.array([o[0] for o in orfs], float)
    nH = np.array([o[1] for o in orfs], float)
    hx = np.array([o[2] for o in orfs], float)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.save(CACHE, np.vstack([L, nH, hx]))
    return L, nH, hx


def matched_proteome(dn_seqs, prot):
    dl = np.array([len(s) for s in dn_seqs]); pl = np.array([len(s) for s in prot])
    o = np.argsort(pl, kind="stable"); ps = [prot[i] for i in o]; pls = pl[o]; nP = len(ps)
    rng = np.random.default_rng(SEED_MATCH); out = []
    for L in dl:
        j = int(np.searchsorted(pls, int(L)))
        cand = [k for k in (j - 1, j, j + 1) if 0 <= k < nP]
        d = np.array([abs(int(pls[k]) - int(L)) for k in cand])
        best = [cand[t] for t in range(len(cand)) if d[t] == d.min()]
        out.append(ps[int(rng.choice(best))])
    return out


def main():
    dn = list(dict.fromkeys(d["prot"] for d in
              load_denovo(sorted(glob.glob(f"{SCRATCH}/gc2023/*_DeNovoGenesNuc.fa")), min_aa=10)))
    prot = load_proteins(f"{SCRATCH}/data/dmel_pep.fa.gz", min_aa=10)
    sets = [("intergenic", *load_intergenic()),
            ("neORF", *seq_arrays(dn)),
            ("length-matched proteome", *seq_arrays(matched_proteome(dn, prot)))]
    rows = []
    for name, L, nH, hx in sets:
        rows += quartile_rows(name, L, nH, hx)
    os.makedirs(f"{_HERE}/results", exist_ok=True)
    out = f"{_HERE}/results/length_quartile_hxh.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    hdr = f"{'set':>24s} {'Q':>3s} {'len(aa)':>10s} {'n':>9s} {'med':>4s} {'obs>=1HxH':>10s} {'shuf-exp':>9s} {'obs/exp':>8s} {'His%':>6s}"
    print(hdr)
    for r in rows:
        lr = "%d-%d" % (r["len_min"], r["len_max"])
        print(f"{r['set']:>24s} {r['quartile']:>3s} {lr:>10s} "
              f"{r['n']:>9,} {r['median_len']:>4d} {r['obs_incidence_pct']:>9.2f}% "
              f"{r['shuffle_exp_pct']:>8.2f}% {r['obs_over_exp']:>8.2f} {r['his_prev_pct']:>5.2f}%")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
