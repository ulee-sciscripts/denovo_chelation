#!/usr/bin/env python3
"""
Regenerate the two comparison sets that were not yet saved as tidy CSVs
(FIGURE_WORK.md sec.5), reusing the EXACT fold/bootstrap code paths of
denovo_hisspacing_analysis.py so everything stays consistent with fig1/fig2.

Produces, under docs/results/data/hisspacing/:
  fig_lenmatched_profile.csv   -- Fig 5: de novo vs length-matched proteome,
                                  His-pair order-fold profile d=1..6 (+CIs, p).
  fig2_hxh.csv / fig2_cxxc.csv -- Fig 2: HxH and CxxC order-fold points (+CIs)
                                  for the chelation-grammar spectrum.
  fig1_stats.csv               -- Figure 1: bootstrap/permutation p for each bracket.

Length-matching method (FIGURE_WORK.md sec.5): for each de novo neORF, match to
the nearest-length WHOLE protein of the comparator (with replacement) -- identical
to denovo_chelation_analysis.py 'whole_protein_match'. Verification anchors printed.
"""
import os, sys, csv, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hs_common import SCRATCH, load_denovo, load_proteins, count_spacings
import denovo_hisspacing_analysis as H   # reuse pooled_fold / boot_ci / boot_diff_p / stars

OUT = os.path.join(HERE, os.pardir, "data", "hisspacing")
CHEL = os.path.join(HERE, os.pardir, "data", "intermediate")
SEED = H.SEED
DMAX = 6


def rows_for_res(prots, res, min_aa=10):
    """[L, count(res), o1..o6] per protein with >=2 of res and L>=min_aa.
    Same layout pooled_fold/boot_ci expect (col0=L, col1=count, col1+d=obs_d)."""
    rows = []
    for p in prots:
        L = len(p); c = p.count(res)
        if c >= 2 and L >= min_aa:
            o = count_spacings(p, DMAX, res)     # d=1..6
            rows.append((L, c, *[int(x) for x in o]))
    return np.array(rows, dtype=np.int64) if rows else np.zeros((0, 8), np.int64)


def tab_from_arrays(L, cnt, o_d, d):
    """Minimal single-spacing table from per-protein arrays (for CSV-sourced sets)."""
    t = np.zeros((len(L), 8), np.int64)
    t[:, 0] = L; t[:, 1] = cnt; t[:, 1 + d] = o_d
    return t


def read_percluster(path):
    L, Hc, Cc, hxh, cxxc = [], [], [], [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            L.append(int(row["length"])); Hc.append(int(row["His_count"]))
            Cc.append(int(row["Cys_count"])); hxh.append(int(row["HxH"]))
            cxxc.append(int(row["CxxC"]))
    return (np.array(L), np.array(Hc), np.array(Cc), np.array(hxh), np.array(cxxc))


def read_fig1_intergenic_d2():
    """Pull the clean-intergenic HxH (d=2) fold, CI and de-novo-vs-clean p from
    the committed fig1_profile.csv (do not recompute the whole genome scan)."""
    with open(os.path.join(OUT, "fig1_profile.csv")) as fh:
        for row in csv.DictReader(fh):
            if int(row["d"]) == 2:
                return (float(row["clean"]), float(row["cl_lo"]), float(row["cl_hi"]),
                        float(row["p_dn_clean"]), row["star"])
    raise RuntimeError("d=2 row not found in fig1_profile.csv")


def main():
    print("loading sets ...", flush=True)
    paths = sorted(glob.glob(f"{SCRATCH}/gc2023/*_DeNovoGenesNuc.fa"))
    den = load_denovo(paths, min_aa=10)
    prots_dn = [d["prot"] for d in den]
    prot = load_proteins(f"{SCRATCH}/dmel_pep.fa.gz", min_aa=10)
    meta = load_proteins(f"{SCRATCH}/dmel_metalbinding.fasta", min_aa=10)
    print(f"  de novo={len(prots_dn)}  proteome={len(prot)}  metal={len(meta)}")

    # ---------- de novo His rows (must reproduce fig1 denovo) ----------
    tab_dn_H = rows_for_res(prots_dn, "H")
    prof_dn = [H.pooled_fold(tab_dn_H, d)[0] for d in range(1, DMAX + 1)]
    print("ANCHOR de novo H profile:", [round(x, 3) for x in prof_dn],
          " target fig1 [0.902,1.764,1.085,1.418,1.120,1.294]")

    # ---------- whole-protein nearest-length matched proteome ----------
    dn_lens = np.array([len(p) for p in prots_dn])
    plen = np.array([len(s) for s in prot])
    order = np.argsort(plen, kind="stable")
    psort = [prot[i] for i in order]
    plsort = plen[order]
    rngm = np.random.default_rng(SEED)
    matched, nP = [], len(psort)
    for L in dn_lens:
        j = int(np.searchsorted(plsort, int(L)))
        cand = [k for k in (j - 1, j, j + 1) if 0 <= k < nP]
        diffs = np.array([abs(int(plsort[k]) - int(L)) for k in cand])
        best = [cand[t] for t in range(len(cand)) if diffs[t] == diffs.min()]
        matched.append(psort[int(rngm.choice(best))])
    tab_mt_H = rows_for_res(matched, "H")
    prof_mt = [H.pooled_fold(tab_mt_H, d)[0] for d in range(1, DMAX + 1)]
    print("ANCHOR matched H profile:", [round(x, 3) for x in prof_mt],
          " (d=2 must reproduce whole_protein_match 0.996)")

    ci_dn = [H.boot_ci(tab_dn_H, d, seed=SEED + 200 + d) for d in range(1, DMAX + 1)]
    ci_mt = [H.boot_ci(tab_mt_H, d, seed=SEED + 700 + d) for d in range(1, DMAX + 1)]
    p_dm = [H.boot_diff_p(tab_dn_H, tab_mt_H, d, nboot=3000, seed=SEED + 800 + d)[0]
            for d in range(1, DMAX + 1)]

    with open(os.path.join(OUT, "fig_lenmatched_profile.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["d", "denovo", "matched", "dn_lo", "dn_hi", "mt_lo", "mt_hi",
                    "p_dn_matched", "star"])
        for i, d in enumerate(range(1, DMAX + 1)):
            w.writerow([d, prof_dn[i], prof_mt[i], ci_dn[i][0], ci_dn[i][1],
                        ci_mt[i][0], ci_mt[i][1], p_dm[i], H.stars(p_dm[i])])
    print("wrote fig_lenmatched_profile.csv")

    # ---------- Fig 2 spectrum points ----------
    tab_meta_H = rows_for_res(meta, "H")
    tab_dn_C = rows_for_res(prots_dn, "C")
    tab_meta_C = rows_for_res(meta, "C")
    tab_mt_C = rows_for_res(matched, "C")

    fdn_h2 = H.pooled_fold(tab_dn_H, 2)[0]
    fmt_h2 = H.pooled_fold(tab_meta_H, 2)[0]
    fdn_c3 = H.pooled_fold(tab_dn_C, 3)[0]
    fmt_c3 = H.pooled_fold(tab_meta_C, 3)[0]
    fmatch_c3 = H.pooled_fold(tab_mt_C, 3)[0]
    print(f"ANCHOR HxH de novo/metal = {fdn_h2:.4f}/{fmt_h2:.4f}"
          f"  (spectrum 1.7637/1.7663)")
    print(f"ANCHOR CxxC de novo/metal/matched = {fdn_c3:.4f}/{fmt_c3:.4f}/{fmatch_c3:.4f}"
          f"  (spectrum 1.2421/4.8634; whole_protein_match 1.148)")

    # random ORF background from the committed per-cluster CSV
    Lr, Hr, Cr, hxhr, cxxcr = read_percluster(f"{CHEL}/chelation_percluster_randorf.csv")
    tab_rand_H = tab_from_arrays(Lr, Hr, hxhr, 2)
    tab_rand_C = tab_from_arrays(Lr, Cr, cxxcr, 3)
    fr_h2 = H.pooled_fold(tab_rand_H, 2)[0]
    fr_c3 = H.pooled_fold(tab_rand_C, 3)[0]
    print(f"ANCHOR randORF HxH/CxxC = {fr_h2:.4f}/{fr_c3:.4f}  (~1.03/1.09)")

    # intergenic (repeat-masked) HxH from committed fig1
    ig_f, ig_lo, ig_hi, p_dn_ig, star_ig = read_fig1_intergenic_d2()

    # CIs for spectrum points (bootstrap over each set's rows)
    ci = lambda tab, d, s: H.boot_ci(tab, d, seed=SEED + s)
    hxh_rows = [
        ("random ORF",              fr_h2, *ci(tab_rand_H, 2, 11)),
        ("intergenic (rpt-masked)", ig_f,  ig_lo, ig_hi),
        ("de novo neORF",           fdn_h2, *ci(tab_dn_H, 2, 12)),
        ("metal-binders",           fmt_h2, *ci(tab_meta_H, 2, 13)),
    ]
    cxxc_rows = [
        ("random ORF",              fr_c3,    *ci(tab_rand_C, 3, 21)),
        ("proteome (len-matched)",  fmatch_c3, *ci(tab_mt_C, 3, 22)),
        ("de novo neORF",           fdn_c3,   *ci(tab_dn_C, 3, 23)),
        ("metal-binders",           fmt_c3,   *ci(tab_meta_C, 3, 24)),
    ]
    with open(os.path.join(OUT, "fig2_hxh.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["order", "series", "fold", "lo", "hi"])
        for i, r in enumerate(hxh_rows, 1):
            w.writerow([i, r[0], r[1], r[2], r[3]])
    with open(os.path.join(OUT, "fig2_cxxc.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["order", "series", "fold", "lo", "hi"])
        for i, r in enumerate(cxxc_rows, 1):
            w.writerow([i, r[0], r[1], r[2], r[3]])

    # brackets: HxH de novo vs intergenic (from fig1); CxxC de novo vs metal (bootstrap)
    p_c_dn_mt, _, _ = H.boot_diff_p(tab_dn_C, tab_meta_C, 3, nboot=3000, seed=SEED + 30)
    with open(os.path.join(OUT, "fig1_stats.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["panel", "comp", "p", "star"])
        w.writerow(["HxH", "denovo_vs_intergenic", p_dn_ig, star_ig])
        w.writerow(["CxxC", "denovo_vs_metal", p_c_dn_mt, H.stars(p_c_dn_mt)])
    print(f"Fig2 brackets: HxH dn-vs-intergenic p={p_dn_ig:.2e} ({star_ig}); "
          f"CxxC dn-vs-metal p={p_c_dn_mt:.2e} ({H.stars(p_c_dn_mt)})")
    print("wrote fig2_hxh.csv, fig2_cxxc.csv, fig1_stats.csv")
    print("DONE.")


if __name__ == "__main__":
    main()
