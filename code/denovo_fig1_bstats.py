#!/usr/bin/env python3
"""
Panel-B (combined figure) comparison p-values for the H-x-H spectrum:
  random ORF vs intergenic (repeat-masked)   -- is intergenic order real vs the null?
  de novo vs metal-binders                    -- the (length-confounded) bulk near-equality
Both are two-sided sequence bootstraps of the difference of pooled d=2 order-folds
(boot_diff_p, 3000 resamples).  Writes ../data/hisspacing/fig1_hxh5_stats.csv.
"""
import os, sys, glob, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from hs_common import SCRATCH, load_denovo, load_proteins, count_spacings
import denovo_hisspacing_analysis as H

OUT = os.path.join(HERE, os.pardir, "data", "hisspacing")
CHEL = os.path.join(HERE, os.pardir, "data", "intermediate"); SEED = H.SEED


def rows_for_res(prots, res, min_aa=10):
    rows = []
    for p in prots:
        L = len(p); c = p.count(res)
        if c >= 2 and L >= min_aa:
            o = count_spacings(p, 6, res); rows.append((L, c, *[int(x) for x in o]))
    return np.array(rows, dtype=np.int64) if rows else np.zeros((0, 8), np.int64)


def tab_from(L, cnt, o_d, d):
    t = np.zeros((len(L), 8), np.int64); t[:, 0] = L; t[:, 1] = cnt; t[:, 1 + d] = o_d
    return t


def main():
    prots_dn = [d["prot"] for d in load_denovo(sorted(glob.glob(f"{SCRATCH}/gc2023/*_DeNovoGenesNuc.fa")), min_aa=10)]
    meta = load_proteins(f"{SCRATCH}/dmel_metalbinding.fasta", min_aa=10)
    tdnH = rows_for_res(prots_dn, "H"); tmeH = rows_for_res(meta, "H")
    Lr, Hr, hxhr = [], [], []
    with open(f"{CHEL}/chelation_percluster_randorf.csv") as fh:
        for r in csv.DictReader(fh):
            Lr.append(int(r["length"])); Hr.append(int(r["His_count"])); hxhr.append(int(r["HxH"]))
    trH = tab_from(np.array(Lr), np.array(Hr), np.array(hxhr), 2)
    print("scanning intergenic His ...", flush=True)
    clean, _ = H.build_intergenic(); tclH = clean.table()
    print("intergenic HxH fold (check ~1.450):", round(H.pooled_fold(tclH, 2)[0], 4))
    p_ri, _, _ = H.boot_diff_p(trH, tclH, 2, nboot=3000, seed=SEED + 61)
    p_dm, _, _ = H.boot_diff_p(tdnH, tmeH, 2, nboot=3000, seed=SEED + 62)
    print(f"random vs intergenic HxH: p={p_ri:.2e} {H.stars(p_ri)}")
    print(f"de novo vs metal HxH:     p={p_dm:.2e} {H.stars(p_dm)}")
    with open(f"{OUT}/fig1_hxh5_stats.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["comp", "p", "star"])
        w.writerow(["random_vs_intergenic", p_ri, H.stars(p_ri)])
        w.writerow(["denovo_vs_metal", p_dm, H.stars(p_dm)])
    print("wrote fig1_hxh5_stats.csv")


if __name__ == "__main__":
    main()
