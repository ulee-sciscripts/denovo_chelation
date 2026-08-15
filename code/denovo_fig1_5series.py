#!/usr/bin/env python3
"""
Build 5-series HxH and CxxC order-fold tables for the COMBINED figure's panels
B and C, in the order: random ORF, proteome (length-matched), intergenic
(repeat-masked), de novo neORF, metal-binders.

Most folds already exist; the two new ones are proteome-length-matched HxH
(nearest-length whole proteins, d=2) and intergenic-clean CxxC (d=3, from a fresh
Cys-spacing scan of the repeat-masked intergenic ORFs).  Writes, next to the
other figure data (kept separate from fig2_hxh.csv / fig2_cxxc.csv so the
standalone Fig 2 is untouched):
    fig1_hxh5.csv   fig1_cxxc5.csv     (cols: order, series, fold, lo, hi)
"""
import sys, glob, csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hs_common import SCRATCH, load_denovo, load_proteins, read_fasta, count_spacings
import denovo_hisspacing_analysis as H

OUT = os.path.join(HERE, os.pardir, "data", "hisspacing")
CHEL = os.path.join(HERE, os.pardir, "data", "intermediate")
SEED = H.SEED
ORD_C = ord("C")


def rows_for_res(prots, res, min_aa=10):
    rows = []
    for p in prots:
        L = len(p); c = p.count(res)
        if c >= 2 and L >= min_aa:
            o = count_spacings(p, 6, res)
            rows.append((L, c, *[int(x) for x in o]))
    return np.array(rows, dtype=np.int64) if rows else np.zeros((0, 8), np.int64)


def tab_from(L, cnt, o_d, d):
    t = np.zeros((len(L), 8), np.int64); t[:, 0] = L; t[:, 1] = cnt; t[:, 1 + d] = o_d
    return t


def fc(tab, d, seed):
    return H.pooled_fold(tab, d)[0], *H.boot_ci(tab, d, seed=SEED + seed)


# ---------- clean-intergenic Cys (C-x-x-C, d=3) scan ----------
def scan_frame_cys(aa, rows):
    st = aa == H.ORD_STAR
    idx = np.flatnonzero(st)
    b = np.empty(idx.size + 2, np.int64); b[0] = -1; b[-1] = aa.size; b[1:-1] = idx
    for k in range(b.size - 1):
        s = b[k] + 1; e = b[k + 1]; L = e - s
        if L < 10:
            continue
        seg = aa[s:e]; cm = seg == ORD_C; c = int(cm.sum())
        if c < 2:
            continue
        cxxc = int(np.count_nonzero(cm[:-3] & cm[3:])) if L > 3 else 0
        rows.append((L, c, 0, 0, cxxc, 0, 0, 0))


def scan_region_cys(dna, rows):
    x = H.BASE2INT[dna]; comp = np.where(x == 4, 4, 3 - x).astype(np.uint8); xr = comp[::-1]
    for src in (x, xr):
        n = src.size
        for frame in range(3):
            m = (n - frame) // 3
            if m < 10:
                continue
            cod = src[frame:frame + 3 * m].reshape(m, 3).astype(np.int64)
            bad = (cod == 4).any(axis=1); ci = cod[:, 0] * 16 + cod[:, 1] * 4 + cod[:, 2]; ci[bad] = 0
            aa = H.AA_LUT[ci].copy(); aa[bad] = H.ORD_STAR
            scan_frame_cys(aa, rows)


def intergenic_cxxc():
    genes, reps = H.load_intervals()
    rows = []
    for c in H.CHRS:
        seq = None
        for _, s in read_fasta(f"{SCRATCH}/dna_{c}.fa.gz"):
            seq = s.upper(); break
        Lc = len(seq); seqb = np.frombuffer(seq.encode(), np.uint8)
        gmask = H.mask_from_intervals(genes[c], Lc); rmask = H.mask_from_intervals(reps[c], Lc)
        clean = (~gmask) & (~rmask)
        for s, e in H.runs_of_true(clean):
            if e - s >= 30:
                scan_region_cys(seqb[s:e], rows)
    tab = np.array(rows, dtype=np.int64)
    return H.pooled_fold(tab, 3)[0], *H.boot_ci(tab, 3, seed=SEED + 55), tab.shape[0]


def read_fig1_clean_d2():
    with open(os.path.join(OUT, "fig1_profile.csv")) as fh:
        for r in csv.DictReader(fh):
            if int(r["d"]) == 2:
                return float(r["clean"]), float(r["cl_lo"]), float(r["cl_hi"])


def main():
    paths = sorted(glob.glob(f"{SCRATCH}/gc2023/*_DeNovoGenesNuc.fa"))
    prots_dn = [d["prot"] for d in load_denovo(paths, min_aa=10)]
    prot = load_proteins(f"{SCRATCH}/dmel_pep.fa.gz", min_aa=10)
    meta = load_proteins(f"{SCRATCH}/dmel_metalbinding.fasta", min_aa=10)

    # whole-protein nearest-length matched proteome (as in Figs 2/3/5)
    dn_lens = np.array([len(p) for p in prots_dn])
    plen = np.array([len(s) for s in prot]); order = np.argsort(plen, kind="stable")
    psort = [prot[i] for i in order]; plsort = plen[order]; nP = len(psort)
    rngm = np.random.default_rng(SEED); matched = []
    for L in dn_lens:
        j = int(np.searchsorted(plsort, int(L)))
        cand = [k for k in (j - 1, j, j + 1) if 0 <= k < nP]
        diffs = np.array([abs(int(plsort[k]) - int(L)) for k in cand])
        best = [cand[t] for t in range(len(cand)) if diffs[t] == diffs.min()]
        matched.append(psort[int(rngm.choice(best))])

    tdnH = rows_for_res(prots_dn, "H"); tdnC = rows_for_res(prots_dn, "C")
    tmeH = rows_for_res(meta, "H");     tmeC = rows_for_res(meta, "C")
    tmtH = rows_for_res(matched, "H");  tmtC = rows_for_res(matched, "C")

    Lr, Hr, Cr, hxhr, cxxcr = ([], [], [], [], [])
    with open(f"{CHEL}/chelation_percluster_randorf.csv") as fh:
        for r in csv.DictReader(fh):
            Lr.append(int(r["length"])); Hr.append(int(r["His_count"])); Cr.append(int(r["Cys_count"]))
            hxhr.append(int(r["HxH"])); cxxcr.append(int(r["CxxC"]))
    trH = tab_from(np.array(Lr), np.array(Hr), np.array(hxhr), 2)
    trC = tab_from(np.array(Lr), np.array(Cr), np.array(cxxcr), 3)

    igH = read_fig1_clean_d2()
    print("scanning intergenic Cys (C-x-x-C) ...", flush=True)
    igC_f, igC_lo, igC_hi, igC_n = intergenic_cxxc()
    print(f"ANCHOR intergenic CxxC fold={igC_f:.4f} CI[{igC_lo:.4f},{igC_hi:.4f}] n={igC_n}")

    hxh = [("random ORF",) + fc(trH, 2, 1),
           ("proteome (len-matched)",) + fc(tmtH, 2, 2),
           ("intergenic (rpt-masked)", *igH),
           ("de novo neORF",) + fc(tdnH, 2, 3),
           ("metal-binders",) + fc(tmeH, 2, 4)]
    cxxc = [("random ORF",) + fc(trC, 3, 11),
            ("proteome (len-matched)",) + fc(tmtC, 3, 12),
            ("intergenic (rpt-masked)", igC_f, igC_lo, igC_hi),
            ("de novo neORF",) + fc(tdnC, 3, 13),
            ("metal-binders",) + fc(tmeC, 3, 14)]

    for name, rowset in (("fig1_hxh5.csv", hxh), ("fig1_cxxc5.csv", cxxc)):
        with open(os.path.join(OUT, name), "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["order", "series", "fold", "lo", "hi"])
            for i, r in enumerate(rowset, 1):
                w.writerow([i, r[0], r[1], r[2], r[3]])
        print(f"wrote {name}: " + "  ".join(f"{r[0].split()[0]}={r[1]:.3f}" for r in rowset))


if __name__ == "__main__":
    main()
