#!/usr/bin/env python3
"""
Recompute Fig 3 (incidence + per-residue density, de novo vs length-matched) using
the NEAREST-LENGTH WHOLE-PROTEIN comparator (FIGURE_WORK.md §5), so Fig 3 uses the
same length-match method as Figs 2 and 5.  Overwrites, with identical schema:
  ../data/hisspacing/fig1_incidence.csv  (tag,denovo,matched,matched_sd,metal,z,p,star)
  ../data/hisspacing/fig1_density.csv     (tag,dn,dn_lo,dn_hi,mt,mt_lo,mt_hi,mwu_p,star)

Metal-binder incidence (the reference "ceiling") is the bulk value from the
chelation spectrum, unchanged.  Anchors printed: de novo HxH incidence ~0.097 /
density ~0.194 must reproduce the committed spectrum values.
"""
import os, sys, csv, glob, json
from math import erfc, sqrt
import numpy as np
from scipy import stats as st

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hs_common import SCRATCH, load_denovo, load_proteins

OUT = os.path.join(HERE, os.pardir, "data", "hisspacing")
CHEL = os.path.join(HERE, os.pardir, "data", "intermediate", "chelation_derived.json")
SEED = 20260812


def per_protein_pairs(seqs, res, d):
    """(obs, L) per protein: obs = ordered same-residue pairs at spacing d."""
    rb = ord(res)
    obs = np.empty(len(seqs), np.int64); L = np.empty(len(seqs), np.int64)
    for i, s in enumerate(seqs):
        n = len(s)
        m = np.frombuffer(s.encode(), np.uint8) == rb
        obs[i] = np.count_nonzero(m[:n - d] & m[d:]) if n > d else 0
        L[i] = n
    return obs, L


def boot_density(obs, L, seed, nb=3000):
    rng = np.random.default_rng(seed)
    obs = obs.astype(float); L = L.astype(float); n = len(obs)
    pt = 100 * obs.sum() / L.sum()
    b = np.empty(nb)
    for i in range(nb):
        ix = rng.integers(0, n, n); b[i] = 100 * obs[ix].sum() / L[ix].sum()
    return pt, float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def boot_incidence_sd(obs, seed, nb=3000):
    rng = np.random.default_rng(seed)
    ind = (obs >= 1).astype(float); n = len(ind)
    b = np.empty(nb)
    for i in range(nb):
        b[i] = ind[rng.integers(0, n, n)].mean()
    return float(np.std(b))


def two_prop_z(x1, n1, x2, n2):
    p1, p2 = x1 / n1, x2 / n2
    p = (x1 + x2) / (n1 + n2)
    se = sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    z = (p1 - p2) / se
    return z, erfc(abs(z) / sqrt(2))


def star(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"


def main():
    den = load_denovo(sorted(glob.glob(f"{SCRATCH}/gc2023/*_DeNovoGenesNuc.fa")), min_aa=10)
    dn_seqs = [d["prot"] for d in den]
    prot = load_proteins(f"{SCRATCH}/dmel_pep.fa.gz", min_aa=10)
    meta = load_proteins(f"{SCRATCH}/dmel_metalbinding.fasta", min_aa=10)
    print(f"de novo={len(dn_seqs)}  proteome={len(prot)}  metal={len(meta)}")

    # nearest-length WHOLE-protein matched set (with replacement), SEED tie-break
    dn_lens = np.array([len(s) for s in dn_seqs])
    plen = np.array([len(s) for s in prot]); order = np.argsort(plen, kind="stable")
    psort = [prot[i] for i in order]; plsort = plen[order]; nP = len(psort)
    rngm = np.random.default_rng(SEED)
    matched = []
    for L in dn_lens:
        j = int(np.searchsorted(plsort, int(L)))
        cand = [k for k in (j - 1, j, j + 1) if 0 <= k < nP]
        diffs = np.array([abs(int(plsort[k]) - int(L)) for k in cand])
        best = [cand[t] for t in range(len(cand)) if diffs[t] == diffs.min()]
        matched.append(psort[int(rngm.choice(best))])

    R = json.load(open(CHEL)); sp = R["spectrum"]["metalbinding"]
    tags = [("HxH", "H", 2), ("CxxC", "C", 3)]

    inc_rows, den_rows = [], []
    for k, (tag, res, d) in enumerate(tags):
        o_dn, L_dn = per_protein_pairs(dn_seqs, res, d)
        o_mt, L_mt = per_protein_pairs(matched, res, d)
        dn_inc = float((o_dn >= 1).mean()); mt_inc = float((o_mt >= 1).mean())
        mt_sd = boot_incidence_sd(o_mt, SEED + 40 + k)
        metal_inc = sp[f"{tag}_incidence"]
        z, p_inc = two_prop_z(int((o_dn >= 1).sum()), len(o_dn),
                              int((o_mt >= 1).sum()), len(o_mt))
        inc_rows.append([tag, dn_inc, mt_inc, mt_sd, metal_inc, z, p_inc, star(p_inc)])

        p1, l1, h1 = boot_density(o_dn, L_dn, SEED + 10 + k)
        p2, l2, h2 = boot_density(o_mt, L_mt, SEED + 20 + k)
        o_me, L_me = per_protein_pairs(meta, res, d)             # metal-binders (bulk reference)
        p3, l3, h3 = boot_density(o_me, L_me, SEED + 30 + k)
        dens_dn = 100.0 * o_dn / L_dn; dens_mt = 100.0 * o_mt / L_mt
        mwu_p = float(st.mannwhitneyu(dens_dn, dens_mt, alternative="two-sided").pvalue)
        den_rows.append([tag, p1, l1, h1, p2, l2, h2, mwu_p, star(mwu_p), p3, l3, h3])
        print(f"{tag}: incid dn={dn_inc:.4f} matched={mt_inc:.4f} metal={metal_inc:.4f} "
              f"z={z:.2f} p={p_inc:.1e} | dens dn={p1:.3f}[{l1:.3f},{h1:.3f}] "
              f"matched={p2:.3f}[{l2:.3f},{h2:.3f}] MWU p={mwu_p:.1e}")

    with open(f"{OUT}/fig1_incidence.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["tag", "denovo", "matched", "matched_sd", "metal", "z", "p", "star"])
        w.writerows(inc_rows)
    with open(f"{OUT}/fig1_density.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["tag", "dn", "dn_lo", "dn_hi", "mt", "mt_lo", "mt_hi",
                                        "mwu_p", "star", "me", "me_lo", "me_hi"])
        w.writerows(den_rows)
    print("ANCHOR de novo HxH incidence 0.0972 / density 0.194 (spectrum) should match above.")
    print("wrote fig1_incidence.csv, fig1_density.csv (whole-protein comparator)")


if __name__ == "__main__":
    main()
