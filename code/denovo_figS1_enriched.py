#!/usr/bin/env python3
"""
Figure S1 (ENRICHED variant): the ELM/SLiM classes most ENRICHED in neORFs relative
to established proteins (neORF > established, BH-FDR q<0.05), the mirror of the
depleted version.  Reads the pre-registered scan (slim_chelation_results.json) --
no re-scan -- selects the top 15 by z (most significant enrichment), and writes
../data/hisspacing/figS1_slim.csv (same schema the R script already uses).
"""
import os, sys, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from hs_common import SCRATCH
OUT = os.path.join(HERE, os.pardir, "data", "hisspacing")
INTER = os.path.join(HERE, os.pardir, "data", "intermediate")
import json


def bh_q(p):
    p = np.asarray(p, float); m = p.size; o = np.argsort(p); q = np.empty(m); prev = 1.0
    for rank, i in enumerate(o[::-1]):
        r = m - rank; q[i] = min(prev, p[i] * m / r); prev = q[i]
    return q


def star(q):
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 5e-2 else "ns"


def main():
    # NOTE: slim_chelation_results.json is the pre-registered ELM/SLiM regex scan
    # produced by the separate slim_chelation_scan.py pipeline (not part of this
    # figure release). It is a DERIVED intermediate, not a raw download, so it is
    # not fetched by fetch_data.sh. The shipped figure CSV (figS1_slim.csv) already
    # contains the selected 15 classes; to re-derive it place the scan output at
    # data/intermediate/slim_chelation_results.json.
    # TODO: ship or regenerate slim_chelation_results.json to run this driver from raw.
    scan = json.load(open(os.path.join(INTER, "slim_chelation_results.json")))
    S = scan["grand"]; n = scan["n_elm"]
    names = np.array(S["reg_names"][:n])
    dn = np.array(S["reg_obs_density"][:n])
    est = np.array(S["reg_density_established"]["null_mean"][:n])
    gc = np.array(S["reg_density_gc"]["null_mean"][:n])
    z = np.array(S["reg_density_established"]["z"][:n])
    fold = np.array(S["reg_density_established"]["fold"][:n])
    q = bh_q(np.array(S["reg_density_established"]["p_norm"][:n]))

    sel = [i for i in np.argsort(-z) if z[i] > 0 and q[i] < 0.05][:15]   # most-significant enrichments
    with open(f"{OUT}/figS1_slim.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["class", "denovo", "established", "randorf", "fold_dn_over_est", "q", "star"])
        for i in sel:
            w.writerow([str(names[i]), float(dn[i]), float(est[i]), float(gc[i]),
                        float(fold[i]), float(q[i]), star(q[i])])
    # H-x-H motif density (regex hits per 100 aa), comparable to the ELM densities,
    # written separately so the R script can pin it as a top reference row.
    CHEL = INTER

    def hxh_density(fname):
        tot = ln = 0
        for r in csv.DictReader(open(f"{CHEL}/{fname}")):
            tot += int(r["HxH"]); ln += int(r["length"])
        return 100.0 * tot / ln
    hd = {"denovo": hxh_density("chelation_percluster_denovo.csv"),
          "established": hxh_density("chelation_percluster_proteome.csv"),
          "randorf": hxh_density("chelation_percluster_randorf.csv")}
    with open(f"{OUT}/figS1_hxh.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["motif", "denovo", "established", "randorf"])
        w.writerow(["H-x-H", hd["denovo"], hd["established"], hd["randorf"]])
    print(f"wrote figS1_hxh.csv: H-x-H density neORF {hd['denovo']:.3f} est {hd['established']:.3f} rand {hd['randorf']:.3f}")

    n_rand_hi = sum(1 for i in sel if gc[i] >= dn[i])
    print(f"wrote figS1_slim.csv: 15 most-enriched ELM classes (all q<0.05, neORF > established)")
    print(f"  of these, random ORF >= neORF in {n_rand_hi}/15 -> enrichment often not de-novo-specific")
    print(f"  aggregate density: neORF {dn.sum():.1f}  established {est.sum():.1f}  random {gc.sum():.1f}")


if __name__ == "__main__":
    main()
