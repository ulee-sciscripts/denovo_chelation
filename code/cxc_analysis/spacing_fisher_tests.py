#!/usr/bin/env python3
"""GBE revision — Fisher's exact tests for bis-histidine SPACING incidence
(fraction of de novo ORFs carrying >=1 motif at spacing d), against two
composition-aware nulls, as an exact-test complement to the permutation
order-fold (data/hisspacing/fig1_profile.csv). R2 asked for exact tests.

Nulls:
  (1) composition-matched random-ORF null (randorf): i.i.d. ORF-like peptides at
      intergenic GC (same n, mean His, mean length as the de novo set). Per-cluster
      motif tables (translate-to-stop).
  (2) dinucleotide-preserving (Altschul-Erikson) shuffle of each de novo CDS,
      translate_full (positions fixed) -- holds EXACT dinucleotide (hence CA)
      content per sequence. R replicate matched null sets; report median OR/p.

Spacings: H-H (d1), H-x-H (d2), H-x-x-H (d3, inert), H-x-x-x-H (d4).
Std-lib + numpy + scipy; deterministic (SEED).
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import os, sys, csv
import numpy as np
from scipy import stats

CODE = os.path.join(_HERE, "..")
INTER = os.path.join(_HERE, "..", "..", "data", "intermediate")
NEORF = f"{_HERE}/../../data/analysis/neORF.fa"
OUT_TSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "spacing_fisher_results.tsv")
sys.path.insert(0, CODE)
from hs_common import read_fasta, translate, translate_full          # noqa: E402
from dinuc_shuffle import altschul_erikson_shuffle                    # noqa: E402

SEED = 20260812
R = 10                                    # replicate dinuc-shuffle null sets
SPAC = [("HH", "H-H", 1), ("HxH", "H-x-H", 2),
        ("HxxH", "H-x-x-H", 3), ("HxxxH", "H-x-x-x-H", 4)]
rows_out = []


def fisher(dw, dtot, cw, ctot):
    OR, p = stats.fisher_exact([[dw, dtot - dw], [cw, ctot - cw]])
    return OR, p


# ---------- (1) random-ORF composition-matched null (per-cluster CSVs) ----------
def load_csv(f):
    return list(csv.DictReader(open(f)))


dn = load_csv(f"{INTER}/chelation_percluster_denovo.csv")
ro = load_csv(f"{INTER}/chelation_percluster_randorf.csv")
mHis = lambda R_: sum(int(r["His_count"]) for r in R_) / len(R_)
mLen = lambda R_: sum(int(r["length"]) for r in R_) / len(R_)
print("=== (1) Fisher's exact: de novo vs composition-matched random-ORF null ===")
print(f"    de novo n={len(dn)} (His {mHis(dn):.2f}, len {mLen(dn):.0f}) | "
      f"randorf n={len(ro)} (His {mHis(ro):.2f}, len {mLen(ro):.0f})")
for col, name, d in SPAC:
    dw = sum(1 for r in dn if int(r[col]) > 0)
    cw = sum(1 for r in ro if int(r[col]) > 0)
    OR, p = fisher(dw, len(dn), cw, len(ro))
    tag = "ns" if p >= .05 else "*"
    print(f"    {name:10s} dn {100*dw/len(dn):5.2f}%  rand {100*cw/len(ro):5.2f}%  "
          f"OR={OR:5.2f}  p={p:.2e}  {tag}")
    rows_out.append(["randorf", name, f"{100*dw/len(dn):.2f}", f"{100*cw/len(ro):.2f}",
                     f"{OR:.3f}", f"{p:.3e}", tag])

# ---------- load de novo CDS for the dinucleotide-shuffle null ----------
seen, cds_list = set(), []
for hdr, seq in read_fasta(NEORF):
    seq = seq.upper()
    prot, cds = translate(seq, truncate_at_stop=True)
    if len(prot) < 10 or prot in seen:
        continue
    seen.add(prot)
    cds_list.append(cds)
valid = [c for c in cds_list if len(c) >= 6 and set(c) <= set("ACGT")]
Nv = len(valid)


def spac_inc_vec(pep):
    a = np.frombuffer(pep.encode(), np.uint8)
    m = a == ord("H")
    L = a.size
    return [int(L > d and bool(np.any(m[:-d] & m[d:]))) for (_, _, d) in SPAC]


# observed incidence (translate_full, positions fixed -- matches the null's method)
obs = np.zeros(len(SPAC), int)
for c in valid:
    for j, v in enumerate(spac_inc_vec(translate_full(c))):
        obs[j] += v

# ---------- (2) dinucleotide-preserving (AE) shuffle null, R replicate sets ----------
print(f"\n=== (2) Fisher's exact: de novo vs dinucleotide-preserving (AE) shuffle "
      f"(R={R} matched null sets, N={Nv}) ===")
rng = np.random.default_rng(SEED + 11)
null_counts = np.zeros((R, len(SPAC)), int)
for r in range(R):
    for c in valid:
        pep = translate_full(altschul_erikson_shuffle(c, rng))
        for j, v in enumerate(spac_inc_vec(pep)):
            null_counts[r, j] += v
    print(f"    replicate {r+1}/{R} done")
for j, (col, name, d) in enumerate(SPAC):
    ORs = np.empty(R); ps = np.empty(R)
    for r in range(R):
        ORs[r], ps[r] = fisher(int(obs[j]), Nv, int(null_counts[r, j]), Nv)
    obs_pct = 100 * obs[j] / Nv
    null_pct = 100 * null_counts[:, j].mean() / Nv
    tag = "ns" if np.median(ps) >= .05 else "*"
    print(f"    {name:10s} obs {obs_pct:5.2f}%  dinuc-null {null_pct:5.2f}%  "
          f"OR med={np.median(ORs):5.2f}  p med={np.median(ps):.2e} "
          f"[{ps.min():.1e}-{ps.max():.1e}]  {tag}")
    rows_out.append(["dinuc_AE", name, f"{obs_pct:.2f}", f"{null_pct:.2f}",
                     f"{np.median(ORs):.3f}", f"{np.median(ps):.3e}", tag])

with open(OUT_TSV, "w") as fh:
    fh.write("null\tspacing\tobs_pct\tnull_pct\tOR\tp\tsig\n")
    for r in rows_out:
        fh.write("\t".join(r) + "\n")
print(f"\nwrote {OUT_TSV}")
