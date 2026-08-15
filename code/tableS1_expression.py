#!/usr/bin/env python3
"""Regenerate Table S1 = per-replicate expression of the 267 mappability-clean non-overlapping
neORFs (correct dm6 coords, H-x-H content, ZMEG labeled). Machine-readable CSV."""
import os, csv, re
base=os.environ.get("ZROOT","../../zinc_analysis")   # analysis data tree
out=os.environ.get("ZS1OUT","../data/TableS1_neORF_expression_267.csv")

# --- counts.txt (featureCounts): geneid -> (len, counts[]); per-sample RPK sum for TPM ---
with open(base+"/star_counts/counts.txt") as f:
    f.readline(); hdr=f.readline().rstrip("\n").split("\t")
    srr=[re.search(r"SRR\d+",c).group(0) for c in hdr[6:]]; n=len(srr)
    rpk=[0.0]*n; rows={}
    for line in f:
        p=line.rstrip("\n").split("\t"); L=int(p[5]); c=list(map(int,p[6:]))
        for i in range(n): rpk[i]+=c[i]/(L/1000.0)
        rows[p[0]]=(L,c)
def tpm_vec(gid):
    L,c=rows[gid]; return [(c[i]/(L/1000.0))/rpk[i]*1e6 for i in range(n)]

# --- sample sheet: run -> (genotype, condition); build ordered 18-sample layout ---
smap={r["run"]:(r["genotype"],r["condition"]) for r in csv.DictReader(open(base+"/sample_sheet.tsv"),delimiter="\t")}
GENOS=["WT","ZnT63-KO","IA2-KO"]; CONDS=["control","Zn","Mn"]
sample_cols=[]   # (label, sample_index)
for g in GENOS:
    for cond in CONDS:
        idxs=[i for i,s in enumerate(srr) if smap.get(s)==(g,cond)]
        for rep,i in enumerate(sorted(idxs),1):
            sample_cols.append((f"TPM_{g}_{cond}_rep{rep}", i))
cond_idx={cond:[i for i,s in enumerate(srr) if smap.get(s) and smap[s][1]==cond] for cond in CONDS}

# --- symbol: only ZMEG is a named candidate in the paper ---
sym={"UA_STRG.15519.1_1_193_357":"ZMEG"}

# --- 267 mappable set ---
mset=[r for r in csv.DictReader(open(base+"/neORF_clean_mappable.tsv"),delimiter="\t")
      if r.get("mappability_flag","pass")=="pass"]
recs=[]
missing=0
for r in mset:
    nm=r["name"]
    if nm not in rows: missing+=1; continue
    tv=tpm_vec(nm)
    means={cond: (sum(tv[i] for i in cond_idx[cond])/len(cond_idx[cond])) for cond in CONDS}
    rec=dict(neORF=nm, symbol=sym.get(nm,""), chrom=r["chrom"], start_dm6=r["start_dm6"], end_dm6=r["end_dm6"],
             strand=r["strand"], len_aa=r["len_aa"], n_hxh=r["n_hxh"], has_hxh=r["has_hxh"], poly_hxh=r["poly_hxh"],
             mean_control_TPM=round(means["control"],3), mean_Zn_TPM=round(means["Zn"],3), mean_Mn_TPM=round(means["Mn"],3))
    for lab,i in sample_cols: rec[lab]=round(tv[i],3)
    recs.append(rec)
recs.sort(key=lambda r: max(r["mean_control_TPM"],r["mean_Zn_TPM"],r["mean_Mn_TPM"]), reverse=True)

cols=["neORF","symbol","chrom","start_dm6","end_dm6","strand","len_aa","n_hxh","has_hxh","poly_hxh",
      "mean_control_TPM","mean_Zn_TPM","mean_Mn_TPM"]+[c for c,_ in sample_cols]
with open(out,"w",newline="") as f:
    wtr=csv.DictWriter(f,fieldnames=cols); wtr.writeheader()
    for r in recs: wtr.writerow(r)

# --- verify ---
print("wrote", out)
print("rows:", len(recs), "| missing-from-counts:", missing)
print("has_hxh==1:", sum(1 for r in recs if r["has_hxh"]=="1"), "| poly_hxh==1:", sum(1 for r in recs if r["poly_hxh"]=="1"))
for nm in ("UA_STRG.15519.1_1_193_357","ZI_STRG.9622.1_3_649_837"):
    rr=[r for r in recs if r["neORF"]==nm]
    if rr: r=rr[0]; print(f"  {r['symbol']}: {r['chrom']}:{r['start_dm6']}-{r['end_dm6']} ctrl={r['mean_control_TPM']} Zn={r['mean_Zn_TPM']} Mn={r['mean_Mn_TPM']} n_hxh={r['n_hxh']} poly={r['poly_hxh']}")
    else: print(f"  {nm}: MISSING")
print("top-3 expressers:", [(r["neORF"][:22],max(r['mean_control_TPM'],r['mean_Zn_TPM'],r['mean_Mn_TPM']),"hxh"+r["has_hxh"]) for r in recs[:3]])
